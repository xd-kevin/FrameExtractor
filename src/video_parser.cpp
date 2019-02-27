//
// Created by shawntang(唐兴) on 2019-02-15.
//

#include "extractor/video_parser.h"

using namespace extractor;
using namespace std;

VideoParser::VideoParser() {
    _nfrm_total = 0;
    _nfrm_given = 0;

    _v_shot_ranges.clear();
}

void VideoParser::release_memory() {
    _v_frm_rgb.clear();
    _v_frm_gray.clear();
}

int VideoParser::read_video(const string &in_video, int step_sz) {

    VideoCapture vr(in_video);

    if(!vr.isOpened()){
        return -1;
    }

    _nfrm_total = vr.get(CV_CAP_PROP_FRAME_COUNT);
    _video_width = vr.get(CV_CAP_PROP_FRAME_WIDTH);
    _video_height = vr.get(CV_CAP_PROP_FRAME_HEIGHT);
    _video_fps = vr.get(CV_CAP_PROP_FPS);
    if(_video_fps != _video_fps)
        _video_fps = 29.97;

    _step_sz = step_sz;

    //Notice: delete some boundary conditions

    _nfrm_total = 0;
    while (true)
    {
        Mat frm;
        vr >> frm;
        if(frm.empty()) break;

        if(_nfrm_total % _step_sz == 0)
            _v_frm_rgb.push_back(frm);
        _nfrm_total++;
    }
    vr.release();
    _nfrm_given = (int) _v_frm_rgb.size();
    _video_sec = (double) _nfrm_given / _video_fps;

    _v_frm_gray.assign(_nfrm_given,Mat());

#pragma omp parallel for
    for( int i = 0; i < _nfrm_given; i++){
        Mat frm_gray;
        cvtColor(_v_frm_rgb[i],frm_gray,CV_BGR2GRAY);
        GaussianBlur(frm_gray,frm_gray,Size(3,3), 0, 0);
        frm_gray.copyTo(_v_frm_gray[i]);
    }

    _v_frm_valid.assign(_nfrm_given, true);

    _X_diff = Mat(_nfrm_given,1,CV_64F,Scalar(0,0,0)[0]);
    _X_ecr = Mat(_nfrm_given, 1, CV_64F, Scalar(0,0,0)[0]);

    return 0;


}

vector<extractor::ShotRange> VideoParser::parse_video(const string &in_video, extractor::parser_params opt) {
    int ret = read_video(in_video, opt.step_sz);

    if(ret<0){
        fprintf( stderr, "VideoParser: Failed to open input video: %s\n", in_video.c_str());
        return vector<extractor::ShotRange>();
    }

    meta.nframes = _nfrm_total;
    meta.width = _video_width;
    meta.height = _video_height;
    meta.fps = _video_fps;
    meta.duration = _video_sec;

    if(!opt.sample){

        if(opt.fltr_begin_sec > .0 || opt.fltr_end_sec >.0)
            filter_heuristic(opt.fltr_begin_sec, opt.fltr_end_sec);
        cout<<"finish filter_heuristic"<<endl;
        if(opt.fltr_lq)
            filter_low_quality();
        cout<<"finish filter low quality"<<endl;
        filter_transition();
        cout<<"finish filter transition"<<endl;
        extract_histo_features();
        cout<<"finish extract_histo_features"<<endl;
        //break up shots if too long, may affect result.
        double min_shot_len_sec = 2.0;
        post_process(min_shot_len_sec);
        cout<<"finish post_process"<<endl;
        release_memory();

        update_shot_range();
        cout<<"update_shot_range"<<endl;
        if(opt.fltr_rdt)
            filter_redundant_and_obtain_subshots();
        cout<<"finish filter redundant and obtain subshots"<<endl;


        return _v_shot_ranges;
    } else{
        subsample_sec();
        return vector<extractor::ShotRange>();
    }
}

void VideoParser::filter_heuristic(double fltr_begin_sec, double fltr_end_sec) {

    int fltr_begin_nfrms = ceil(fltr_begin_sec * _video_fps / (double)_step_sz);
    int fltr_end_nfrms = ceil(fltr_end_sec * _video_fps / (double) _step_sz);

    for(int i = 0; i < fltr_begin_nfrms; i++)
        mark_invalid(_v_frm_valid, i);
    for(int i = 0; i < fltr_end_nfrms; i++)
        mark_invalid(_v_frm_valid,_nfrm_total-1-i);
}

void VideoParser::filter_low_quality(double thrsh_bright, double thrsh_sharp, double thrsh_uniform) {

    int nfrm_nperc = (int) (0.15*_nfrm_given); // filter at most 15%

    vector<double> v_brightness(_nfrm_given, 0.0);
    vector<double> v_sharpness(_nfrm_given, 0.0);
    vector<double> v_uniformity(_nfrm_given, 0.0);

#pragma omp parallel for
    for(int i = 0; i < _nfrm_given; i++){
        v_brightness[i] = extractor::calc_brightness(_v_frm_rgb[i]);
        v_sharpness[i] = extractor::calc_sharpness(_v_frm_gray[i]);
        v_uniformity[i] = extractor::calc_uniformity(_v_frm_gray[i]);
    }

    vector<size_t> v_srt_idx;
    vector<double> v_srt_val;


    extractor::sort(v_sharpness,v_srt_val,v_srt_idx);
    for(int i =0; i < nfrm_nperc; i++)
        if(v_srt_val[i] <= thrsh_sharp)
            mark_invalid(_v_frm_valid, v_srt_idx[i]);
    extractor::sort(v_brightness,v_srt_val,v_srt_idx);
    for(int i = 0; i < nfrm_nperc; i++)
        if(v_srt_val[i] <= thrsh_bright)
            mark_invalid(_v_frm_valid, v_srt_idx[i]);

    extractor::sort(v_uniformity, v_srt_val, v_srt_idx);
    for(int i = 0; i < nfrm_nperc; i++)
        if(v_srt_val[_nfrm_given-i-1] >= thrsh_uniform)
            mark_invalid(_v_frm_valid,v_srt_idx[_nfrm_given-i-1]);
}

void VideoParser::filter_transition(double thrsh_diff, double thrsh_ecr) {

    int nfrm_nperc = (int) (0.10*_nfrm_given);

    vector<double> v_diff(_nfrm_given, 0.0);
    vector<double> v_ecr(_nfrm_given, 0.0);

    vector<size_t> v_srt_idx;
    vector<double> v_srt_val;

    int img_sz = _v_frm_gray[0].cols * _v_frm_gray[0].rows;
#pragma omp parallel for
    for(int i=1; i < _nfrm_given -1; i++){
        v_diff[i] = (double) (cv::norm((_v_frm_rgb[i] - _v_frm_rgb[i-1])) + cv::norm((_v_frm_rgb[i] - _v_frm_rgb[i+1])))/(2.0*img_sz);
        _X_diff.at<double>(i) = v_diff[i];
    }
    //cout<<"finish first for"<<endl;

    {
        int dl_sz = 5;
        Mat dl_elm = getStructuringElement(MORPH_CROSS, Size(2*dl_sz+1, 2*dl_sz+1), Point(dl_sz,dl_sz));

        vector<Mat> v_edge;
        v_edge.resize(_nfrm_given,Mat());

        vector<cv::Mat> v_edge_dl;
        v_edge_dl.resize(_nfrm_given,Mat());


#pragma omp parallel for
        //cout<<"enter in 2nd for"<<endl;
        for(int i=0; i < _nfrm_given; i++)
        {
            Mat tmp;
            double theta = threshold(_v_frm_gray[i], tmp, 0, 255, CV_THRESH_BINARY|CV_THRESH_OTSU);
            Canny(_v_frm_gray[i], v_edge[i], theta, 1.2*theta);
            dilate(v_edge[i], v_edge_dl[i],dl_elm);
            v_edge[i] -= 254;
            v_edge_dl[i] -= 254;
        }
#pragma omp parallel for
        //cout<<"enter in 3rd for"<<endl;
        for(int i =1; i< _nfrm_given; i++)
        {
            double rho_out,rho_in;
            //cout<<max(1e-6,sum(v_edge[i-1])[0])<<endl;
            rho_out = 1.0 - min(1.0, sum(v_edge[i-1].mul(v_edge_dl[i]))[0]/max(1e-6,sum(v_edge[i-1])[0]));
            rho_in  = 1.0 - min(1.0, sum(v_edge_dl[i-1].mul(v_edge[i]))[0]/max(1e-6,sum(v_edge[i-1])[0]));
            //cout<<"i:"<<i<<"rho_out:"<<rho_out<<"rho_in:"<<rho_in<<endl;
            v_ecr[i] = max(rho_out, rho_in);

            _X_ecr.at<double>(i) = v_ecr[i];
            //cout<<_X_ecr.at<double>(i) <<endl;
        }
        //cout<<"end for"<<endl;
    }

    extractor::sort(v_diff, v_srt_val,v_srt_idx);
    for(int i = 0; i < nfrm_nperc; i++)
        if(v_srt_val[_nfrm_given-i-1] >= thrsh_diff)
            mark_invalid(_v_frm_valid, v_srt_idx[_nfrm_given-i-1]);

    extractor::sort(v_ecr,v_srt_val,v_srt_idx);
    for(int i =0; i < nfrm_nperc; i++)
        if(v_srt_val[_nfrm_given-i-1] >= thrsh_ecr)
            mark_invalid(_v_frm_valid,v_srt_idx[_nfrm_given-i-1]);
}

void VideoParser::extract_histo_features(int pyr_level, bool omit_filtered, int nbin_color, int nbin_edge_ori,
                                         int nbin_edge_mag) {
    int npatches = 0;
    for(int l = 0; l < pyr_level; l++)
        npatches += pow(4,l);

    int nbin_edge = nbin_edge_mag + nbin_edge_ori;

    Mat X_color_hist = Mat(npatches*3*nbin_color, _nfrm_given, CV_32F, Scalar(0,0,0)[0]);
    Mat X_edge_hist = Mat(npatches*nbin_edge, _nfrm_given, CV_32F, Scalar(0,0,0)[0]);

#pragma omp parallel for
    for(int i =0; i < _nfrm_given; i++)
    {
        if(omit_filtered && !_v_frm_valid[i]) continue;

        Mat color_hist;
        extractor::calc_pyr_color_hist(_v_frm_rgb[i], color_hist, nbin_color, pyr_level);
        color_hist.copyTo(X_color_hist.col(i));

        Mat edge_hist;
        extractor::calc_pyr_edge_hist(_v_frm_gray[i], edge_hist, nbin_edge_ori, nbin_edge_mag, pyr_level);
        edge_hist.copyTo(X_edge_hist.col(i));
    }

    X_color_hist = X_color_hist.t();
    X_edge_hist = X_edge_hist.t();

    hconcat(X_color_hist,X_edge_hist,_X_feat);
}

void VideoParser::post_process(double min_shot_sec) {
    int start_idx = -1, end_idx = -1, shotlen = -1;
    int min_shot_len = min_shot_sec * _video_fps / _step_sz;
    int max_shot_len = 3 * min_shot_len;

    for (size_t i = 0; i < _v_frm_valid.size(); ++i) {
        if(start_idx < 0 && _v_frm_valid[i])
            start_idx = i;
        if(start_idx >=0 && (!_v_frm_valid[i] || i+1 == _v_frm_valid.size()))
        {
            end_idx = i;
            shotlen = end_idx - start_idx + 1;
            if(shotlen >= max_shot_len)
            {
                int njumps = floor(shotlen / min_shot_len);
                vector<int> jump;
                vector<double> v_diff;
                for (int i = start_idx; i <= end_idx ; ++i) {
                    v_diff.push_back(_X_diff.at<double>(i));
                }
                sbd_heuristic(v_diff, jump, njumps, min_shot_len);
                for (size_t k = 0; k < jump.size(); ++k) {
                    mark_invalid(_v_frm_valid, start_idx + jump[k]-1);
                    mark_invalid(_v_frm_valid, start_idx + jump[k]);
                }
            }
            start_idx = end_idx = -1;
        }
    }
}

void VideoParser::update_shot_range(int min_shot_len) {
    _v_shot_ranges.clear();

    int sb0 = -1, sb1 = -1;
    for (int i = 0; i < _nfrm_given; ++i) {
        if(_v_frm_valid[i]){
            if(sb0 < 0)
                sb0 = i;
            sb1 = i;
        }

        if(sb0 >= 0 && sb1 >=0 && (!_v_frm_valid[i] || i+1 == _nfrm_given))
        {
            extractor::ShotRange r(sb0, sb1);
            if(r.length() > min_shot_len)
                _v_shot_ranges.push_back(r);
            else{
                for(int j=sb0; j<= sb1; j++)
                    mark_invalid(_v_frm_valid, j);
            }
            sb0=sb1=-1;
        }
    }
}

void VideoParser::filter_redundant_and_obtain_subshots() {


    if(_v_shot_ranges.empty())
        update_shot_range();

    int nfrm_valid = get_nfrm_valid();

    if(nfrm_valid == 0)
        return;

    Mat km_data( nfrm_valid, _X_feat.cols, _X_feat.type());
    vector<int> v_idxmap(nfrm_valid, 0);

    int row = 0;
    for(int i = 0; i<_nfrm_given; i++)
    {
        if(_v_frm_valid[i]) {
            _X_feat.row(i).copyTo(km_data.row(row));
            v_idxmap[row] = i;
            row++;
        }
    }

    int ncluster = min(nfrm_valid/2, (int)_v_shot_ranges.size());
    int per_kfrms = 0;

    Mat km_lbl;
    Mat km_ctr;
    extractor::perform_kmeans(km_data, km_lbl, km_ctr, ncluster);

    vector<int> v_frm_clusterid(_nfrm_given, -1);
    for (int i = 0; i < km_lbl.rows ; ++i)
        v_frm_clusterid[v_idxmap[i]] = km_lbl.at<int>(i);

    int shot_range_num = 0;

    for (size_t shotid = 0; shotid < _v_shot_ranges.size(); ++shotid) {
        int sb0 = _v_shot_ranges[shotid].start;
        int sb1 = _v_shot_ranges[shotid].end;

        int ssb0 = -1, ssb1 = -1, lbl = -1;

        for (int j = sb0; j <= sb1; ++j) {
            if(_v_frm_valid[j]){
                if(ssb0<0){
                    ssb0 = j;
                    lbl = v_frm_clusterid[j];
                }
                ssb1 = j;
            }
            if(ssb0 >= 0 && (v_frm_clusterid[j] != lbl || j == sb1))
            {
                extractor::Range r(ssb0, ssb1);
                _v_shot_ranges[shotid].v_range.push_back(r);
                shot_range_num += 1;
                ssb0 = ssb1 = lbl = -1;
            }
        }
        sb0 = sb1 =-1;
    }


    per_kfrms = ceil( _video_sec / shot_range_num);


    for(size_t shotid = 0; shotid < _v_shot_ranges.size(); ++shotid){
        for(size_t sshotid = 0; sshotid < _v_shot_ranges[shotid].v_range.size();sshotid++){
            extractor::Range r = _v_shot_ranges[shotid].v_range[sshotid];
            vector<double> diff_vec;
            vector<size_t> idx_map;
           // cout<<"shotid:"<<shotid<<"r.start:"<<r.start<<",r.end"<<r.end<<endl;
            for (int i = r.start; i < r.end ; ++i){
                diff_vec.push_back(_X_diff.at<double>(i));
                idx_map.push_back(i);
            }
            vector<double> diff_srt;
            vector<size_t> diff_srt_idx;
            extractor::sort(diff_vec,diff_srt,diff_srt_idx);


            for (int j = 0; j < per_kfrms && j < diff_srt_idx.size(); ++j) {
               // cout<<idx_map[diff_srt_idx[j]]<<endl;
                r.v_idx.push_back(idx_map[diff_srt_idx[j]]);
                _v_shot_ranges[shotid].v_idx.push_back(idx_map[diff_srt_idx[j]]);
            }
        }
    }

}

void VideoParser::sbd_heuristic(vector<double> v_diff, vector<int> &jump, int njumps, int min_shot_len) {

    vector<size_t> v_srt_idx;
    vector<double> v_srt_val;
    extractor::sort(v_diff, v_srt_val, v_srt_idx);
    for (int i = (int)v_srt_val.size() -1; i >= 0; --i) {
        bool add = true;
        if(((unsigned)v_srt_idx[i] + 1 < min_shot_len) ||
           ((unsigned)v_diff.size() - v_srt_idx[i] <min_shot_len))
            add = false;
        else{
            for (size_t j = 0; j < jump.size(); ++j) {
                int len = abs(jump[j] - (int)v_srt_idx[i]) + 1;
                if(len < min_shot_len){
                    add = false;
                    break;
                }
            }
        }
        if(add){
            jump.push_back((int)v_srt_idx[i]);
        }
        if((int) jump.size() == njumps)
            break;
    }
}

void VideoParser::subsample_sec() {

    int n_frm = (int) _video_sec;
    int fps = (int) _video_fps;
    vector<int> valid_idx(n_frm, 0);
    default_random_engine e(time(0));
    unsigned int beg = 0;
    unsigned int end = fps-1;

    for (int i = 0; i < n_frm; ++i) {
        uniform_int_distribution<unsigned> u(beg, end);
        int rnd_idx = u(e);
        while(true){
            if(_v_frm_valid[rnd_idx])
            {
                valid_idx[i] = rnd_idx;
                break;
            }
            else
                rnd_idx = u(e);
        }
        beg += fps;
        end += fps;
    }


    for (int j = 0; j < _nfrm_given; ++j) {
        vector<int>::iterator it = find(valid_idx.begin(), valid_idx.end(),j);
        if(it != valid_idx.end())
            continue;
        mark_invalid(_v_frm_valid,j);
    }
}

void VideoParser::mark_invalid(vector<bool> &vec, int idx, int wnd_sz) {
    int vec_len = (int) vec.size();
    for(int i = max(0,idx-wnd_sz); i<= min(vec_len-1, idx+wnd_sz); i++)
        vec[i] = false;
}

int VideoParser::get_nfrm_valid() {
    return (int) accumulate(_v_frm_valid.begin(),_v_frm_valid.end(),0);
}