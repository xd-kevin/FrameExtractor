#include "extractor/video_parser.h"
#include "extractor/extractor.h"
#include "extractor/file_helper.h"

using namespace std;
using namespace cv;


int main(int argc, char** argv) {

    if(argc <2)
        usage();

    extractor_params opt;
    extractor_parse_params(argc, argv, opt);

    vector<int> v_frm_idx;
    run_extractor(opt, v_frm_idx);

    return 0;
}

void run_extractor(extractor_params& opt, vector<int>& v_frm_idx){
    if(!extractor::file_exsits(opt.in_video)){
        fprintf(stderr, "File not exist: %s\n", opt.in_video.c_str());
        return;
    }
    v_frm_idx.clear();

    extractor::VideoParser parser;
    vector<extractor::ShotRange> v_shot_range;

    extractor::parser_params parser_opt;
    parser_opt.step_sz = opt.step_sz;
    parser_opt.fltr_begin_sec = (opt.fltr_begin_sec<0)?max(0.5, 0.05*parser.meta.duration):opt.fltr_begin_sec;
    parser_opt.fltr_end_sec = (opt.fltr_end_sec<0)?max(0.5,0.10*parser.meta.duration):opt.fltr_end_sec;
    parser_opt.sample = opt.sample;
    cout<<"enter in run_extractor"<<endl;

    v_shot_range = parser.parse_video(opt.in_video, parser_opt);

    char strbuf[256];
    string filename = extractor::get_filename(std::string(opt.in_video));
    if(opt.sample){
        for (int i = 0; i < parser._v_frm_valid.size(); ++i) {
            if(parser._v_frm_valid[i])
                v_frm_idx.push_back(i);
        }
    }
    else{
        for (int shotid = 0; shotid < v_shot_range.size(); ++shotid) {
            for (int i = 0; i < v_shot_range[shotid].v_idx.size(); ++i) {
                v_frm_idx.push_back(v_shot_range[shotid].v_idx[i]);
                //cout<<v_shot_range[shotid].v_idx[i]<<endl;
            }
        }
    }

    VideoCapture vr(opt.in_video);
    int njpg = 0;
    int frm_idx = 0;

    if(opt.sample){
        while(njpg < (int)v_frm_idx.size()){
            Mat frm;
            vr >> frm;
            if(frm.empty()) break;
            vector<int>::iterator ret;
            ret = find(v_frm_idx.begin(),v_frm_idx.end(),frm_idx);
            if(ret != v_frm_idx.end()){
                sprintf(strbuf, "%s/%s_%02d.jpg", opt.out_dir.c_str(), filename.c_str(), frm_idx);
                imwrite(strbuf,frm);
                njpg++;
            }
            frm_idx++;
        }
    }
    else{
        while(njpg < (int)v_frm_idx.size()){
            Mat frm;
            vr >> frm;
            if(frm.empty()) break;
            for (int shotid = 0; shotid < v_shot_range.size(); ++shotid) {

                vector<int>::iterator ret=find(v_shot_range[shotid].v_idx.begin(),
                        v_shot_range[shotid].v_idx.end(),frm_idx);
                if(ret != v_shot_range[shotid].v_idx.end()){
                    //cout<<frm_idx<<endl;
                    sprintf(strbuf,"%s/%s_%02d-%02d.jpg",opt.out_dir.c_str(), filename.c_str(), shotid, frm_idx);
                    //cout<<strbuf<<endl;
                    imwrite(strbuf,frm);
                    njpg++;
                }
            }
            frm_idx++;
        }
    }
    vr.release();

}