//
// Created by shawntang(唐兴) on 2019-02-15.
//

#ifndef FRAMEEXTRACTOR_VIDEO_PARSER_H
#define FRAMEEXTRACTOR_VIDEO_PARSER_H

#include <string>
#include <vector>
#include <math.h>
#include <numeric>
#include <random>
#include <ctime>
#include <algorithm>

#if defined(_OPENMP)
#include <omp.h>
#endif

#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>



#include "extractor/shot_range.h"
#include "extractor/image_metric.h"
#include "extractor/sort.h"
#include "extractor/hist_opencv.h"
#include "extractor/gapstat.h"

using namespace std;
using namespace cv;

namespace extractor {

    struct video_metadata{
        int nframes;
        int width;
        int height;
        double duration;
        double fps;
    };

    struct parser_params{
        int step_sz;
        double fltr_begin_sec;
        double fltr_end_sec;
        bool fltr_lq;
        bool fltr_rdt;
        bool sample;

        parser_params():
                step_sz(1),
                fltr_begin_sec(0),
                fltr_end_sec(0),
                fltr_lq(true),
                fltr_rdt(true),
                //still(false),
                sample(false)
        {};
    };

    inline void print_video_metadata(const string filename, const extractor::video_metadata m){
        printf("%s\n seconds=%.2f, nframes=%d, fps=%.2f, resolution=%dx%d\n",
                filename.c_str(), m.duration,m.nframes,m.fps,m.width,m.height);
    };

    class VideoParser{
    public:
        VideoParser();
        vector<extractor::ShotRange> parse_video(const string& in_video, extractor::parser_params opt);
        int get_nfrm_valid();

    private:
        int read_video(const string& in_video,int step_sz=1);
        void filter_heuristic(double fltr_begin_sec = 0, double fltr_end_sec = 0);
        void filter_low_quality(double thrsh_bright = 0.075, double thrsh_sharp = 0.08, double thrsh_uniform = 0.80);
        void filter_transition(double thrsh_diff=0.50, double thrsh_ecr=0.10);
        void extract_histo_features(int pyr_level = 2, bool omit_filtered = true, int nbin_color_hist=128,
                int nbin_edge_ori_hist=8, int nbin_edge_mag_hist = 8);
        void post_process(double min_shot_sec = 2.0);
        void update_shot_range(int min_shot_len = 5);
        void filter_redundant_and_obtain_subshots();
        void sbd_heuristic(vector<double> v_diff, vector<int>& jump, int njumps, int min_shot_len);
        void mark_invalid( vector<bool>& vec, int idx, int wnd_sz=0);
        void subsample_sec();
        void release_memory();

    public:
        extractor::video_metadata meta;
        vector<bool> _v_frm_valid;

    private:
        int _step_sz;
        int _nfrm_total;
        int _nfrm_given;
        int _video_width;
        int _video_height;
        double _video_fps;
        double _video_sec;

        Mat _X_feat;
        Mat _X_diff;
        Mat _X_ecr;

        vector<Mat> _v_frm_rgb;
        vector<Mat> _v_frm_gray;
        vector<extractor::ShotRange> _v_shot_ranges;

    };
}

#endif //FRAMEEXTRACTOR_VIDEO_PARSER_H
