//
// Created by shawntang(唐兴) on 2019-02-20.
//

#ifndef FRAMEEXTRACTOR_GAPSTAT_H
#define FRAMEEXTRACTOR_GAPSTAT_H

#include "opencv2/opencv.hpp"
#include "opencv2/core/core.hpp"
#include "opencv2/imgproc/imgproc.hpp"

using namespace cv;

namespace extractor{

    inline void perform_kmeans(const Mat& km_data, Mat& km_lbl, Mat& km_ctr, int ncluster,
            int km_attempts = 1, int km_max_cnt = 1000, double km_eps = 0.0001){
        if(km_data.rows ==1){
            km_lbl = Mat::zeros(1,1,km_lbl.type());
            cv::reduce(km_data, km_ctr, 0, CV_REDUCE_AVG);
        }
        else{
            int km_k = min(ncluster, km_data.rows);
            TermCriteria km_opt = TermCriteria(CV_TERMCRIT_ITER|CV_TERMCRIT_EPS, km_max_cnt, km_eps);
            kmeans(km_data, km_k, km_lbl, km_opt, km_attempts, KMEANS_PP_CENTERS,km_ctr);
        }
    };


}
#endif //FRAMEEXTRACTOR_GAPSTAT_H
