//
// Created by shawntang(唐兴) on 2019-02-14.
//

#ifndef FRAMEEXTRACTOR_IMAGE_METRIC_H
#define FRAMEEXTRACTOR_IMAGE_METRIC_H

#include <cmath>
#include <vector>
#include <math.h>


#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/opencv.hpp>

#include "extractor/hist_opencv.h"

using namespace std;
using namespace cv;

namespace extractor{

    #define VALIDATE(x) (std::isnan(x)) ? 0 : x


    inline double calc_sharpness(const Mat& gray_img)
    {
        Mat img;
        gray_img.convertTo(img, CV_32FC1);
        img *= 1./255;

        Mat dx, dy;
        Sobel(img, dx, img.type(), 1, 0, 3);
        Sobel(img, dy, img.type(), 0, 1, 3);
        magnitude(dx, dy, dx);

        int npixels = gray_img.rows * gray_img.cols;

        return VALIDATE(cv::sum(dx)[0] / npixels);
    }

    inline double calc_brightness(const Mat& img)
    {
        vector<Mat> bgr;
        cv::split(img, bgr);
        for(size_t j = 0; j < bgr.size(); j++)
            bgr[j].convertTo(bgr[j],CV_32F);

        Mat img_pb = (0.2126*bgr[2] + 0.7152*bgr[1] + 0.0722*bgr[0])/255.0;

        return VALIDATE(mean(img_pb)[0]);
    }

    inline double calc_uniformity(const Mat& gray_img, int nbins=128)
    {
        double val = 0.0;

        Mat ghist;
        extractor::calc_gray_hist(gray_img, ghist, nbins);
        if(cv::sum(ghist)[0] == 0){
            val = 1.0;
        }
        else{
            Mat hist_sorted;
            cv::sort(ghist, hist_sorted, CV_SORT_EVERY_COLUMN|CV_SORT_DESCENDING);
            hist_sorted /= cv::sum(hist_sorted)[0];
            for(int j = 0; j < (int)(0.05*nbins); j++)
                val += (double) hist_sorted.at<float>(j);
        }
        return VALIDATE(val);
    }


}
#endif //FRAMEEXTRACTOR_IMAGE_METRIC_H
