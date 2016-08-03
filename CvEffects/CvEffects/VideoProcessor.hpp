// Yet anther C++ implementation of EVM, based on OpenCV and Qt. 
// Copyright (C) 2014  Joseph Pan <cs.wzpan@gmail.com>
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
// 
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
// 02110-1301 USA
// 

#pragma once

#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <deque>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include "SpatialFilter.hpp"

enum spatialFilterType {LAPLACIAN, GAUSSIAN};
enum temporalFilterType {IIR, IDEAL};

class VideoProcessor {

public:

    VideoProcessor(int frameInterval);

    // set spatial filter
    void setSpatialFilter(spatialFilterType type);

    // set temporal filter
    void setTemporalFilter(temporalFilterType type);

    // motion magnification
    void motionMagnify(cv::Mat& input);

    // color magnification
    void colorMagnify(cv::Mat& input);


private:    

    int Interval;
    
    // video frame rate
    double rate;
    // number of processed frames
    long fnumber;
    // current level of pyramid
    int curLevel;
    // spatial filter type
    spatialFilterType spatialType;
    // temporal filter type
    temporalFilterType temporalType;
    // level numbers of image pyramid
    int levels;
    // amplification factor
    float alpha;
    // cut-off wave length
    float lambda_c;    
    // low cut-off
    float fl;
    // high cut-off
    float fh;
    // chromAttenuation
    float chromAttenuation;
    // delta
    float delta;
    // extraggon factor
    float exaggeration_factor;
    // lambda
    float lambda;

    // low pass filters for IIR
    std::deque<cv::Mat> lowpass1;
    std::deque<cv::Mat> lowpass2;


    // spatial filtering
    bool spatialFilter(const cv::Mat &src, std::deque<cv::Mat> &pyramid);

    // temporal filtering
    void temporalFilter(const cv::Mat &src,
                        cv::Mat &dst);

    // temporal IIR filtering
    void temporalIIRFilter(const cv::Mat &src,
                        cv::Mat &dst);

    // temporal ideal bandpass filtering
    void temporalIdealFilter(const cv::Mat &src,
                             cv::Mat &dst);

    // amplify motion
    void amplify(const cv::Mat &src, cv::Mat &dst);

    // attenuate I, Q channels
    void attenuate(cv::Mat &src, cv::Mat &dst);

    // concat images into a large Mat
    void concat(const std::deque<cv::Mat> &frames, cv::Mat &dst);

    // de-concat the concatnate image into frames
    void deConcat(const cv::Mat &src, const cv::Size &frameSize, std::deque<cv::Mat> &frames);

    // create an ideal bandpass processor
    void createIdealBandpassFilter(cv::Mat &filter, double fl, double fh, double rate);
};


