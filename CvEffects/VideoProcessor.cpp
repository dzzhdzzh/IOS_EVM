

#include "VideoProcessor.hpp"

VideoProcessor::VideoProcessor(int frameInterval)
: Interval(frameInterval)
, Frate(30)
, fnumber(0)
, visualflag(false)
, hrnumber(400)
, levels(4)
, alpha(200)
, lambda_c(80)
, Fl(40.0/60.0)
, Fh(120.0/60.0)
, chromAttenuation(1)
, delta(0)
, exaggeration_factor(2.0)
, lambda(0)
{}

/**
 * spatialFilter	-	spatial filtering an image
 *
 * @param src		-	source image
 * @param pyramid	-	destinate pyramid
 */
bool VideoProcessor::spatialFilter(const cv::Mat &src, std::deque<cv::Mat> &pyramid)
{
    switch (spatialType) {
        case LAPLACIAN:     // laplacian pyramid
            return buildLaplacianPyramid(src, levels, pyramid);
            break;
        case GAUSSIAN:      // gaussian pyramid
            return buildGaussianPyramid(src, levels, pyramid);
            break;
        default:
            return false;
            break;
    }
}

/**
 * temporalFilter	-	temporal filtering an image
 *
 * @param src	-	source image
 * @param dst	-	destinate image
 */
void VideoProcessor::temporalFilter(const cv::Mat &src,
                                    cv::Mat &dst,cv::Mat &dst_norm)
{
    switch(temporalType) {
        case IIR:       // IIR bandpass filter
            temporalIIRFilter(src, dst);
            break;
        case IDEAL:     // Ideal bandpass filter
            temporalIdealFilter(src, dst, dst_norm);
            break;
        default:
            break;
    }
    return;
}

/**
 * temporalIIRFilter	-	temporal IIR filtering an image
 *                          (thanks to Yusuke Tomoto)
 * @param pyramid	-	source image
 * @param filtered	-	filtered result
 *
 */
void VideoProcessor::temporalIIRFilter(const cv::Mat &src,
                                       cv::Mat &dst)
{
    cv::Mat temp1 = (1-Fh)*lowpass1[curLevel] + Fh*src;
    cv::Mat temp2 = (1-Fl)*lowpass2[curLevel] + Fl*src;
    lowpass1[curLevel] = temp1;
    lowpass2[curLevel] = temp2;
    dst = lowpass1[curLevel] - lowpass2[curLevel];
}

/**
 * temporalIdalFilter	-	temporal IIR filtering an image pyramid of concat-frames
 *                          (Thanks to Daniel Ron & Alessandro Gentilini)
 *
 * @param pyramid	-	source pyramid of concatenate frames
 * @param filtered	-	concatenate filtered result
 *
 */
void VideoProcessor::temporalIdealFilter(const cv::Mat &src,
                                         cv::Mat &dst, cv::Mat &dst_norm)
{
    cv::Mat channels[3];
    
    // split into 3 channels
    cv::split(src, channels);
    
    for (int i = 0; i < 3; ++i){
        
        cv::Mat current = channels[i];  // current channel
        cv::Mat tempImg;
        
        int width = cv::getOptimalDFTSize(current.cols);
        int height = cv::getOptimalDFTSize(current.rows);
        
        cv::copyMakeBorder(current, tempImg,
                           0, 0,
                           0, width - current.cols,
                           cv::BORDER_CONSTANT, cv::Scalar::all(0));
        
        cv::Mat planes[] = {cv::Mat_<float>(tempImg), cv::Mat::zeros(tempImg.size(), CV_32F)};
        cv::Mat complexI;
        merge(planes, 2, complexI);
//        cv::dft(complexI, complexI, cv::DFT_ROWS|cv::DFT_COMPLEX_OUTPUT|cv::DFT_SCALE);
//        cv::split(complexI, planes);

        // do the DFT
        cv::dft(complexI, complexI, cv::DFT_ROWS | cv::DFT_SCALE);
        
        // construct the filter
        cv::Mat filter = tempImg.clone();
        createIdealBandpassFilter(filter, Fl, Fh, Frate);
        //std::cout<<filter.rowRange(0, 1)<<std::endl;
        // apply filter
        cv::Mat filterplanes[] = {cv::Mat_<float>(filter), cv::Mat_<float>(filter)};
        cv::Mat complexfilter,roifft;
        cv::Mat avgfft;
        //cv::Mat avgfft(floor(imgGaussianSize.width*0.7)*floor(imgGaussianSize.height*0.7), complexI.cols, CV_32FC1);
        merge(filterplanes, 2, complexfilter);
        cv::mulSpectrums(complexI, complexfilter, complexI, 0);

        if(i==2)
        {
            cv::split(complexI, planes);
            cv::magnitude(planes[0], planes[1], planes[0]);
            cv::Mat magI = planes[0];
            
            //std::cout<<magI<<std::endl;
            /*
            for(int k=0;k<magI.cols;k++)
            {
                roifft = magI.col(k).clone();
                roifft = roifft.reshape(0,imgGaussianSize.height);
                cv::Rect roi = cv::Rect(0.15*imgGaussianSize.width, 0.15*imgGaussianSize.height, imgGaussianSize.width*0.7, imgGaussianSize.height*0.7);
                roifft = roifft(roi).clone();
                roifft = roifft.reshape(0,roifft.cols*roifft.rows).clone();
                cv::Mat line = avgfft.col(k);
                roifft.copyTo(line);
            }
            */
            //std::cout<<avgfft<<std::endl;
            reduce(magI,avgfft, 0, CV_REDUCE_AVG);
            cv::Point maxLoc;
            double maxVal;
            cv::minMaxLoc(avgfft.colRange(0, avgfft.cols / 2.0 - 1), NULL, &maxVal, NULL, &maxLoc);
            //std::cout<< maxVal <<"  "<<maxLoc.x<<"   "<<60.0 / avgfft.cols * Frate * (maxLoc.x)<<std::endl;
            if(maxVal<0)//0.07)
                visualflag = false;
            else
                visualflag = true;
        }
        
        // do the inverse DFT on filtered image
        cv::idft(complexI, complexI, cv::DFT_ROWS);
        cv::split(complexI, planes);
        // copy back to the current channel
        planes[0](cv::Rect(0, 0, current.cols, current.rows)).copyTo(channels[i]);
    }
    // merge channels
    cv::merge(channels, 3, dst);
    
    // normalize the filtered image
    cv::normalize(dst, dst_norm, 0, 1, CV_MINMAX);
}

/**
 * amplify	-	ampilfy the motion
 *
 * @param filtered	- motion image
 */
void VideoProcessor::amplify(const cv::Mat &src, cv::Mat &dst)
{
    float currAlpha;
    switch (spatialType) {
        case LAPLACIAN:
            //compute modified alpha for this level
            currAlpha = lambda/delta/8 - 1;
            currAlpha *= exaggeration_factor;
            if (curLevel==levels || curLevel==0)     // ignore the highest and lowest frequency band
                dst = src * 0;
            else
                dst = src * cv::min(alpha, currAlpha);
            break;
        case GAUSSIAN:
//            cv::Mat tmp;
            dst = src.clone() * alpha;
            break;
        default:
            break;
    }
}

/**
 * attenuate	-	attenuate I, Q channels
 *
 * @param src	-	source image
 * @param dst   -   destinate image
 */
void VideoProcessor::attenuate(cv::Mat &src, cv::Mat &dst)
{
    cv::Mat planes[3];
    cv::split(src, planes);
    planes[0] = planes[0] * chromAttenuation;
    planes[1] = planes[1] * chromAttenuation;
    cv::merge(planes, 3, dst);
}


/**
 * concat	-	concat all the frames into a single large Mat
 *              where each column is a reshaped single frame
 *
 * @param frames	-	frames of the video sequence
 * @param dst		-	destinate concatnate image
 */
void VideoProcessor::concat(const std::deque<cv::Mat> &frames,
                            cv::Mat &dst)
{
    cv::Size frameSize = frames.at(0).size();
    cv::Mat temp(frameSize.width*frameSize.height, frames.size(), CV_32FC3);
    for (int i = 0; i < frames.size(); ++i) {
        // get a frame if any
        cv::Mat out = frames.at(i);
        // reshape the frame into one column
        cv::Mat reshaped = out.reshape(3, out.cols*out.rows).clone();
        cv::Mat line = temp.col(i);
        // save the reshaped frame to one column of the destinate big image
        reshaped.copyTo(line);
    }
    temp.copyTo(dst);
}

/**
 * deConcat	-	de-concat the concatnate image into frames
 *
 * @param src       -   source concatnate image
 * @param framesize	-	frame size
 * @param frames	-	destinate frames
 */
void VideoProcessor::deConcat(const cv::Mat &src,
                              const cv::Size &frameSize,
                              std::vector<cv::Mat> &frames)
{
    for (int i = 0; i < src.cols; ++i) {    // get a line if any
        cv::Mat line = src.col(i).clone();
        cv::Mat reshaped = line.reshape(3, frameSize.height).clone();
        frames.push_back(reshaped);
    }
}

/**
 * createIdealBandpassFilter	-	create a 1D ideal band-pass filter
 *
 * @param filter    -	destinate filter
 * @param fl        -	low cut-off
 * @param fh		-	high cut-off
 * @param rate      -   sampling rate(i.e. video frame rate)
 */
void VideoProcessor::createIdealBandpassFilter(cv::Mat &filter, double fl, double fh, double rate)
{
    int width = filter.cols;
    int height = filter.rows;
    
    fl = fl * width / rate;
    fh = fh * width / rate;
    
    double response;
    
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            // filter response
            //if (j >= fl && j <= fh)
            if ((j >= fl && j <= fh) || ((j >= width-fh-1) && (j<= width-fl-1)))
                response = 1.0f;
            else
                response = 0.0f;
            filter.at<float>(i, j) = response;
        }
    }
}


/**
 * setSpatialFilter	-	set the spatial filter
 *
 * @param type	-	spatial filter type. Could be:
 *					1. LAPLACIAN: laplacian pyramid
 *					2. GAUSSIAN: gaussian pyramid
 */
void VideoProcessor::setSpatialFilter(spatialFilterType type)
{
    spatialType = type;
}

/**
 * setTemporalFilter	-	set the temporal filter
 *
 * @param type	-	temporal filter type. Could be:
 *					1. IIR: second order(IIR) filter
 *					2. IDEAL: ideal bandpass filter
 */
void VideoProcessor::setTemporalFilter(temporalFilterType type)
{
    temporalType = type;
}


/**
 * colorMagnify	-	color magnification
 *
 */
void VideoProcessor::colorMagnify(cv::Mat& input)
{
    Frate = 1.0 / (((double)cv::getTickCount() - (double) t_name) / (double) cv::getTickFrequency());
    t_name = cv::getTickCount();
    //std::cout << "frame rate: " << Frate << std::endl;
    cv::Mat rgb, output;
    cv::Mat channelA(input.rows, input.cols, CV_8UC1, 1);
    cv::extractChannel(input, channelA, 3);
    // set filter
    setSpatialFilter(GAUSSIAN);
    setTemporalFilter(IDEAL);
    
    cvtColor(input, rgb, CV_RGBA2RGB);
    // motion image
    cv::Mat motion;
    
    // concatenate image of all the down-sample frames
    cv::Mat videoMat;
    // concatenate filtered image
    cv::Mat filtered;
    cv::Mat filtered_norm;
    cv::Mat R;
    //std::cout << "new: " << Frate << std::endl;
    if (frames.size() < Interval - 1) {
        rgb.convertTo(output, CV_32FC3);
        frames.push_back(output.clone());
        // spatial filtering
        std::deque<cv::Mat> pyramid;
        spatialFilter(output, pyramid);
        downSampledFrames.push_back(pyramid.at(levels-1));
    }
    else{
        rgb.convertTo(output, CV_32FC3);
        frames.push_back(output.clone());
        // spatial filtering
        std::deque<cv::Mat> pyramid;
        spatialFilter(output, pyramid);
        downSampledFrames.push_back(pyramid.at(levels-1));
        concat(downSampledFrames, videoMat);
        
        imgGaussianSize = pyramid.at(levels-1).size();
        
        // push the HR Frames queue
        //if(fnumber==0)
            hrFrames.push_back(pyramid.at(levels-1));
        if(hrFrames.size() > hrnumber)
        {
            startHR = 1;
            hrFrames.pop_front();
            cv::Mat tempMat;
            concat(hrFrames, tempMat);
            cv::extractChannel(tempMat, R, 2);
            //hrMat = R;
            cv::extractChannel(tempMat, hrMat, 1);
            hrMat = hrMat - R;
        }
        else
            startHR = 0;

        
        // 3. temporal filtering
        temporalFilter(videoMat, filtered, filtered_norm);
        
        // 4. amplify color motion
        amplify(filtered_norm, filtered_norm);
        
        filterednormFrames.clear();
        filteredFrames.clear();
        
        // 5. de-concat the filtered image into filtered frames
        deConcat(filtered_norm, downSampledFrames.at(0).size(), filterednormFrames);
        
        
        // up-sample the motion image
        upsamplingFromGaussianPyramid(filterednormFrames.back(), levels, motion);
        
        // resize(filteredFrames.back(), motion, frames.back().size());
        output = frames.back() + motion;
        
        if(visualflag)
        {
            double minVal, maxVal;
            minMaxLoc(output, &minVal, &maxVal); //find minimum and maximum intensities
            output.convertTo(output, CV_8UC3, 255.0/(maxVal - minVal),-minVal * 255.0/(maxVal - minVal));
        }
        else
            frames.back().convertTo(output,CV_8UC3);
        
        std::vector<cv::Mat> out;
        cv::split(output, out);
        out.push_back(channelA);
        cv::merge(out, input);
        frames.pop_front();
        downSampledFrames.pop_front();
    }
    fnumber++;
}


void VideoProcessor::fft1DMag(cv::Mat& src, cv::Mat& dst){
    cv::Mat padded;                            //expand input image to optimal size
    int n = cv::getOptimalDFTSize( src.cols ); // on the border add zero pixels
    copyMakeBorder(src, padded, 0, 0, 0, n - src.cols, cv::BORDER_CONSTANT, cv::Scalar::all(0));
    cv::Mat planes[] = {cv::Mat_<float>(padded), cv::Mat::zeros(padded.size(), CV_32F)};
    cv::Mat complexI;
    merge(planes, 2, complexI);
    cv::dft(complexI, complexI, cv::DFT_ROWS|cv::DFT_COMPLEX_OUTPUT);
    
    // construct the filter
    cv::Mat filter = padded.clone();
    createIdealBandpassFilter(filter, Fl, Fh, Frate);
    // apply filter
    cv::Mat filterplanes[] = {cv::Mat_<float>(filter), cv::Mat_<float>(filter)};
    cv::Mat complexfilter;
    merge(filterplanes, 2, complexfilter);
    cv::mulSpectrums(complexI, complexfilter, complexI, 0);
    
    cv::split(complexI, planes);
    cv::magnitude(planes[0], planes[1], planes[0]);
    cv::Mat magI = planes[0];
    reduce(magI,dst, 0, CV_REDUCE_AVG);
    //std::cout<<dst<<std::endl;
    //dst.at<float>(0, 0) = 0;
}


double VideoProcessor::heartRate(){
    cv::Mat rowMean;
    //    std::cout << Frate << std::endl;
    if (startHR)
    {
        fft1DMag(hrMat, rowMean);
//        std::cout << rowMean << std::endl;
        cv::Point maxLoc;
        double maxVal;
        cv::minMaxLoc(rowMean.colRange(0, rowMean.cols / 2.0 - 1), NULL, &maxVal, NULL, &maxLoc);
        //std::cout<<HR<<" "<<maxVal<<" "<<maxLoc.x<<std::endl;
        double HR_current;
        HR_current = 60.0 / rowMean.cols * Frate * (maxLoc.x);
        hrResult.push_back(HR_current);
        if (hrResult.size()>30)
        {
            hrResult.pop_front();
            double tempsum = 0;
            for(int k=0;k<hrResult.size();k++)
                tempsum += hrResult[k];
            HR = tempsum/hrResult.size();
            std::cout << maxVal <<"  "<<maxLoc.x << "   "<<HR_current<<"  "<<HR<<std::endl;
            return HR;
        }
        else
            return 0.0;
    }
    else
        return 0.0;
}



