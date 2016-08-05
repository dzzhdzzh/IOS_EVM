


#include "VideoProcessor.hpp"

VideoProcessor::VideoProcessor(int frameInterval)
: Interval(frameInterval)
, Frate(30)
, fnumber(0)
, levels(4)
, alpha(200)
, lambda_c(80)
, Fl(40.0/60.0)
, Fh(80.0/60.0)
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
                                    cv::Mat &dst)
{
    switch(temporalType) {
        case IIR:       // IIR bandpass filter
            temporalIIRFilter(src, dst);
            break;
        case IDEAL:     // Ideal bandpass filter
            temporalIdealFilter(src, dst);
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
                                         cv::Mat &dst)
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
        
        // apply filter
        cv::Mat filterplanes[] = {cv::Mat_<float>(filter), cv::Mat_<float>(filter)};
        cv::Mat complexfilter;
        merge(filterplanes, 2, complexfilter);
        cv::mulSpectrums(complexI, complexfilter, complexI, 0);
        
        // do the inverse DFT on filtered image
        cv::idft(complexI, complexI, cv::DFT_ROWS);
        cv::split(complexI, planes);
        // copy back to the current channel
        planes[0](cv::Rect(0, 0, current.cols, current.rows)).copyTo(channels[i]);
    }
    // merge channels
    cv::merge(channels, 3, dst);
    
    // normalize the filtered image
//    cv::normalize(dst, dst, 0, 1, CV_MINMAX);
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
    planes[1] = planes[1] * chromAttenuation;
    planes[2] = planes[2] * chromAttenuation;
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
    cv::Mat temp(frameSize.width*frameSize.height, Interval, CV_32FC3);
    for (int i = 0; i < Interval; ++i) {
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
                              std::deque<cv::Mat> &frames)
{
    for (int i = 0; i < Interval; ++i) {    // get a line if any
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
            if (j >= fl && j <= fh)
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
 * motionMagnify	-	eulerian motion magnification
 *
 */
//void VideoProcessor::motionMagnify(cv::Mat& input)
// {
//	// set filter
//	setSpatialFilter(LAPLACIAN);
//	setTemporalFilter(IIR);
//
//	// output frame
//	cv::Mat output;
//
//	// motion image
//	cv::Mat motion;
//
//	std::deque<cv::Mat> pyramid;
//	std::deque<cv::Mat> filtered;
//
//
//	while (true) {
//
//		cv::resize(input, input, cv::Size(320, 180), 0, 0,   cv::INTER_LINEAR);
////		imshow("input", input);
//		input.convertTo(input, CV_32FC3, 1.0/255.0f);
//
//		// 1. convert to Lab color space
//		cv::cvtColor(input, input, CV_BGR2Lab);
//
//		// 2. spatial filtering one frame
//		cv::Mat s = input.clone();
//		spatialFilter(s, pyramid);
//
//		// 3. temporal filtering one frame's pyramid
//		// and amplify the motion
//		if (fnumber == 0){      // is first frame
//			lowpass1 = pyramid;
//			lowpass2 = pyramid;
//			filtered = pyramid;
//		} else {
//			for (int i=0; i<levels; ++i) {
//				curLevel = i;
//				temporalFilter(pyramid.at(i), filtered.at(i));
//			}
//
//			// amplify each spatial frequency bands
//			// according to Figure 6 of paper
//			cv::Size filterSize = filtered.at(0).size();
//			int w = filterSize.width;
//			int h = filterSize.height;
//
//			delta = lambda_c/8.0/(1.0+alpha);
//			// the factor to boost alpha above the bound
//			// (for better visualization)
//			exaggeration_factor = 2.0;
//
//			// compute the representative wavelength lambda
//			// for the lowest spatial frequency band of Laplacian pyramid
//			lambda = sqrt(w*w + h*h)/3;  // 3 is experimental constant
//
//			for (int i=levels; i>=0; i--) {
//				curLevel = i;
//
//				amplify(filtered.at(i), filtered.at(i));
//
//				// go one level down on pyramid
//				// representative lambda will reduce by factor of 2
//				lambda /= 2.0;
//			}
//		}
//
//		// 4. reconstruct motion image from filtered pyramid
//		reconImgFromLaplacianPyramid(filtered, levels, motion);
//
//		// 5. attenuate I, Q channels
//		attenuate(motion, motion);
//
//		// 6. combine source frame and motion image
//		if (fnumber > 0)    // don't amplify first frame
//			s += motion;
//
//		// 7. convert back to rgb color space and CV_8UC3
//		output = s.clone();
//		cv::cvtColor(output, output, CV_Lab2BGR);
//		output.convertTo(output, CV_8UC3, 255.0, 1.0/255.0);
//
//		// write the frame to the temp file
//		// tempWriter.write(output);
//
//		// update process
////		std::string msg= "Processing...";
////		imshow("output", output);
////		cv::waitKey(1);
//		fnumber++;
//	}
//}

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
    //std::cout << "new: " << Frate << std::endl;
    if (frames.size() < Interval - 1) {
        rgb.convertTo(output, CV_32FC3,1./255.);
        frames.push_back(output.clone());
        // spatial filtering
        std::deque<cv::Mat> pyramid;
        spatialFilter(output, pyramid);
        downSampledFrames.push_back(pyramid.at(levels-1));
    }
    else{
        startHR = 1;
        rgb.convertTo(output, CV_32FC3,1./255.);
        frames.push_back(output.clone());
        // spatial filtering
        std::deque<cv::Mat> pyramid;
        spatialFilter(output, pyramid);
        downSampledFrames.push_back(pyramid.at(levels-1));
        concat(downSampledFrames, videoMat);
        
        // 3. temporal filtering
        temporalFilter(videoMat, filtered);
        
        // 4. amplify color motion
        amplify(filtered, filtered);
        //        hrMat = filtered.clone();
        cv::Mat R;
        cv::extractChannel(filtered, R, 2);
        cv::extractChannel(filtered, hrMat, 1);
        hrMat = hrMat - R;
        filteredFrames.clear();
        
        //std::cout<<cv::mean(R)[0]<<std::endl;
        
        // 5. de-concat the filtered image into filtered frames
        deConcat(filtered, downSampledFrames.at(0).size(), filteredFrames);
        //        hrMat = filteredFrames.back().clone();
        // up-sample the motion image
        upsamplingFromGaussianPyramid(filteredFrames.back(), levels, motion);
        
//        resize(filteredFrames.back(), motion, frames.back().size());
        output = frames.back() + motion;
        double minVal, maxVal;
        minMaxLoc(output, &minVal, &maxVal); //find minimum and maximum intensities
        output.convertTo(output, CV_8UC3, 255.0);
        std::vector<cv::Mat> out;
        cv::split(output, out);
        out.push_back(channelA);
        cv::merge(out, input);
        frames.pop_front();
        downSampledFrames.pop_front();
    }
}


void fft1DMag(cv::Mat& src, cv::Mat& dst){
    cv::Mat padded;                            //expand input image to optimal size
    int n = cv::getOptimalDFTSize( src.cols ); // on the border add zero pixels
    copyMakeBorder(src, padded, 0, 0, 0, 512 - src.cols, cv::BORDER_CONSTANT, cv::Scalar::all(0));
    cv::Mat planes[] = {cv::Mat_<float>(padded), cv::Mat::zeros(padded.size(), CV_32F)};
    cv::Mat complexI;
    merge(planes, 2, complexI);
    cv::dft(complexI, complexI, cv::DFT_ROWS|cv::DFT_COMPLEX_OUTPUT|cv::DFT_SCALE);
    cv::split(complexI, planes);
    cv::magnitude(planes[0], planes[1], planes[0]);
    cv::Mat magI = planes[0];
    reduce(magI,dst, 0, CV_REDUCE_AVG);
//    dst.at<float>(0, 0) = 0;
}


double VideoProcessor::heartRate(){
    cv::Mat rowMean;
    //    std::cout << Frate << std::endl;
    if (startHR){
        fft1DMag(hrMat, rowMean);
//        std::cout << rowMean << std::endl;
        cv::Point maxLoc;
        double maxVal;
        cv::minMaxLoc(rowMean.colRange(0, rowMean.cols / 2.0 - 1), NULL, &maxVal, NULL, &maxLoc);
        //std::cout<<maxVal<<" "<<maxLoc.x<<std::endl;
        HR = 60.0 / rowMean.cols * Frate * (maxLoc.x);
        std::cout << HR << std::endl;
        return HR;
    }
    else
        return 0.0;
}



