


#include "VideoProcessor.hpp"

VideoProcessor::VideoProcessor(int frameInterval)
: Interval(frameInterval)
, rate(50)
, fnumber(0)
, levels(4)
, alpha(500)
, lambda_c(80)
, fl(40.0/60.0)
, fh(80.0/60.0)
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
	cv::Mat temp1 = (1-fh)*lowpass1[curLevel] + fh*src;
	cv::Mat temp2 = (1-fl)*lowpass2[curLevel] + fl*src;
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
		   0, height - current.rows,
		   0, width - current.cols,
		   cv::BORDER_CONSTANT, cv::Scalar::all(0));

		// do the DFT
		cv::dft(tempImg, tempImg, cv::DFT_ROWS | cv::DFT_SCALE);

		// construct the filter
		cv::Mat filter = tempImg.clone();
		createIdealBandpassFilter(filter, fl, fh, rate);

		// apply filter
		cv::mulSpectrums(tempImg, filter, tempImg, cv::DFT_ROWS);

		// do the inverse DFT on filtered image
		cv::idft(tempImg, tempImg, cv::DFT_ROWS | cv::DFT_SCALE);

		// copy back to the current channel
		tempImg(cv::Rect(0, 0, current.cols, current.rows)).copyTo(channels[i]);
	}
	// merge channels
	cv::merge(channels, 3, dst);

	// normalize the filtered image
	cv::normalize(dst, dst, 0, 1, CV_MINMAX);
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
		dst = src * alpha;
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
	cv::Mat temp(frameSize.width*frameSize.height, Interval-1, CV_32FC3);
	for (int i = 0; i < Interval - 1; ++i) {
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
	for (int i = 0; i < Interval - 1; ++i) {    // get a line if any
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

	fl = 2 * fl * width / rate;
	fh = 2 * fh * width / rate;

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
    cv::Mat rgb;
    std::vector<cv::Mat> spl;
    cv::split(input,spl);
	// set filter
	setSpatialFilter(GAUSSIAN);
	setTemporalFilter(IDEAL);
    cvtColor(input, rgb, CV_RGBA2RGB);
	// output frame
	cv::Mat output;
	// motion image

	cv::Mat motion;
	// temp image
	cv::Mat temp;

	// filtered frames
	// std::deque<cv::Mat> filteredFrames;

	// concatenate image of all the down-sample frames
	cv::Mat videoMat;
	// concatenate filtered image
	cv::Mat filtered;
// std:cout << input.cols << std::endl;
    if (frames.size() < Interval) {
        rgb.convertTo(temp, CV_32FC3);
        frames.push_back(temp.clone());
        // spatial filtering
        std::deque<cv::Mat> pyramid;
        spatialFilter(temp, pyramid);
        downSampledFrames.push_back(pyramid.at(levels-1));
    }
    else{
        rgb.convertTo(temp, CV_32FC3);
        frames.push_back(temp.clone());
        // spatial filtering
        std::deque<cv::Mat> pyramid;
        spatialFilter(temp, pyramid);
        downSampledFrames.push_back(pyramid.at(levels-1));
			concat(downSampledFrames, videoMat);

			// 3. temporal filtering
			temporalFilter(videoMat, filtered);

			// 4. amplify color motion
			amplify(filtered, filtered);

			std::deque<cv::Mat> filteredFrames;
			// 5. de-concat the filtered image into filtered frames
			deConcat(filtered, downSampledFrames.at(0).size(), filteredFrames);

			// fnumber = 0;
			// up-sample the motion image        
			upsamplingFromGaussianPyramid(filteredFrames.at(Interval / 2), levels, motion);
			resize(motion, motion, frames.at(Interval / 2).size());
			temp = frames.back() + motion;
			output = temp.clone();
			double minVal, maxVal;
			minMaxLoc(output, &minVal, &maxVal); //find minimum and maximum intensities
			output.convertTo(output, CV_8UC3, 255.0/(maxVal - minVal), -minVal * 255.0/(maxVal - minVal));
        
            std::vector<cv::Mat> spl_new;
            cv::split(output,spl_new);
            std::vector<cv::Mat> out = {spl_new[0],spl_new[1],spl_new[2],spl[3]};
            cv::merge(out,input);

            frames.pop_front();
			downSampledFrames.pop_front();
        
	}
}


