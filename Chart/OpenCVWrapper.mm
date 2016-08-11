#include "OpenCVWrapper.h"

#import <opencv2/highgui/cap_ios.h>
#import "VideoProcessor.hpp"


using namespace cv;
// Class extension to adopt the delegate protocol
@interface CvVideoCameraWrapper () <CvVideoCameraDelegate>
{
    Ptr<VideoProcessor> changeDetector;
    double hr;
    int original;
}
@end
@implementation CvVideoCameraWrapper
{
    ViewController * viewController;
    UIImageView * imageView;
    CvVideoCamera * videoCamera;
    
}

- (double) getHR
{
    return hr;
}

- (void)setOriginal:(int)o
{
    original = o;
}

-(id)initWithController:(ViewController*)c andImageView:(UIImageView*)iv
{
    viewController = c;
    imageView = iv;
    
    videoCamera = [[CvVideoCamera alloc] initWithParentView:imageView];
    videoCamera.delegate = self;
    videoCamera.defaultAVCaptureDevicePosition = AVCaptureDevicePositionFront;
    videoCamera.defaultAVCaptureSessionPreset = AVCaptureSessionPreset352x288;
    videoCamera.defaultAVCaptureVideoOrientation = AVCaptureVideoOrientationPortrait;
    videoCamera.defaultFPS = 20;
    videoCamera.rotateVideo = true;
    [videoCamera start];
    changeDetector = new VideoProcessor(100);
    return self;
}
// This #ifdef ... #endif is not needed except in special situations
#ifdef __cplusplus
- (void)processImage:(Mat&)image
{
    Mat im_copy = image.clone();
    if (original == 1)
        changeDetector->colorMagnify(im_copy);
    else
        changeDetector->colorMagnify(image);
    hr = changeDetector->heartRate();
    //std::cout<<hr<<std::endl;
}
#endif
 
@end