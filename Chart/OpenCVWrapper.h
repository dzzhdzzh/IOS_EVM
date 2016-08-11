#import <UIKit/UIKit.h>
#import <Foundation/Foundation.h>

@class ViewController;

@interface CvVideoCameraWrapper : NSObject

- (id) initWithController:(ViewController*)c andImageView:(UIImageView*)iv;
- (double) getHR;
- (void)setOriginal:(int)o;


@end