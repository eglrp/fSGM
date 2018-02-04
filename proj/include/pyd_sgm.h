#ifndef __PYD_SGM_H__
#define __PYD_SGM_H__

#include <opencv2/opencv.hpp>
using namespace cv;

class PydSGM
{
public:
    //main routine to caclulate optical flow for image1/image2
    //return a WXHx2 optical flow vector map
    Mat compute(Mat& I1, Mat& I2);
};


#endif __PYD_SGM_H__