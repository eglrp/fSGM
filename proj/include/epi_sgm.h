#ifndef __EPI_SGM_H__
#define __EPI_SGM_H__
#include <opencv2/opencv.hpp>
using namespace cv;

class EpiSGM
{
public:
    //main routine to caclulate optical flow for image1/image2
    //return a WXHx2 optical flow vector map
    Mat compute(Mat& I1, Mat& I2);
};

#endif