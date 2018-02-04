#include "epi_sgm.h"

Mat EpiSGM::compute(Mat& I1, Mat& I2)
{
    return Mat::zeros(I1.rows, I1.cols, CV_32FC2);
}
