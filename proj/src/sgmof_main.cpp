#include <stdio.h>
#include <opencv2/opencv.hpp>
#include <utility>
#include "sgmof_main.h"
#include "epi_sgm.h"
#include "pyd_sgm.h"
#include "utils.h"

using namespace cv;

const String keys =
    "{help h usage ? |      | Call SGMOF to do optical flow for two images, Usage:\n ./SGMOF I1<first image> I2<second image> [-o]=<output flow file name> [-m]=0(epiSGM)/1(pydSGM)\n }"
    "{@I1 image1     |<none>| first image  }"
    "{@I2 image2     |<none>| second image   }"
    "{o outFile      |flow.png| output flow file (in KITTI format)}"
    "{m mode         |0     | epiSGM(0)/pydSGM mode(1) }"
	"{c calibFile    |calib.txt| calibration file, must have when mode = 0}"
	"{b benchmark    |0     | 0/1 for specifying Kitti2012/kitti2015 benchmark, used in EpiSGM only (calibration file has different format for 2012/2015)}"
    "{p passNum      |2     | number of SGM passes   }"
    "{d enableDiagonal |    | enable diagnoal directions in SGM }"
    "{V vzIndex      |      | enable vz-index in epipolar SGM    }"
    "{N pydNum       |5     |number of pyramidal level in PydSGM }"
;

//main entry to call Epi/Pyd SGM OF
int main(int argc, char** argv )
{
    CommandLineParser parser(argc, argv, keys);
    parser.about("SGM OF v0.0.1");

    if (parser.has("help"))
    {
        parser.printMessage();
        exit(0);
    }

    String image1FileName = parser.get<String>("image1");
    String image2FileName = parser.get<String>("image2");
    int mode = parser.get<int>("mode");
	String calibFileName = parser.get<String>("calibFile");
	int benchmark = parser.get<int>("benchmark");

    if (!parser.check())
    {
        parser.printErrors();
        exit(1);
    }

    Mat I1, I2, flow;
    I1 = imread(image1FileName, IMREAD_UNCHANGED);
    I2 = imread(image2FileName, IMREAD_UNCHANGED);
    if (I1.data == NULL || I2.data == NULL) {
        std::cout << "Open image failed..." << std::endl;
        exit(1);
    }
    else if (I1.size != I2.size) {
        std::cout << "Size of image1/2 must match" << std::endl;
    }
    
    if (mode == 0) {
		//check if calibration file provided
		Mat P = read_calib_file(calibFileName, benchmark == 1);
		//get intrinsic matrix
		Mat K = P(Range(0, 3), Range(0, 3));

        //call EpiSGM OF to calculate optical flow
        EpiSGM epiSGM;
        flow = epiSGM.compute(I1, I2);
    }
    else {
        //call PydSGM OF to calculate optical flow
        PydSGM pydSGM;
        flow = pydSGM.compute(I1, I2);
    }

    //write optical flow

    return 0;
}