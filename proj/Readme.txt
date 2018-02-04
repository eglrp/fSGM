This is the C++ algorithm project for Epipolar SGM/Pyramidal SGM OF methods. The project is based on OpenCV, the opencv library is a pre-requisite. the code structure for this project is as following:
include/ has the main function/class definitions;
src/ has the main implementations for these two methods
external/ contains the external codes which is used in this project

To build the solution (cmake is required to generate the build solution)
1. change the opencv path to your local path (should contains OpenCVConfig.cmake, OpenCVConfig-version.cmake under that path) in LIST(APPEND CMAKE_PREFIX_PATH "C:/opencv/opencv-3.2.0/build")
2. build solution

    Linux: 
    mkdir build; cd build; cmake ../; make;
    
    Windows:
    mkdir build; cd build; cmake ../; then open VS solution and build all