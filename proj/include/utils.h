#ifndef __UTILS_H__
#define __UTILS_H__
#include <opencv2/opencv.hpp>
#include <string>
#include <fstream>
using namespace cv;
using namespace std;
#define M_PI 3.1415926

class FlowImage {
public:
	// default constructor
	FlowImage() {
		data_ = 0;
		width_ = 0;
		height_ = 0;
	}

	// construct flow image from png file
	FlowImage(const std::string file_name) {
		readFlowField(file_name);
	}

	// copy constructor
	FlowImage(const FlowImage &F) {
		width_ = F.width_;
		height_ = F.height_;
		data_ = (float*)malloc(width_*height_ * 3 * sizeof(float));
		memcpy(data_, F.data_, width_*height_ * 3 * sizeof(float));
	}

	// construct flow field from data
	FlowImage(const float* data, const int32_t width, const int32_t height) : width_(width), height_(height) {
		data_ = (float*)malloc(width*height * 3 * sizeof(float));
		memcpy(data_, data, width*height * 3 * sizeof(float));
	}

	// construct empty (= all pixels invalid) flow field of given width / height
	FlowImage(const int32_t width, const int32_t height) : width_(width), height_(height) {
		data_ = (float*)malloc(width*height * 3 * sizeof(float));
		for (int32_t i = 0; i<width*height * 3; i++)
			data_[i] = 0;
	}

	// deconstructor
	virtual ~FlowImage() {
		if (data_) { free(data_); data_ = 0; }
	}
	// read flow field from png file
	void read(const std::string file_name) {
		if (data_) { free(data_); data_ = 0; }
		readFlowField(file_name);
	}

	// write flow field to png file
	void write(const std::string file_name) {
		writeFlowField(file_name);
	}

	// write flow field to false color map using the Middlebury colormap
	void writeColor(const std::string file_name, float max_flow = -1.0f) {
		if (max_flow <= 1.0f)
			max_flow = std::max(maxFlow(), 1.0f);
		writeFalseColors(file_name, max_flow);
	}

	// get maximal optical flow magnitude
	float maxFlow() {
		float max_flow = 0;
		for (int32_t u = 0; u<width_; u++)
			for (int32_t v = 0; v<height_; v++)
				if (isValid(u, v) && getFlowMagnitude(u, v)>max_flow)
					max_flow = getFlowMagnitude(u, v);
		return max_flow;
	}
	// assignment operator, copies contents of F
	FlowImage& operator= (const FlowImage &F) {
		if (this != &F) {
			if (F.width_ != width_ || F.height_ != height_) {
				free(data_);
				width_ = F.width_;
				height_ = F.height_;
				data_ = (float*)malloc(width_*height_ * 3 * sizeof(float));
			}
			memcpy(data_, F.data_, width_*height_ * 3 * sizeof(float));
		}
		return *this;
	}

	// get optical flow u-component at given pixel
	inline float getFlowU(const int32_t u, const int32_t v) {
		return data_[3 * (v*width_ + u) + 0];
	}

	// get optical flow v-component at given pixel
	inline float getFlowV(const int32_t u, const int32_t v) {
		return data_[3 * (v*width_ + u) + 1];
	}

	// check if optical flow at given pixel is valid
	inline bool isValid(const int32_t u, const int32_t v) {
		return data_[3 * (v*width_ + u) + 2]>0.5;
	}

	// get optical flow magnitude at given pixel 
	inline float getFlowMagnitude(const int32_t u, const int32_t v) {
		float fu = getFlowU(u, v);
		float fv = getFlowV(u, v);
		return sqrt(fu*fu + fv*fv);
	}

	// set optical flow u-component at given pixel
	inline void setFlowU(const int32_t u, const int32_t v, const float val) {
		data_[3 * (v*width_ + u) + 0] = val;
	}

	// set optical flow v-component at given pixel
	inline void setFlowV(const int32_t u, const int32_t v, const float val) {
		data_[3 * (v*width_ + u) + 1] = val;
	}

	// set optical flow at given pixel to valid / invalid
	inline void setValid(const int32_t u, const int32_t v, const bool valid) {
		data_[3 * (v*width_ + u) + 2] = valid ? 1 : 0;
	}

	Mat FlowImage::errorImage(FlowImage &F_noc, FlowImage &F_occ, bool log_colors = false);

	// direct access to private variables
	float*  data() { return data_; }
	int32_t width() { return width_; }
	int32_t height() { return height_; }


private:
	float* data_;
	int width_;
	int height_;
private:
	//read KITTI format png flow file to 2D flow Matrix
	void readFlowField(const std::string fileName);
	
	//write 2D flow file to KITTI format png
	void writeFlowField(const std::string fileName);

	inline void hsvToRgb(float h, float s, float v, float &r, float &g, float &b) {
		float c = v*s;
		float h2 = 6.0*h;
		float x = c*(1.0 - fabs(fmod(h2, 2.0) - 1.0));
		if (0 <= h2&&h2<1) { r = c; g = x; b = 0; }
		else if (1 <= h2&&h2<2) { r = x; g = c; b = 0; }
		else if (2 <= h2&&h2<3) { r = 0; g = c; b = x; }
		else if (3 <= h2&&h2<4) { r = 0; g = x; b = c; }
		else if (4 <= h2&&h2<5) { r = x; g = 0; b = c; }
		else if (5 <= h2&&h2 <= 6) { r = c; g = 0; b = x; }
		else if (h2>6) { r = 1; g = 0; b = 0; }
		else if (h2<0) { r = 0; g = 1; b = 0; }
	}

	void writeFalseColors(const std::string file_name, const float max_flow) {
		float n = 8; // multiplier
		Mat image(height_, width_, CV_8UC3);
		for (int32_t v = 0; v<height_; v++) {
			uchar* tgt = image.ptr<uchar>(v);
			for (int32_t u = 0; u<width_; u++) {
				float r = 0, g = 0, b = 0;
				if (isValid(u, v)) {
					float mag = getFlowMagnitude(u, v);
					float dir = atan2(getFlowV(u, v), getFlowU(u, v));
					float h = fmod(dir / (2.0*M_PI) + 1.0, 1.0);
					float s = std::min(std::max(mag*n / max_flow, 0.0f), 1.0f);
					float v = std::min(std::max(n - s, 0.0f), 1.0f);
					hsvToRgb(h, s, v, r, g, b);
				}

				tgt[3 * u] = b*255.0f;
				tgt[3 * u + 1] = g*255.0f;
				tgt[3 * u + 2] = r*255.0f;
			}
		}
		imwrite(file_name, image);
	}
};

//read KITTI calibration file, return the projection matrix for cam0
Mat read_calib_file(string fileName, bool isKITTI2015 = false);


#endif