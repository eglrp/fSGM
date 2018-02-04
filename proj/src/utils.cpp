#include "utils.h"
#include "log_colormap.h"
void FlowImage::readFlowField(const std::string fileName)
{
	Mat flowRaw = imread(fileName, IMREAD_UNCHANGED);

	if (flowRaw.data == NULL) {
		cout << "can't open flow file " << fileName << endl;
		exit(1);
	}
	else {

		int nr = flowRaw.rows; // number of rows  
		int sc = flowRaw.channels();
		if (sc < 3)
		{
			cout << "not a valid KITTI format flow file " << endl;
			exit(1);
		}

		int nc = flowRaw.cols; // total number of elements per line
		width_ = nc;
		height_ = nr;
		data_ = (float*)malloc(width_*height_ * 3 * sizeof(float));

		if (flowRaw.isContinuous()) {
			nc = nc*nr;
			nr = 1;
		}

		for (int v = 0; v < nr; v++) {
			ushort* src = flowRaw.ptr<ushort>(v);

			for (int u = 0; u < nc; u++) {
				if (src[u*sc]) {
					setFlowU(u, v, (src[u*sc + 2] - 32768.0f) / 64.0f);
					setFlowV(u, v, (src[u*sc + 1] - 32768.0f) / 64.0f);
					setValid(u, v, true);
				}
				else {
					setFlowU(u, v, 0);
					setFlowV(u, v, 0);
					setValid(u, v, false);
				}
			}
		}
	}
}

void FlowImage::writeFlowField(const std::string fileName) {

	Mat flowRaw(height_, width_, CV_16UC3);
	int tc = flowRaw.channels();

	for (int32_t v = 0; v<height_; v++) {
		ushort* tgt = flowRaw.ptr<ushort>(v);
		for (int32_t u = 0; u<width_; u++) {
			if (isValid(u, v)) {
				tgt[tc*u + 2] = (uint16_t)std::max(std::min(getFlowU(u, v)*64.0f + 32768.0f, 65535.0f), 0.0f);
				tgt[tc*u + 1] = (uint16_t)std::max(std::min(getFlowV(u, v)*64.0f + 32768.0f, 65535.0f), 0.0f);
				tgt[tc*u]	  = 1;
			}
			else {
				tgt[tc*u + 2] = 0;
				tgt[tc*u + 1] = 0;
				tgt[tc*u] = 0;
			}

		}
	}

	imwrite(fileName, flowRaw);
}

// compute error map of flow field, given the non-occluded and occluded
// ground truth optical flow maps. stores result as color png image.
Mat FlowImage::errorImage(FlowImage &F_noc, FlowImage &F_occ, bool log_colors) {
	Mat image(height(), width(), CV_8UC3);
	for (int32_t v = 1; v<height() - 1; v++) {
		uchar* tgt = image.ptr<uchar>(v);
		for (int32_t u = 1; u<width() - 1; u++) {
			if (F_occ.isValid(u, v)) {
	
				if (log_colors) {
					float dfu = getFlowU(u, v) - F_occ.getFlowU(u, v);
					float dfv = getFlowV(u, v) - F_occ.getFlowV(u, v);
					float f_err = sqrt(dfu*dfu + dfv*dfv);
					float f_mag = F_occ.getFlowMagnitude(u, v);
					float n_err = std::min(f_err / 3.0, 20.0*f_err / f_mag);
					for (int32_t i = 0; i<10; i++) {
						if (n_err >= LC[i][0] && n_err<LC[i][1]) {
							tgt[3 * u + 2] = (uint8_t)LC[i][2];
							tgt[3 * u + 1] = (uint8_t)LC[i][3];
							tgt[3 * u] = (uint8_t)LC[i][4];
						}
					}
					if (!F_noc.isValid(u, v)) {
						tgt[3 * u + 2] *= 0.5;
						tgt[3 * u + 1] *= 0.5;
						tgt[3 * u] *= 0.5;
					}
				}
				else {
					float dfu = getFlowU(u, v) - F_occ.getFlowU(u, v);
					float dfv = getFlowV(u, v) - F_occ.getFlowV(u, v);
					float f_err = std::min<float>(sqrt(dfu*dfu + dfv*dfv), 5.0) / 5.0;
					tgt[3 * u + 2] = (uint8_t)(f_err*255.0);
					tgt[3 * u + 1] = (uint8_t)(f_err*255.0);
					tgt[3 * u] = (uint8_t)(f_err*255.0);
					if (!F_noc.isValid(u, v)) {
						tgt[3 * u + 1] = 0;
						tgt[3 * u] = 0;
					}
				}
				for (int32_t v2 = v - 1; v2 <= v + 1; v2++) {
					uchar* tgt2 = image.ptr<uchar>(v2);
					for (int32_t u2 = u - 1; u2 <= u + 1; u2++) {
						tgt2[3 * u2 + 2] = tgt[3 * u + 2];
						tgt2[3 * u2 + 1] = tgt[3 * u + 1];
						tgt2[3 * u2] = tgt[3 * u];
					}
				}
			}
		}
	}
	return image;
}

Mat read_calib_file(string fileName, bool isKITTI2015)
{
	Mat P(3, 4, CV_32F);
	ifstream f(fileName);
	if (!f.is_open())
	{
		cout << "can't open calibration file " << fileName << endl;
		exit(1);
	}
	if (!isKITTI2015) {
		string tmp;
		f >> tmp;
		for (int j = 0; j < 3; j++)
		{
			float_t* data = P.ptr<float_t>(j);
			f >> data[0] >> data[1] >> data[2] >> data[3];

		}
	}
	else {
		string line;
		getline(f, line);
		getline(f, line);
		for (int i = 0; i < 7; i++) 
			getline(f, line);

		string tmp;
		f >> tmp;

		for (int j = 0; j < 3; j++)
		{
			float_t* data = P.ptr<float_t>(j);
			f >> data[0] >> data[1] >> data[2] >> data[3];

		}

	}
	f.close();
	//cout << P << endl;
	return P;
}