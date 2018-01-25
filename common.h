#ifndef _COMMON_H_
#define _COMMON_H_
#include <algorithm>
typedef unsigned char PathCost;
typedef unsigned char PixelType;
typedef unsigned char CostType;

#define MAX_PATH_COST 255
#define SUBPIXEL_PRECISION 8
void census(PixelType* img, unsigned * cen, int width, int height, int halfWin);

inline int clamp(int val, int minVal, int maxVal) {return std::min<int>(maxVal, std::max<int>(minVal, val));}
#endif