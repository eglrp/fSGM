#include "mex.h"
#include <nmmintrin.h>
/*
 * calc_cost_pyd.c 
 *
 * Construct a cost volume for image 1/image 2 with fixed search region. (search center could be initialized by a given 
 * mv map
 * 
 *
 * The calling syntax is:
 *
 *		C = calc_cost_pyd(I1, I2, preMv, halfSearchWinSize, aggSize)
 *
*/

void census(unsigned char* img, unsigned * cen, int width, int height, int halfWin)
{
    
    
    for (int y = 0; y< height; y++) {
    	for(int x= 0; x< width; x++) {
            
            unsigned censusCode = 0;
            unsigned char centerValue = img[x +  width*y];
            for (int offsetY = -halfWin; offsetY <= halfWin; ++offsetY) {
				for (int offsetX = -halfWin; offsetX <= halfWin; ++offsetX) {
                    int y2 = y + offsetY;
                    int x2 = x + offsetX;
                    
                    y2 = y2 < 0? 0 : (y2 > height-1? height-1:y2);
                    x2 = x2 < 0? 0 : (x2 > width-1? width-1:x2);
					if (img[x + offsetX + width*(y + offsetY)] >= centerValue)
                        censusCode += 1;
					censusCode = censusCode << 1;
				}
			}
			cen[x + y*width] = censusCode;
        }
    }
}


/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    unsigned char *I1;             /* pointer to Input image I1 */
    unsigned char *I2;             /* pointer to Input image I2 */
    double *preMv;          /* pointer to initial search position */
    
    mwSize width;               /* rows (width) assume I1/I2 are permuted before passing in */
    mwSize height;               /* cols (height)*/
    
    I1 = (unsigned char*)mxGetData(prhs[0]);
    I2 = (unsigned char*)mxGetData(prhs[1]);
    
    
    preMv = mxGetPr(prhs[2]);
    
    width = mxGetM(prhs[0]);
    height = mxGetN(prhs[0]);
    
    double halfSearchWinSize = mxGetScalar(prhs[3]);
    double aggSize= mxGetScalar(prhs[4]);
    int dMax = (4*halfSearchWinSize+1)*(2*halfSearchWinSize+1);
   
    
    /* create the output matrix */
    const mwSize dims[]={width, height, dMax};
    plhs[0] = mxCreateNumericArray(3, dims, mxUINT8_CLASS, mxREAL);
    unsigned char* C = (unsigned char*) mxGetData(plhs[0]);
    unsigned* cen1 = (unsigned*)malloc(width * height * sizeof(unsigned));
    unsigned *cen2 = (unsigned*)malloc(width * height * sizeof(unsigned));

    census(I1, cen1, width, height, 2);
    census(I2, cen2, width, height, 2);
    
    int winRadiusY = halfSearchWinSize;
    int winRadiusX = 2*halfSearchWinSize;
    int winRadiusAgg = (int)aggSize/2;
    //mexPrintf("width: %d, height: %d, dMax: %d, winRadiusAgg: %d\n", width, height, dMax, winRadiusAgg);
    
    int mvWidth = mxGetM(prhs[2]);
    int mvHeight = mxGetN(prhs[2])/2;
    
    double* pMvx = preMv;
    double* pMvy = preMv + mvWidth*mvHeight;
	int winPixels = (2 * winRadiusAgg + 1)*(2 * winRadiusAgg + 1);
    for (int offy = -winRadiusY; offy <= winRadiusY; offy++) {
        for(int offx = -winRadiusX; offx <= winRadiusX; offx ++) {
            int d = (offx + winRadiusX)* (2*winRadiusY + 1) + offy + winRadiusY;
            //mexPrintf("d: %d\n", d);
            unsigned char* pC = C + d * width * height;
            
            for (int y = 0; y< height; y++) {
                for(int x= 0; x< width; x++) {
                    
                    unsigned costSum = 0;
                    double mvx = pMvx[mvWidth*y + x];
                    double mvy = pMvy[mvWidth*y + x];
                            
                    for (int aggy = -winRadiusAgg; aggy <= winRadiusAgg; aggy ++) {
                        for (int aggx = -winRadiusAgg; aggx <=winRadiusAgg; aggx++){
                            
                            int y1 = y + aggy;
                            int x1 = x + aggx;
                            
                            y1 = y1 < 0 ? 0 : (y1 > height-1? height-1: y1);
                            x1 = x1 < 0 ? 0 : (x1 > width -1? width -1 :x1);

                            unsigned cenCode1 = cen1[width*y1 + x1];
                            
                            int y2 = 1.0*(offy + y1) + mvy + 0.5;
                            int x2 = 1.0*(offx + x1) + mvx + 0.5;
                            y2 = y2 < 0 ? 0 :(y2 > height-1? height -1:y2);
                            x2 = x2 < 0 ? 0 :(x2 > width -1? width -1 :x2);

                            unsigned cenCode2 = cen2[width*y2 + x2];
                            
                            int censusCost = _mm_popcnt_u32((cenCode1^cenCode2));
                            costSum += censusCost;
                        }
                    }
                    pC[y*width + x] = (1.0 * costSum / winPixels) + 0.5;
                }
            }
        }
    }
}