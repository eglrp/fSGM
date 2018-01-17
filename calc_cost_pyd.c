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
typedef unsigned char PathCost;
typedef unsigned char PixelType;
#define MAX_PATH_COST 255

void census(PixelType* img, unsigned * cen, int width, int height, int halfWin)
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
					if (img[x2 + width*y2] >= centerValue)
                        censusCode += 1;
					censusCode = censusCode << 1;
				}
			}
			cen[x + y*width] = censusCode;
        }
    }
}

//perform a single step to calculate path cost for current pixel position
inline PathCost sgm_step(PathCost* L, //current path cost
	PathCost* Lpre, //previous path cost
	PathCost LpreMin,
	PathCost* C, //cost map
    double dx, double dy, int searchWinX, int searchWinY, 
    int P1, int P2)
{
	PathCost minPathCost = MAX_PATH_COST;

    for (int sx = 0; sx < searchWinX; sx ++) {
        for (int sy = 0; sy < searchWinY; sy ++) {

            int ypre = sy + dy  + 0.5;
            int xpre = sx + dx  + 0.5;

			PathCost min1 = LpreMin + P2;
			PathCost min2 = LpreMin + P2;
			PathCost min3 = LpreMin + P2;
			PathCost bestCost = min3;

            // ||d-d'|| = 0
            if(xpre >= 0 && xpre < searchWinX && ypre>=0 && ypre <searchWinY) {
                int dtemp = xpre * searchWinY + ypre;
                min1 = Lpre[dtemp];
            }
            // ||d-d'|| < r
			const int r = 2;
            for(int k = -r; k<= r; k++) {
                for(int m = -r; m<= r; m++) {
                    if(m == 0 && k == 0)
                        continue;

                    int ty = ypre + k;
                    int tx = xpre + m;

                    if(tx >= 0 && tx < searchWinX && ty >= 0 && ty < searchWinY)
                    {
                        int dtemp = tx * searchWinY + ty;
                        min2 = min(min2, Lpre[dtemp] + P1);
                    }
                }
            }

			bestCost = min(bestCost, min1);
			bestCost = min(bestCost, min2);

            int d = sx * searchWinY + sy;
			mxAssert(bestCost > LpreMin, "bestCost Must > LpreMin\n");
            L[d] = C[d] + bestCost - LpreMin;
			minPathCost = min(L[d], minPathCost);
        }
    }

    return minPathCost;
}
/* sgm on 3-D cost volume
 * Output:
 * bestD is the output best index along the third dimension
 * minC is the corresponding cost along with best index
 * mvSub is the output subpixel position for mvx/mvy
 *
 * Input:
 * C: 3-d cost volume
 * width/height/dMax: width/height/dMax(third dimension) of C
 * mvPre: previous level's the mv map
 * mvWidth/mvHeight: width/height of mvPre
 * searchWinX: search window size at x-direction
 * searchWinY: search window size at y-direction
 * P1/P2: small/large penalty
 * subpixelRefine: enable/disable subpixel position estimation
 *
 */
void sgm2d(unsigned* bestD, unsigned* minC, double* mvSub, 
		PixelType* I1, unsigned char* C, int width, int height, int dMax,
        double* mvPre, int mvWidth, int mvHeight, 
        int searchWinX, int searchWinY, int P1, int P2, int subpixelRefine )
{
    
	//allocate path cost buffers
    PathCost* L1 = (PathCost*) mxMalloc (sizeof(PathCost) * 2 * dMax);			//Left -> Right direction
    PathCost* L2 = (PathCost*) mxMalloc (sizeof(PathCost) * 2 * width * dMax);  //top-left -> bottom right direction
	PathCost* L3 = (PathCost*) mxMalloc (sizeof(PathCost) * 2 * width * dMax);	//up -> bottom direction
	PathCost* L4 = (PathCost*) mxMalloc (sizeof(PathCost) * 2 * width * dMax);  //top-right->bottom left direction
    unsigned * Sp = (unsigned *) mxMalloc (sizeof(unsigned ) * width * height * dMax); //sum of path cost from all directions
	memset(Sp, 0, sizeof(unsigned)*width*height*dMax);
    PathCost minL1;
    
	PathCost* minL3 = (unsigned char*) mxMalloc(sizeof(unsigned char)* width);
    
    double* pMvx = mvPre;
    double* pMvy = mvPre + mvWidth * mvHeight;
    
    const int pathCostPerRowEntry = width * dMax;
	bool adpativeP2 = false;
	int ystart = 0;
	int yend = height ;
	int ystep = 1;

	int xstart = 0;
	int xend = width ;
	int xstep = 1;

	for (int pass = 0; pass <= 1; pass++) {
		if (pass == 1) {
			ystart = height - 1;
			yend = -1;
			ystep = -1;

			xstart = width - 1;
			xend = -1;
			xstep = -1;
		}

		PathCost* ptrL1Pre = L1;
		PathCost* ptrL1Cur = L1 + dMax;
		PathCost* ptrL3Pre = L3;
		PathCost* ptrL3Cur = L3 + pathCostPerRowEntry;

		for (int y = ystart; y != yend; y += ystep) {

			unsigned* Spptr = Sp + y*pathCostPerRowEntry;
			unsigned char* Cptr = C + y*pathCostPerRowEntry;

			for (int x = xstart; x != xend; x += xstep) {

				if (x == xstart) {
					minL1 = MAX_PATH_COST;
					memcpy(ptrL1Cur, Cptr + x*dMax, sizeof(PathCost)*dMax);
				}

				if (y == ystart) {
					minL3[x] = MAX_PATH_COST;
					memcpy(ptrL3Cur, Cptr + x*dMax, sizeof(PathCost)*dMax);
				}


				if (x != xstart) {
					//hint map may have different size with image, must set width to mvWidth, otherwise will have 45degree error propagation issue 
					//when 2nd pyd processing
					double dx = pMvx[y*mvWidth + x] - pMvx[y*mvWidth + x - xstep];
					double dy = pMvy[y*mvWidth + x] - pMvy[y*mvWidth + x - xstep];
					PixelType pixCur = I1[width*y + x];
					PixelType pixPre = I1[width*y + x - xstep];
					
					minL1 = sgm_step(ptrL1Cur,			//current path cost
						ptrL1Pre,						//previous path cost
						minL1,
						Cptr + x*dMax,					//cost map
						dx, dy, searchWinX, searchWinY, P1, adpativeP2 ? (abs((int)pixCur - (int)pixPre)> 30 ? P2 / 3 : P2) : P2);
				}


				if (y != ystart) {
					double dx = pMvx[y*mvWidth + x] - pMvx[(y - ystep)*mvWidth + x];
					double dy = pMvy[y*mvWidth + x] - pMvy[(y - ystep)*mvWidth + x];

					PixelType pixCur = I1[width*y + x];
					PixelType pixPre = I1[width*(y-ystep) + x];

					minL3[x] = sgm_step(ptrL3Cur + x*dMax,//current path cost
						ptrL3Pre + x*dMax,				//previous path cost
						minL3[x],
						Cptr + x*dMax,					//cost map
						dx, dy, searchWinX, searchWinY, P1, adpativeP2 ? (abs((int)pixCur - (int)pixPre)> 30 ? P2 / 3 : P2) : P2);
				}


				for (int d = 0; d < dMax; d++) {
					Spptr[x*dMax + d] +=  ptrL1Cur[d] + ptrL3Cur[x*dMax + d];
				}

				//swap buffer pointer for left->right direction
				PathCost* tmp = ptrL1Pre;
				ptrL1Pre = ptrL1Cur;
				ptrL1Cur = tmp;

			}

			//swap buffer pointer for top->bottom direction
			PathCost* tmp = ptrL3Pre;
			ptrL3Pre = ptrL3Cur;
			ptrL3Cur = tmp;
		}
	}
    
    for(int y = 0; y< height; y++) {
		unsigned* SpPtr = Sp + y*pathCostPerRowEntry;
        for (int x = 0; x <width; x++) {
            
            unsigned minCost = SpPtr[x*dMax];
            unsigned minIdx = 0;
            for (int d = 1; d<dMax; d++) {
                
                if(SpPtr[x*dMax + d] < minCost) {
                    minCost = SpPtr[x*dMax + d];
                    minIdx = d;
                }
            }
            minC[y*width +x] = minCost;
            bestD[y*width +x] = minIdx;
        }
    }
    
    double * ptrMvSubMvx = mvSub;
    double * ptrMvSubMvy = mvSub + width*height;
    
    if(subpixelRefine) {
        /*
         * do subpixel quadratic interpolation:
         *fit parabola into (x1=d-1, y1=C[d-1]), (x2=d, y2=C[d]), (x3=d+1, y3=C[d+1])
         *then find minimum of the parabola
         */
        for(int y = 0; y< height; y++) {
            unsigned* SpPtr = Sp + y*pathCostPerRowEntry;
            
            for (int x = 0; x <width; x++) {
                unsigned bestIdx = bestD[y*width + x];
                
                double c0 = (double) (SpPtr[x*dMax + bestIdx]);
                
                int dx = bestIdx / searchWinY;
                int dy = bestIdx % searchWinY;
                
                if(dy > 0 && dy < searchWinY - 1) {
                    double cLeft = (double)(SpPtr[x*dMax + bestIdx - 1]);
                    double cRight  = (double)(SpPtr[x*dMax + bestIdx + 1]);
                    
                    if (cRight < cLeft)
                        ptrMvSubMvy[y*width + x] = (cRight-cLeft)/(c0 - cLeft)/2.0;
                    else
                        ptrMvSubMvy[y*width + x] = (cRight-cLeft)/(c0 - cRight)/2.0;
                    
                } else {
                    ptrMvSubMvy[y*width + x] = 0;
                }
                
                if(dx >0 && dx < searchWinX-1) {
                    double cLeft = (double)(SpPtr[x*dMax + bestIdx - searchWinY]);
                    double cRight  = (double)(SpPtr[x*dMax + bestIdx + searchWinY]);
                    
                    if (cRight < cLeft)
                        ptrMvSubMvx[y*width + x] = (cRight-cLeft)/(c0 - cLeft)/2.0;
                    else
                        ptrMvSubMvx[y*width + x] = (cRight-cLeft)/(c0 - cRight)/2.0;
                    
                } else {
                    ptrMvSubMvx[y*width + x] = 0;
                }
            }
        }
        
    }
       
       
    mxFree(L1);
    mxFree(L2);
    mxFree(L3);
    mxFree(L4);
    mxFree(Sp);
    mxFree(minL3);
}
        
void calc_cost(unsigned char* C, 
	const unsigned* cen1, const unsigned* cen2, int width, int height,
	const double* preMv, int mvWidth, int mvHeight, 
	int winRadiusAgg, int winRadiusX, int winRadiusY)
{
	double* pMvx = preMv;
	double* pMvy = preMv + mvWidth*mvHeight;
	int winPixels = (2 * winRadiusAgg + 1)*(2 * winRadiusAgg + 1);
    
	for (int offy = -winRadiusY; offy <= winRadiusY; offy++) {
		for (int offx = -winRadiusX; offx <= winRadiusX; offx++) {
			int d = (offx + winRadiusX)* (2 * winRadiusY + 1) + offy + winRadiusY;
			//mexPrintf("d: %d\n", d);
			unsigned char* pC = C + d * width * height;

			for (int y = 0; y< height; y++) {
				for (int x = 0; x< width; x++) {

					unsigned costSum = 0;
					double mvx = pMvx[mvWidth*y + x];
					double mvy = pMvy[mvWidth*y + x];

					for (int aggy = -winRadiusAgg; aggy <= winRadiusAgg; aggy++) {
						for (int aggx = -winRadiusAgg; aggx <= winRadiusAgg; aggx++) {

							int y1 = y + aggy;
							int x1 = x + aggx;

                            if(y1 < 0 || y1 > height-1 || x1 <0 || x1 > width-1) {
                                costSum += 5;  //add constant cost if not valid current pixel position
                                continue;
                            }
                            
							//y1 = y1 < 0 ? 0 : (y1 > height - 1 ? height - 1 : y1);
							//x1 = x1 < 0 ? 0 : (x1 > width - 1 ? width - 1 : x1);

							unsigned cenCode1 = cen1[width*y1 + x1];

							int y2 = 1.0*(offy + y1) + mvy + 0.5;
							int x2 = 1.0*(offx + x1) + mvx + 0.5;
                            
                            if(y2 < 0 || y2 > height-1 || x2 <0 || x2 > width-1) {
                                costSum += 5; //add constant cost if not valid reference pixel position
                                continue;
                            }
							//y2 = y2 < 0 ? 0 : (y2 > height - 1 ? height - 1 : y2);
							//x2 = x2 < 0 ? 0 : (x2 > width - 1 ? width - 1 : x2);

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
   
    int subPixelRefine = mxGetScalar(prhs[5]);
    
    /* create the output matrix */
    const mwSize dims[]={width, height, dMax};		//output: costvolume
	const mwSize dims2[] = { width, height };		//output: best index map
	const mwSize dims3[] = { width, height, 2 };	//output: subpixel position

    plhs[0] = mxCreateNumericArray(3, dims, mxUINT8_CLASS, mxREAL);
    plhs[1] = mxCreateNumericArray(2, dims2, mxUINT32_CLASS, mxREAL);
    plhs[2] = mxCreateNumericArray(2, dims2, mxUINT32_CLASS, mxREAL);
	plhs[3] = mxCreateNumericArray(3, dims3, mxDOUBLE_CLASS, mxREAL);
    unsigned char* C = (unsigned char*) mxGetData(plhs[0]);
    unsigned * bestD = (unsigned*) mxGetData(plhs[1]);
	unsigned* minC = (unsigned*)mxGetData(plhs[2]);
	double* mvSub = mxGetPr(plhs[3]);

    unsigned* cen1 = (unsigned*)mxMalloc(width * height * sizeof(unsigned));
    unsigned *cen2 = (unsigned*)mxMalloc(width * height * sizeof(unsigned));

    census(I1, cen1, width, height, 2);
    census(I2, cen2, width, height, 2);
    
    int winRadiusY = halfSearchWinSize;
    int winRadiusX = 2*halfSearchWinSize;
    int winRadiusAgg = (int)aggSize/2;
    mexPrintf("width: %d, height: %d, dMax: %d, winRadiusAgg: %d\n", width, height, dMax, winRadiusAgg);
    
    int mvWidth = mxGetM(prhs[2]);
    int mvHeight = mxGetN(prhs[2])/2;
    
    
    //construct cost volume
	calc_cost(C, cen1, cen2, width, height, preMv, mvWidth, mvHeight,
		winRadiusAgg, winRadiusX, winRadiusY);
    
    
    int P1 = 6;
    int P2 = 64;
    
    int totalElementNum = width*height*dMax;
    unsigned char* Cre = (unsigned char*) mxMalloc(totalElementNum*sizeof(unsigned char));
    for(int d = 0; d< dMax; d++) { 
        for (int y = 0; y< height; y++) {
            for(int x = 0; x <width; x++) {
                int indOrig = d*width*height + (y*width + x);
                int indTarget = y* width *dMax + x*dMax + d;
                mxAssert(indOrig < totalElementNum && indTarget < totalElementNum, "incorrect indexing");
                Cre[indTarget] = C[indOrig];
            }
        }
    }

	//perform sgm
    sgm2d(bestD, minC, mvSub, 
        I1, Cre, width, height, dMax, 
        preMv, mvWidth, mvHeight, 
        winRadiusX*2 +1, winRadiusY*2+1, P1,  P2, subPixelRefine);
     
    mxFree(Cre);
    mxFree(cen1);
    mxFree(cen2);
    
}