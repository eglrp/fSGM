#include "mex.h"
#include <nmmintrin.h>
/*
 * calc_cost_pyd_sgm.c 
 * Perfrom cost volume construct and sgm for pyramidal sgm OF method. 
 * 
 * Description: 
 * Construct a cost volume for image 1/image 2 with fixed search region. (search center could be initialized by a given 
 * mv map, and then perform sgm on the generated Cost volume
 * the output of this program is a 2D index map, each element is the corresponding best cost index for each pixel 
 * 
 *
 * The calling syntax is:
 *
 *      [C, minIdx, minC, mvSub] = calc_cost_sgm(I1, I2, preMv, halfSearchWinSize, aggSize, subPixelRefine)
 *     
 * Input:
 * I1/I2 are input images
 * preMv is mv map from previous pyramidal level, must be have same or large size with I1/I2
 * halfSearchWinSize is the half search windows size in vertical direction. it is doubled in horizontal direction
 * aggSize is the aggregation window size. typically 5x5 is good
 * subPixelRefine: enable sub-pixel position calculation. if set mvSub contains the subpixel location of current level.  
 *
 * Output:
 * C: generated Cost volume
 * minIdx: output index (disparity) by sgm
 * minC: the corresponding sum of the path cost w.r.t. minIdx
 * mvSub: subpixel localtion
*/
typedef unsigned char PathCost;
typedef unsigned char PixelType;
typedef unsigned char CostType;

#define MAX_PATH_COST 255

typedef struct _costEntry
{
	int mvx;
	int mvy;
	int cost;
} CostEntry;

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

inline void sgm_step(CostEntry* L, //current path cost
					 CostEntry* Lpre, //previous path cost
					 CostEntry* C, //cost map
					 int dMax, 
					 int P1, int P2)
{
	PathCost minPathCost = MAX_PATH_COST;
	PathCost LpreMin = Lpre[dMax].cost; //get minimum value of pre path cost

	for (int d = 0; d< dMax; d++) {
		int mvx = C[d].mvx;
		int mvy = C[d].mvy;

		PathCost min1 = LpreMin + P2;
		PathCost min2 = LpreMin + P2;
		PathCost min3 = LpreMin + P2;
		PathCost bestCost = min3;

		for (int d2 = 0; d2 < dMax; d2++) {

			int mvxPre = Lpre[d2].mvx;
			int mvyPre = Lpre[d2].mvy;

			if(mvx == mvxPre && mvy == mvyPre) // ||d-d'|| = 0
				min1 = Lpre[d2].cost;
			else if(abs(mvx - mvxPre) <= 2 && abs(mvy - mvyPre) <= 2) {// ||d-d'|| < r
				min2 = min(min2, Lpre[d2].cost + P1);
			} 
		}

		bestCost = min(bestCost, min1);
		bestCost = min(bestCost, min2);
		L[d].cost = (C[d].cost + bestCost) - LpreMin;
		L[d].mvx = mvx;
		L[d].mvy = mvy;
		minPathCost = min(L[d].cost, minPathCost);
	}

	L[dMax].cost = minPathCost; //set minimum value of current path cost
}


inline int adaptive_P2(int P2, int pixCur, int pixPre) {
    const int threshold = 50;
    
    return (abs(pixCur - pixPre) > threshold ? P2 / 8 : P2);
}

/* sgm on 3-D cost volume
 * Output:
 * minC is the corresponding cost along with best index
 * mvSub is the output subpixel position for mvx/mvy
 *
 * Input:
 * C: 3-d cost volume
 * width/height/dMax: width/height/dMax(third dimension) of C
 * mvPre: previous level's the mv map
 * mvWidth/mvHeight: width/height of mvPre
 * P1/P2: small/large penalty
 *
 */
  
void sgm2d(unsigned* minC, double* mvSub, 
        PixelType* I1, CostEntry* C, int width, int height, int dMax,
        int P1, int P2, int subpixelRefine )
{
    //allocate path cost buffers. dMax cost entries + 1 minimun cost entry
    CostEntry* L1 = (CostEntry*) mxMalloc (sizeof(CostEntry) * 2 * (dMax + 1));            //Left -> Right direction
    CostEntry* L2 = (CostEntry*) mxMalloc (sizeof(CostEntry) * 2 * width * (dMax + 1));  //top-left -> bottom right direction
    CostEntry* L3 = (CostEntry*) mxMalloc (sizeof(CostEntry) * 2 * width * (dMax + 1));    //up -> bottom direction
    CostEntry* L4 = (CostEntry*) mxMalloc (sizeof(CostEntry) * 2 * width * (dMax + 1));  //top-right->bottom left direction
    unsigned * Sp = (unsigned *) mxMalloc (sizeof(unsigned ) * width * height * dMax); //sum of path cost from all directions
    memset(Sp, 0, sizeof(unsigned)*width*height*dMax);

	double* flowX = mvSub;
	double* flowY = mvSub + width*height;

    const int pathCostEntryPerPixel = (dMax + 1); //dMax + 1 minimun
    const int pathCostEntryPerRow = width * pathCostEntryPerPixel;
    const int costPerRowEntry = width*dMax;
    
    const bool adpativeP2 = false;
    const int totalPass = 2;
    const bool enableDiagnalPath = false;

    int ystart = 0;
    int yend = height ;
    int ystep = 1;

    int xstart = 0;
    int xend = width ;
    int xstep = 1;

    for (int pass = 0; pass < totalPass; pass++) {
        if (pass == 1) {
            ystart = height - 1;
            yend = -1;
            ystep = -1;

            xstart = width - 1;
            xend = -1;
            xstep = -1;
        }

        CostEntry* ptrL1Pre = L1;
        CostEntry* ptrL1Cur = L1 + dMax + 1;
        CostEntry* ptrL3PreRow = L3;
        CostEntry* ptrL3CurRow = L3 + pathCostEntryPerRow;
        CostEntry* ptrL2PreRow = L2;
        CostEntry* ptrL2CurRow = L2 + pathCostEntryPerRow;
        CostEntry* ptrL4PreRow = L4;
        CostEntry* ptrL4CurRow = L4 + pathCostEntryPerRow;

        for (int y = ystart; y != yend; y += ystep) {

            unsigned* ptrSp = Sp + y*costPerRowEntry;
            CostEntry* ptrC = C + y*costPerRowEntry;

            for (int x = xstart; x != xend; x += xstep) {

                CostEntry* ptrL3Cur = ptrL3CurRow + x*(dMax + 1);
                CostEntry* ptrL3Pre = ptrL3PreRow + x*(dMax + 1);

                CostEntry* ptrL2Cur = ptrL2CurRow + x*(dMax + 1);
                CostEntry* ptrL2Pre = ptrL2PreRow + (x - xstep)*(dMax + 1);

                CostEntry* ptrL4Cur = ptrL4CurRow + x*(dMax + 1);
                CostEntry* ptrL4Pre = ptrL4PreRow + (x + xstep)*(dMax + 1);

                CostEntry* ptrCCur = ptrC + x*dMax;

                if (x == xstart) {
                    memcpy(ptrL1Cur, ptrCCur, sizeof(CostEntry)*dMax);
                    ptrL1Cur[dMax].cost = 0;

                    if (enableDiagnalPath) {
                        memcpy(ptrL2Cur, ptrCCur, sizeof(CostEntry)*dMax);
                        ptrL2Cur[dMax].cost = 0;
                    }
                }

                if (y == ystart) {
                    memcpy(ptrL3Cur, ptrCCur, sizeof(CostEntry)*dMax);
                    ptrL3Cur[dMax].cost = 0;

                    if (enableDiagnalPath) {
                        memcpy(ptrL2Cur, ptrCCur, sizeof(CostEntry)*dMax);
                        ptrL2Cur[dMax].cost = 0;

                        memcpy(ptrL4Cur, ptrCCur, sizeof(CostEntry)*dMax);
                        ptrL4Cur[dMax].cost = 0;
                    }
                }

                if (x == xend) {
                    if (enableDiagnalPath) {
                        memcpy(ptrL4Cur, ptrCCur, sizeof(CostEntry)*dMax);
                        ptrL4Cur[dMax].cost = 0;
                    }
                }

                if (x != xstart) {
                    //hint map may have different size with image, must set width to mvWidth, otherwise will have 45degree error propagation issue 
                    //when 2nd pyd processing
                    PixelType pixCur = I1[width*y + x];
                    PixelType pixPre = I1[width*y + x - xstep];
                    
                    sgm_step(ptrL1Cur,              //current path cost
                        ptrL1Pre,                   //previous path cost
                        ptrCCur,                    //cost map
                        dMax, P1, adpativeP2 ? adaptive_P2(P2, pixCur, pixPre) : P2);
                }


                if (y != ystart) {
                    PixelType pixCur = I1[width*y + x];
                    PixelType pixPre = I1[width*(y-ystep) + x];

                    sgm_step(ptrL3Cur,              //current path cost
                        ptrL3Pre,                   //previous path cost
                        ptrCCur,                    //cost map
                        dMax, P1, adpativeP2 ? adaptive_P2(P2, pixCur, pixPre) : P2);
                }

                if (enableDiagnalPath) {
                    if (x != xstart && y != ystart) {

                        PixelType pixCur = I1[width*y + x];
                        PixelType pixPre = I1[width*(y - ystep) + x - xstep];

                        sgm_step(ptrL2Cur,          //current path cost
                            ptrL2Pre,               //previous path cost
                            ptrCCur,                //cost map
                            dMax, P1, adpativeP2 ? adaptive_P2(P2, pixCur, pixPre) : P2);
                    }

                    if (x != xend && y != ystart) {

                        PixelType pixCur = I1[width*y + x];
                        PixelType pixPre = I1[width*(y - ystep) + x + xstep];

                        sgm_step(ptrL4Cur,          //current path cost
                            ptrL4Pre,               //previous path cost
                            ptrCCur,                //cost map
                            dMax, P1, adpativeP2 ? adaptive_P2(P2, pixCur, pixPre) : P2);

                    }
                }

                for (int d = 0; d < dMax; d++) {
                    ptrSp[x*dMax + d] += ptrL1Cur[d].cost + ptrL3Cur[d].cost;
                    if (enableDiagnalPath) {
                        ptrSp[x*dMax + d] += ptrL2Cur[d].cost + ptrL4Cur[d].cost;
                    }
                }

                //swap buffer pointer for left->right direction
                CostEntry* tmp = ptrL1Pre;
                ptrL1Pre = ptrL1Cur;
                ptrL1Cur = tmp;
            }

            //swap buffer pointer for top->bottom direction
            CostEntry* tmp = ptrL3PreRow;
            ptrL3PreRow = ptrL3CurRow;
            ptrL3CurRow = tmp;

            if (enableDiagnalPath) {
                //swap buffer pointer for top left->bottom rightdirection
                tmp = ptrL2PreRow;
                ptrL2PreRow = ptrL2CurRow;
                ptrL2CurRow = tmp;

                //swap buffer pointer for top right->bottom leftdirection
                tmp = ptrL4PreRow;
                ptrL4PreRow = ptrL4CurRow;
                ptrL4CurRow = tmp;
            }
        }
    }
    
    //mexPrintf("path cost aggregate done!\n");
    for(int y = 0; y< height; y++) {
        unsigned* SpPtr = Sp + y*costPerRowEntry;
		
        for (int x = 0; x <width; x++) {
            CostEntry* ptrC = C + width* dMax*y + dMax*x;
            unsigned minCost = SpPtr[x*dMax];
            unsigned minIdx = 0;
            for (int d = 1; d<dMax; d++) {
                
                if(SpPtr[x*dMax + d] < minCost) {
                    minCost = SpPtr[x*dMax + d] ;
				//if(ptrC[d].cost < minCost) {
				//	minCost = ptrC[d].cost;
                    minIdx = d;
                }
            }
            minC[y*width +x] = minCost;
            flowX[y*width + x] = ptrC[minIdx].mvx;
			flowY[y*width + x] = ptrC[minIdx].mvy;
        }
    }
		

    mxFree(L1);
    mxFree(L2);
    mxFree(L3);
    mxFree(L4);
    mxFree(Sp);
}
        
void calc_cost(CostEntry* C,
	const unsigned* cen1, const unsigned* cen2, int width, int height,
	const double* preMv, int mvWidth, int mvHeight,
	int winRadiusAgg, int winRadiusX, int winRadiusY, 
	int hintXradius,  int hintYradius, int dMax)
{
	double* pMvx = preMv;
	double* pMvy = preMv + mvWidth*mvHeight;
	int winPixels = (2 * winRadiusAgg + 1)*(2 * winRadiusAgg + 1);
	const CostType defaultCost = 5;

	const int step = 8;
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			int hintIdx = 0;

			CostEntry* ptrC = C + dMax*width*y + dMax*x;

			int d = 0;

			for (int dy = -step*hintYradius; dy <= step*hintYradius; dy += step) {
				for (int dx = -step*hintXradius; dx <= step*hintXradius; dx += step) {
					int yn = min(mvHeight-1, max(0, y + dy));
					int xn = min(mvWidth-1,  max(0, x + dx));

					double mvx = pMvx[mvWidth*yn + xn];
					double mvy = pMvy[mvWidth*yn + xn];
					
					
					for (int offx = -winRadiusX; offx <= winRadiusX; offx++) {
						for (int offy = -winRadiusY; offy <= winRadiusY; offy++) {
				
							unsigned costSum = 0;

							for (int aggy = -winRadiusAgg; aggy <= winRadiusAgg; aggy++) {
								for (int aggx = -winRadiusAgg; aggx <= winRadiusAgg; aggx++) {

									int y1 = y + aggy;
									int x1 = x + aggx;

									if(y1 < 0 || y1 > height-1 || x1 <0 || x1 > width-1) {
										costSum += defaultCost;  //add constant cost if not valid current pixel position
										continue;
									}

									unsigned cenCode1 = cen1[width*y1 + x1];

									int y2 = (offy + y1) + mvy ;
									int x2 = (offx + x1) + mvx ;

									if(y2 < 0 || y2 > height-1 || x2 <0 || x2 > width-1) {
										costSum += defaultCost; //add constant cost if not valid reference pixel position
										continue;
									}

									unsigned cenCode2 = cen2[width*y2 + x2];

									int censusCost = _mm_popcnt_u32((cenCode1^cenCode2));
									costSum += censusCost;
								}
							}
							
							ptrC[d].cost = (1.0 * costSum / winPixels) + 0.5;
							ptrC[d].mvx = mvx + offx;
							ptrC[d].mvy = mvy + offy;
							
							d++;
						}
					}
					hintIdx++;
				} //dx
			} //dy

			mxAssert(d == dMax, "incorrect candidates per pixel!\n");
		}
	}
}
/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    PixelType *I1;             /* pointer to Input image I1 */
    PixelType *I2;             /* pointer to Input image I2 */
    double *preMv;          /* pointer to initial search position */
    
    mwSize width;               /* rows (width) assume I1/I2 are permuted before passing in */
    mwSize height;               /* cols (height)*/
    
    I1 = (PixelType*)mxGetData(prhs[0]);
    I2 = (PixelType*)mxGetData(prhs[1]); 
    
    preMv = mxGetPr(prhs[2]);
    width = mxGetM(prhs[0]);
    height = mxGetN(prhs[0]);
    
    double halfSearchWinSize = mxGetScalar(prhs[3]);
    double aggSize= mxGetScalar(prhs[4]);
    int subPixelRefine = mxGetScalar(prhs[5]);

    int P1 = mxGetScalar(prhs[6]);
    int P2 = mxGetScalar(prhs[7]);
    
    /* create the output matrix */

    const mwSize dims2[] = { width, height };       //output: best index map
    const mwSize dims3[] = { width, height, 2 };    //output: subpixel position
    plhs[0] = mxCreateNumericArray(2, dims2, mxUINT32_CLASS, mxREAL);
	plhs[1] = mxCreateNumericArray(3, dims3, mxDOUBLE_CLASS, mxREAL);

    unsigned* minC = (unsigned*)mxGetData(plhs[0]);
    double* flowResult = mxGetPr(plhs[1]);

    unsigned* cen1 = (unsigned*)mxMalloc(width * height * sizeof(unsigned));
    unsigned *cen2 = (unsigned*)mxMalloc(width * height * sizeof(unsigned));

    census(I1, cen1, width, height, 2);
    census(I2, cen2, width, height, 2);
    
    int winRadiusY = halfSearchWinSize;
    int winRadiusX = halfSearchWinSize;
    int winRadiusAgg = (int)aggSize/2;

	const int hintYradius = 1;
	const int hintXradius = 1;
	const int candidatePerHint = (2*winRadiusX+1)*(2*winRadiusY+1);
	const int hintN = (2*hintXradius+1)*(2*hintYradius+1);

	int dMax = hintN*candidatePerHint ;

    mexPrintf("width: %d, height: %d, dMax: %d, winRadiusAgg: %d\n", width, height, dMax, winRadiusAgg);
    
    int mvWidth = mxGetM(prhs[2]);
    int mvHeight = mxGetN(prhs[2])/2;
    
	CostEntry* C2 = (CostEntry*) mxMalloc(width * height * sizeof(CostEntry) * dMax);
	//construct cost volume
	calc_cost(C2, cen1, cen2,  width, height, preMv, mvWidth, mvHeight,
		winRadiusAgg, winRadiusX, winRadiusY,
		hintXradius, hintYradius, dMax);

	//perform sgm
	sgm2d(minC, flowResult, 
		I1, C2, width, height, dMax, 
		P1, P2, subPixelRefine);

	mxFree(C2);
    mxFree(cen1);
    mxFree(cen2);
    
}