#include "mex.h"
#include <nmmintrin.h>
#include "common.h"
/*
 * calc_cost_sgm.cpp 
 * Perform cost volume construct and SGM for epipolar sgm OF method. 
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

    

//perform a single step to calculate path cost for current pixel position
inline void sgm_step(PathCost* L, //current path cost
    PathCost* Lpre, //previous path cost
    PathCost* C, //cost map
    int dMax, 
    int P1, int P2)
{
    PathCost minPathCost = MAX_PATH_COST;
    PathCost LpreMin = Lpre[dMax]; //get minimum value of pre path cost
    for (int d = 0; d < dMax; d ++) {

        //mxAssert((int)LpreMin + P2 < 256);
		// d= d'
        PathCost min1 = Lpre[d];

		//|d-d'| <= 1
        PathCost min2 = LpreMin + P2;
		if(d > 0) min2 = std::min<PathCost>(min2, Lpre[d-1] + P1);
		if(d < dMax-1) min2 = std::min<PathCost>(min2, Lpre[d+1] + P1);

		//|d-d'| >= 2
        PathCost min3 = LpreMin + P2;
        PathCost bestCost = min3;


        bestCost = min(bestCost, min1);
        bestCost = min(bestCost, min2);

        mxAssert(C[d] + bestCost >= LpreMin, "bestCost Must > LpreMin\n");

        L[d] = (C[d] + bestCost) - LpreMin;
        minPathCost = min(L[d], minPathCost);
        
    }

    L[dMax] = minPathCost; //set minimum value of current path cost
}

inline int adaptive_P2(int P2, int pixCur, int pixPre) {
    const int threshold = 50;
    
    return (abs(pixCur - pixPre) > threshold ? P2 / 8 : P2);
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
void sgm2d(unsigned* bestD, unsigned* minC,
        PixelType* I1, CostType* C, int width, int height, int dMax,
        int P1, int P2, int subpixelRefine )
{
    //allocate path cost buffers. dMax cost entries + 1 minimun cost entry
    PathCost* L1 = (PathCost*) mxMalloc (sizeof(PathCost) * 2 * (dMax + 1));            //Left -> Right direction
    PathCost* L2 = (PathCost*) mxMalloc (sizeof(PathCost) * 2 * width * (dMax + 1));  //top-left -> bottom right direction
    PathCost* L3 = (PathCost*) mxMalloc (sizeof(PathCost) * 2 * width * (dMax + 1));    //up -> bottom direction
    PathCost* L4 = (PathCost*) mxMalloc (sizeof(PathCost) * 2 * width * (dMax + 1));  //top-right->bottom left direction
    unsigned * Sp = (unsigned *) mxMalloc (sizeof(unsigned ) * width * height * dMax); //sum of path cost from all directions
    memset(Sp, 0, sizeof(unsigned)*width*height*dMax);

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

        PathCost* ptrL1Pre = L1;
        PathCost* ptrL1Cur = L1 + dMax + 1;
        PathCost* ptrL3PreRow = L3;
        PathCost* ptrL3CurRow = L3 + pathCostEntryPerRow;
        PathCost* ptrL2PreRow = L2;
        PathCost* ptrL2CurRow = L2 + pathCostEntryPerRow;
        PathCost* ptrL4PreRow = L4;
        PathCost* ptrL4CurRow = L4 + pathCostEntryPerRow;

        for (int y = ystart; y != yend; y += ystep) {

            unsigned* ptrSp = Sp + y*costPerRowEntry;
            CostType* ptrC = C + y*costPerRowEntry;

            for (int x = xstart; x != xend; x += xstep) {

                PathCost* ptrL3Cur = ptrL3CurRow + x*(dMax + 1);
                PathCost* ptrL3Pre = ptrL3PreRow + x*(dMax + 1);

                PathCost* ptrL2Cur = ptrL2CurRow + x*(dMax + 1);
                PathCost* ptrL2Pre = ptrL2PreRow + (x - xstep)*(dMax + 1);

                PathCost* ptrL4Cur = ptrL4CurRow + x*(dMax + 1);
                PathCost* ptrL4Pre = ptrL4PreRow + (x + xstep)*(dMax + 1);

                CostType* ptrCCur = ptrC + x*dMax;

                if (x == xstart) {
                    memcpy(ptrL1Cur, ptrCCur, sizeof(PathCost)*dMax);
                    ptrL1Cur[dMax] = 0;

                    if (enableDiagnalPath) {
                        memcpy(ptrL2Cur, ptrCCur, sizeof(PathCost)*dMax);
                        ptrL2Cur[dMax] = MAX_PATH_COST;
                    }
                }

                if (y == ystart) {
                    memcpy(ptrL3Cur, ptrCCur, sizeof(PathCost)*dMax);
                    ptrL3Cur[dMax] = 0;

                    if (enableDiagnalPath) {
                        memcpy(ptrL2Cur, ptrCCur, sizeof(PathCost)*dMax);
                        ptrL2Cur[dMax] = 0;

                        memcpy(ptrL4Cur, ptrCCur, sizeof(PathCost)*dMax);
                        ptrL4Cur[dMax] = 0;
                    }
                }

                if (x == xend) {
                    if (enableDiagnalPath) {
                        memcpy(ptrL4Cur, ptrCCur, sizeof(PathCost)*dMax);
                        ptrL4Cur[dMax] = 0;
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
                    ptrSp[x*dMax + d] += ptrL1Cur[d] + ptrL3Cur[d];
                    if (enableDiagnalPath) {
                        ptrSp[x*dMax + d] += ptrL2Cur[d] + ptrL4Cur[d];
                    }
                }

                //swap buffer pointer for left->right direction
                PathCost* tmp = ptrL1Pre;
                ptrL1Pre = ptrL1Cur;
                ptrL1Cur = tmp;
            }

            //swap buffer pointer for top->bottom direction
            PathCost* tmp = ptrL3PreRow;
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
    
    for(int y = 0; y< height; y++) {
        unsigned* SpPtr = Sp + y*costPerRowEntry;
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
    
    if(subpixelRefine) {
        /*
         * do subpixel quadratic interpolation:
         *fit parabola into (x1=d-1, y1=C[d-1]), (x2=d, y2=C[d]), (x3=d+1, y3=C[d+1])
         *then find minimum of the parabola
         */
        for(int y = 0; y< height; y++) {
            unsigned* SpPtr = Sp + y*costPerRowEntry;
            
            for (int x = 0; x <width; x++) {

				unsigned* ptrSpCur = SpPtr + dMax*x;
                unsigned bestIdx = bestD[y*width + x];
                

				if (bestIdx > 1 && bestIdx < dMax) {
					double c_1 = double(ptrSpCur[bestIdx-1]);
					double c = double(ptrSpCur[bestIdx]);
					double c1 = double(ptrSpCur[bestIdx+1]);
					double bestSubIdx = bestIdx;
					if (c1 < c_1)
						bestSubIdx = bestSubIdx + (c1-c_1)/(c - c_1)/2.0;
					else
						bestSubIdx = bestSubIdx + (c1-c_1)/(c - c1)/2.0;
					
					bestD[y*width + x] = bestSubIdx * (1<<SUBPIXEL_PRECISION);
				}
            }
        }
        
    }
       
       
    mxFree(L1);
    mxFree(L2);
    mxFree(L3);
    mxFree(L4);
    mxFree(Sp);
}
        
void calc_cost(unsigned char* C, 
			 PixelType* I1, PixelType* I2, int width, int height,
			int dMax, double vMax, double* pixelPosD0, double* normlizeDirection, double* offsetFromPosD0)
{
	const int aggWinRadius = 2;
	const int cenWinRadius = 2;
	unsigned* cen1 = (unsigned*)mxMalloc(width * height * sizeof(unsigned));
	unsigned *cen2 = (unsigned*)mxMalloc(width * height * sizeof(unsigned));

	census(I1, cen1, width, height, cenWinRadius);
	census(I2, cen2, width, height, cenWinRadius);

    int winPixels = (2 * aggWinRadius + 1)*(2 * aggWinRadius + 1);

	double* normlizeDirectionX = normlizeDirection;
	double* normlizeDirectionY = normlizeDirection + width*height;

	double* offsetFromPosD0X = offsetFromPosD0;
	double* offsetFromPosD0Y = offsetFromPosD0 + width*height;

	double* pixelPosD0X = pixelPosD0;
	double* pixelPosD0Y = pixelPosD0 + width*height;

    for (int y = 0; y< height; y++) {
        for (int x = 0; x< width; x++) {
            CostType* ptrC = C + y*dMax*width + dMax*x;

			double refPosD0X = pixelPosD0X[y*width + x];
			double refPosD0Y = pixelPosD0Y[y*width + x];
            
            int d = 0;

			/*
			%get 2D offset of disparity from 0 to maximun   
				offset = vzInd.*O(j, i).*repmat(dnorm, 1, dMax+1);
			pt = pd0 + offset;
			pt = round(pt);

			pt(1,:) = min(cols, max(1, pt(1,:)));
			pt(2,:) = min(rows, max(1, pt(2,:)));
			ind2 = sub2ind([rows, cols], pt(2,:), pt(1,:));

			cenCur = cen1Flat(:, ind);
			cenRef = cen2Flat(:, ind2);

			diff = cenCur ~= cenRef;

			cost = sum(diff);
			CFlat(:, ind) = cost'; 
            */


        }
    }


	mxFree(cen1);
	mxFree(cen2);
}
/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    PixelType *I1;             /* pointer to Input image I1 */
    PixelType *I2;             /* pointer to Input image I2 */
    
    mwSize width;              /* rows (width) assume I1/I2 are permuted before passing in */
    mwSize height;             /* cols (height)*/
    
    I1 = (PixelType*)mxGetData(prhs[0]);
    I2 = (PixelType*)mxGetData(prhs[1]);

	int dMax = mxGetScalar(prhs[2]);
	double vMax = mxGetScalar(prhs[3]);
	double* pixelPosD0 = mxGetPr(prhs[4]);
	double* normlizeDirection = mxGetPr(prhs[5]);
	double* offsetFromPosD0 = mxGetPr(prhs[6]);

	int P1 = 6;//mxGetScalar(prhs[8]);
	int P2 = 64;//mxGetScalar(prhs[9]);

    width = mxGetM(prhs[0]);
    height = mxGetN(prhs[0]);
    
    /* create the output matrix */
    const mwSize dims2[] = { width, height };       //output: best index map
    const mwSize dims3[] = { width, height };    //output: cost corresponds to best Index map

    plhs[0] = mxCreateNumericArray(2, dims2, mxUINT32_CLASS, mxREAL);
    plhs[1] = mxCreateNumericArray(2, dims2, mxUINT32_CLASS, mxREAL);

    unsigned * bestD = (unsigned*) mxGetData(plhs[0]);
    unsigned* minC = (unsigned*)mxGetData(plhs[1]);

	//allocate temporal buffers
	CostType* C = mxMalloc(width * height * dMax* sizeof(CostType));
    //construct cost volume
    calc_cost(C, I1, I2, width, height, dMax, vMax, pixelPosD0, normlizeDirection, offsetFromPosD0);

    //perform sgm
    sgm2d(bestD, minC,
        I1, C, width, height, dMax, 
        P1,  P2);

    mxFree(C);
}