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
        unsigned char* C, int width, int height, int dMax, 
        double* mvPre, int mvWidth, int mvHeight, 
        int searchWinX, int searchWinY, int P1, int P2, int subpixelRefine )
{
    
    //bestD = (unsigned*) malloc(sizeof(unsigned)  * width * height);
    minC = (unsigned*) malloc(sizeof(unsigned)  * width * height);
    mvSub = (double*) malloc(sizeof(double)  * width * height * 2);
    
    unsigned char* L1 = (unsigned char*) malloc (sizeof(unsigned char) * width * height * dMax); //
    unsigned char* L2 = (unsigned char*) malloc (sizeof(unsigned char) * width * height * dMax); //
    unsigned char* L3 = (unsigned char*) malloc (sizeof(unsigned char) * width * height * dMax); //
    unsigned char* L4 = (unsigned char*) malloc (sizeof(unsigned char) * width * height * dMax); //
    unsigned char* Sp = (unsigned char*) malloc (sizeof(unsigned char) * width * height * dMax); //
    unsigned char minL[4];
    
    unsigned char* minL3 = (unsigned char*) malloc(sizeof(unsigned char)* width);
    
    double* pMvx = mvPre;
    double* pMvy = mvPre + mvWidth * mvHeight;
    
    int pass = 0;
    int ystart = 0;
    int yend = height-1;
    int ystep = 1;

    int xstart = 0;
    int xend = width-1;
    int xstep = 1;
    
    if(pass == 1) {
        ystart = height-1;
        yend = 0;
        ystep = -1;

        xstart = width -1;
        xend = 0;
        xstep = -1;
    } 
    
    int pathCostPerRowEntry = width * dMax;
    
    int* LUTY = (int*) malloc(dMax * sizeof(int));
    int* LUTX = (int*) malloc(dMax * sizeof(int));
    
    for (int d=0; d < dMax; d++) {
        int y = d % searchWinY;
        int x = d / searchWinY;
        LUTX[d] = x;   
        LUTY[d] = y;
    }
    
    
    for(int y = ystart; y <= yend; y += ystep) {
        unsigned char* L1ptr = L1 + y*pathCostPerRowEntry;
        unsigned char* L2ptr = L2 + y*pathCostPerRowEntry;
        unsigned char* L3ptr = L3 + y*pathCostPerRowEntry;
        unsigned char* L4ptr = L4 + y*pathCostPerRowEntry;
        unsigned char* Cptr = C + y*pathCostPerRowEntry;
            
        for (int x= xstart; x <= xend ; x += xstep) {
            
            if(x == xstart) {
                minL[0] = minL[1] = 255;
                for(int d = 0; d < dMax; d++) {
                    
                    L1ptr[xstart*dMax + d] = Cptr[xstart*dMax + d];
                    L2ptr[xstart*dMax + d] = Cptr[xstart*dMax + d];
                    
                    if(L1ptr[d] < minL[0]) {
                        minL[0] = L1ptr[d];
                    }
                    
                    if(L2ptr[d] < minL[1]) {
                        minL[1] = L2ptr[d];
                    }
                }
            } 
            
            if(y == ystart) {
                minL[2] = minL3[x] = 255;
                
                for(int d = 0; d < dMax; d++) {
                    
                    L3ptr[ystart* pathCostPerRowEntry + x*dMax + d] = Cptr[ystart* pathCostPerRowEntry + x*dMax + d];
                    L4ptr[ystart* pathCostPerRowEntry + x*dMax + d] = Cptr[ystart* pathCostPerRowEntry + x*dMax + d];
                    
                    if(L3ptr[d] < minL[2]) {
                        minL3[x] = L3ptr[d];
                    }
                    
                    if(L4ptr[d] < minL[3]) {
                        minL[3] = L4ptr[d];
                    }
                }
            }
            
            
            if(x > xstart) {
                double dx = pMvx[y*width + x] - pMvx[y*width + x - xstep];
                double dy = pMvy[y*width + x] - pMvy[y*width + x - xstep];
                
                int dOff = dx * searchWinY + dy;
                int minL1 = 255;
                for (int d = 0; d < dMax; d++) {
                    int sy = LUTY[d];
                    int sx = LUTX[d];
                    unsigned bestC = 255;
                    
                    for (int d2 = 0; d2 < dMax; d2++) {
                        int dpre = d2 - dOff;
                        int sy2, sx2;
                        unsigned L1pre_d = 255;
                        if(dpre >0 && dpre <dMax) {
                            L1pre_d = L1ptr[(x- xstep)*dMax + dpre];
                            sy2 = LUTY[dpre];
                            sx2 = LUTX[dpre];
                        }else {
                            sy2 = sy + 3;
                            sx2 = sx + 3;
                        }
                        
                        int penalty = P2;
                        if(sx2 == sx && sy2 == sy){
                            penalty = 0;
                        } else if(abs(sx2-sx) <= 2&& abs(sy2-sy) <=2) {
                            penalty = P1;
                        }
                        
                        bestC = min(bestC, L1pre_d + penalty);
                        
                    }
                    L1ptr[x*dMax + d] = Cptr[x*dMax + d] + bestC - minL[0];
                    minL1 = min(L1ptr[x*dMax + d], minL1);
                }
                
                minL[0] = minL1;
            }
            
            if(y > ystart) {
                double dx = pMvx[y*width + x] - pMvx[(y-ystep)*width + x];
                double dy = pMvy[y*width + x] - pMvy[(y-ystep)*width + x];
                
                int dOff = dx * searchWinY + dy;
                int minL = 255;
                for (int d = 0; d < dMax; d++) {
                    int sy = LUTY[d];
                    int sx = LUTX[d];
                    unsigned bestC = 255;
                    
                    for (int d2 = 0; d2 < dMax; d2++) {
                        int dpre = d2 - dOff;
                        int sy2, sx2;
                        unsigned L3pre_d = 255;
                        
                        
                        if(dpre >0 && dpre <dMax) {
                            L3pre_d = L3ptr[-ystep*pathCostPerRowEntry + x*dMax + dpre];
                            sy2 = LUTY[dpre];
                            sx2 = LUTX[dpre];
                        }else {
                            sy2 = sy + 3;
                            sx2 = sx + 3;
                        }
                        
                        int penalty = P2;
                        if(sx2 == sx && sy2 == sy){
                            penalty = 0;
                        } else if(abs(sx2-sx) <= 2&& abs(sy2-sy) <=2) {
                            penalty = P1;
                        }
                        
                        bestC = min(bestC, L3pre_d + penalty);
                        
                    }
                    
                    L3ptr[x*dMax + d] = Cptr[x*dMax + d] + bestC - minL3[x];
                    minL = min(L3ptr[x*dMax + d], minL);
                }
                
                minL3[x] = minL;
            }
            
            if(x > xstart && y > ystart) {
            }
            
        }
    }
    
    
    for(int y = 0; y< height; y++) {
        unsigned char* L1ptr = L1 + y*pathCostPerRowEntry;
        unsigned char* L3ptr = L3 + y*pathCostPerRowEntry;
        for (int x = 0; x <width; x++) {
            
            unsigned minCost = 255*4;
            unsigned minIdx = 0;
            for (int d = 0; d<dMax; d++) {
                unsigned costSum = L1ptr[x*dMax+d] + L3ptr[x*dMax+d];
                if(costSum < minCost) {
                    minCost = costSum;
                    minIdx = d;
                }
            }
            
            minC[y*width +x ] = minCost;
            bestD[y*width +x] = minIdx;
            
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
    const mwSize dims2[]={width, height};
    plhs[1] = mxCreateNumericArray(2, dims2, mxUINT32_CLASS, mxREAL);
    
    unsigned char* C = (unsigned char*) mxGetData(plhs[0]);
    unsigned * bestD = (unsigned*) mxGetData(plhs[1]);
        
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
    

    unsigned* minC;
    double* mvSub;
    int P1 = 6;
    int P2 = 64;
    
    unsigned char* Cre = (unsigned char*) malloc(width*height*dMax*sizeof(unsigned char));
    for(int d = 0; d< dMax; d++) { 
        for (int y = 0; y< height; y++) {
            for(int x = 0; x <width; x++) {
                Cre[y* width *dMax + x*dMax + d] = C[d*width*height + (y*width + x)];
            }
        }
    }
    sgm2d(bestD, minC, mvSub, 
        Cre, width, height, dMax, 
        preMv, mvWidth, mvHeight, 
        winRadiusX*2 +1, winRadiusY*2+1, P1,  P2, 0);
}