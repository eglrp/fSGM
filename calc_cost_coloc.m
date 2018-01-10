function C = calc_cost_coloc(I1, I2, preMv, halfSearchWinSize, aggSize)

if nargin < 3
    halfSearchWinSize = 2;
end

assert(isequal(size(I1), size(I2)));

cenWinSize = 2;
cen1 = census(I1, cenWinSize);
cen2 = census(I2, cenWinSize);

rows = size(I1, 1);
cols = size(I1, 2);
numPixels = rows*cols;
dMax = (4*halfSearchWinSize+1)*(2*halfSearchWinSize+1);
maxHDCost = (2*cenWinSize+1).^2 - 1;
C = uint8(maxHDCost*ones(rows, cols, 2*halfSearchWinSize+1));
cen1Flat = false(size(cen1, 3), numPixels);
cen2Flat = false(size(cen2, 3), numPixels);
CFlat = uint8(maxHDCost*ones(dMax, numPixels));

for i = 1:size(cen1, 3)
    cen1Flat(i, :) = reshape(cen1(:,:,i), 1, []);
    cen2Flat(i, :) = reshape(cen2(:,:,i), 1, []);
end
assert(numel(cen1Flat) == numel(cen1));

tic;

[offx, offy] = meshgrid(-2*halfSearchWinSize:2*halfSearchWinSize, -halfSearchWinSize:halfSearchWinSize);

offx = reshape(offx, 1, []);
offy = reshape(offy, 1, []);

numAggPixels = aggSize * aggSize;
ofx = repmat(offx, numAggPixels, 1);
ofy = repmat(offy, numAggPixels, 1);

halfAggSize = floor(aggSize/2);
[aggOffx, aggOffy] = meshgrid(-halfAggSize:halfAggSize, -halfAggSize:halfAggSize);

%main loop
for ind = 1:numPixels
    
    [j, i] = ind2sub([rows, cols], ind);
    
    curPixelPosY = min(rows, max(1, j + aggOffy));
    curPixelPosX = min(cols, max(1, i + aggOffx));
    
    indCur = sub2ind([rows, cols] , curPixelPosY(:), curPixelPosX(:));
    
    mvxPre = preMv(j, i, 1);
    mvyPre = preMv(j, i, 2);
    
    cX = repmat(curPixelPosX(:), 1, dMax);
    cY = repmat(curPixelPosY(:), 1, dMax);
    
    tx = min(cols, max(1, round(ofx + cX + mvxPre)));
    ty = min(rows, max(1, round(ofy + cY + mvyPre)));
    
    ind2 = sub2ind([rows, cols], ty, tx);
    
    cenRef = cen2Flat(:, ind2);
    cenRef = reshape(cenRef, maxHDCost, [], dMax);
    cenCur = cen1Flat(:, indCur);

    diff = cenCur ~= cenRef;
    cost = sum(sum(diff));
    CFlat(:, ind) = reshape(cost, [], 1)/numAggPixels; 
 
end


for w = 1:dMax
    C(:,:, w) = reshape(CFlat(w, :), rows, cols);
    
end


toc;
