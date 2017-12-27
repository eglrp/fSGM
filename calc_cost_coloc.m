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
C = 65535 * ones(rows, cols, 2*halfSearchWinSize+1);
cen1Flat = zeros(size(cen1, 3), numPixels);
cen2Flat = zeros(size(cen2, 3), numPixels);
CFlat = 65535*ones(dMax, numPixels);

for i = 1:size(cen1, 3)
    cen1Flat(i, :) = reshape(cen1(:,:,i), 1, []);
    cen2Flat(i, :) = reshape(cen2(:,:,i), 1, []);
end
assert(numel(cen1Flat) == numel(cen1));

tic;

[offx, offy] = meshgrid(-2*halfSearchWinSize:2*halfSearchWinSize, -halfSearchWinSize:halfSearchWinSize);

for ind = 1:numPixels
    
    [j, i] = ind2sub([rows, cols], ind);
    
    %get 2D offset of disparity from 0 to maximun   
    mvxPre = preMv(j, i, 1);
    mvyPre = preMv(j, i, 2);
    targetX = min(cols, max(1, round(offx + i + mvxPre)));
    targetY = min(rows, max(1, round(offy + j + mvyPre)));
    
    
    ind2 = sub2ind([rows, cols], reshape(targetY, [], 1), reshape(targetX, [], 1));

    cenCur = cen1Flat(:, ind);
    cenRef = cen2Flat(:, ind2);

    diff = cenCur ~= cenRef;

    cost = sum(diff);
    CFlat(:, ind) = cost'; 
 
end


for w = 1:dMax
    C(:,:, w) = reshape(CFlat(w, :), rows, cols);
    
end

C = imboxfilt3(C, [aggSize aggSize 1], 'Padding', 'replicate');

toc;
