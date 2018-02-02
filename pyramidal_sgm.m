function [ mvCurLevel, mvPyd , minC] = pyramidal_sgm( I0, I1, numPyd )
% Calculate flow from I0 to I1 using Pyramidal SGM Optical flow Method
% 
%   Detailed explanation goes here
% main parameters
    if(nargin < 3)
        %default pyramidal levels
        % the maximun search range is 2^(numPyd-1) *
        % [-2*verSearchHalfWinSize 2*verSearchHalfWinSize] in horizontal
        % direction; [-verSearchHalfWinSize verSearchHalfWinSize] in
        % vertical direction
        numPyd = 5;
    end
    
    P1 = 6;
    P2 = 32;
    aggHalfWinSize = 2;                %half aggregation window size
    verSearchHalfWinSize = 5;   %half search window size in vertical direction
    horSearchHalfWinSize = 5;   %half search window size in vertical direction
    enableDiagonal = 1;
    totalPass = 2;
    adaptiveP2 = 0;
    
    I0pyd{1} = I0;
    I1pyd{1} = I1;
    mvPyd = {numPyd};
    %create image pyramid
    for l = 2:numPyd
        I0pyd{l} = impyramid(I0pyd{l-1}, 'reduce');
        I1pyd{l} = impyramid(I1pyd{l-1}, 'reduce');
    end
  
    %initialize previous level's mv result as zero
    mvPreLevel = zeros(size(I0pyd{numPyd}, 1), size(I0pyd{numPyd}, 2), 2);
    
    % loop pyramidal levels
    for l = numPyd:-1:1
        rowl = size(I0pyd{l}, 1);
        coll = size(I0pyd{l}, 2);
        mvCurLevel = zeros(rowl, coll, 2);
        
        tic;
        % permute the matrix for row-major processing order in c code. 
        I1gray = rgb2gray(permute(I0pyd{l}, [2, 1, 3]));
        I2gray = rgb2gray(permute(I1pyd{l}, [2, 1, 3]));
        mvPrePermuted = permute(mvPreLevel, [2, 1, 3]);
        
        %construct cost volume and SGM
        subpixelRefine = l == 1;
        [minIdx, minC, mvSub] = calc_pyd_cost_sgm(I1gray, I2gray, mvPrePermuted, horSearchHalfWinSize, verSearchHalfWinSize, aggHalfWinSize, subpixelRefine, ...
                P1, P2, enableDiagonal, totalPass, adaptiveP2);
        minIdx = minIdx' + 1; 
        mvSub = permute(mvSub, [2, 1, 3]);
        toc;
        
        % recover mv from idx
        [r, c] = ind2sub([2*verSearchHalfWinSize+1, 2*horSearchHalfWinSize + 1], minIdx(:));

        mvx = c - horSearchHalfWinSize - 1;
        mvy = r - verSearchHalfWinSize - 1;

        mvCurLevel(:,:,1) = reshape(mvx, [rowl, coll]);       
        mvCurLevel(:,:,2) = reshape(mvy, [rowl, coll]);
        mvCurLevel = mvCurLevel + mvPreLevel(1:rowl, 1:coll, :) + mvSub;

        mvPyd{l} = mvCurLevel;
        if (l > 1)
            %pass to next level, need to upscale mv map size and also the
            %mv magnitude.               
%             mvCurLevel(:,:,1) = medfilt2(mvCurLevel(:,:,1));
%             mvCurLevel(:,:,2) = medfilt2(mvCurLevel(:,:,2));
            mvPreLevel = 2*imresize(mvCurLevel, 2, 'nearest');
        end

    end
end

