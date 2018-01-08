function [ mvCurLevel, mvPyd , minC] = pyramidal_sgm( I0, I1, numPyd )
%PYRAMIDAL_SGM Summary of this function goes here
%   Detailed explanation goes here

    I0pyd{1} = I0;
    I1pyd{1} = I1;
    mvPyd = {numPyd};

    
    
    %create image pyramid
    for l = 2:numPyd
        I0pyd{l} = impyramid(I0pyd{l-1}, 'reduce');
        I1pyd{l} = impyramid(I1pyd{l-1}, 'reduce');
    end
    
    mvPreLevel = zeros(size(I0pyd{numPyd}, 1), size(I0pyd{numPyd}, 2), 2);
    P1 = 6;
    P2 = 64;
    % loop pyramidal levels
    for l = numPyd:-1:1
        verSearchHalfWinSize = 4;
        rowl = size(I0pyd{l}, 1);
        coll = size(I0pyd{l}, 2);
        mvCurLevel = zeros(rowl, coll, 2);
        
        %construct cost volume
        aggSize = 5;

        
        tic;
        I1gray = rgb2gray(permute(I0pyd{l}, [2, 1, 3]));
        I2gray = rgb2gray(permute(I1pyd{l}, [2, 1, 3]));
        mvPrePermuted = permute(mvPreLevel, [2, 1, 3]);
        C = calc_cost_pyd(I1gray, I2gray, mvPrePermuted, verSearchHalfWinSize, aggSize);
        C = permute(C, [2, 1, 3]);
        toc;
%         C = calc_cost_coloc(I0pyd{l}, I1pyd{l}, mvPreLevel, verSearchHalfWinSize, aggSize);

        %2d sgm
  

        [minIdx, minC, mvSub] = sgm2d(C, 2*verSearchHalfWinSize+1, 4*verSearchHalfWinSize + 1, P1, P2, 1);

        %WTA
%         [minC, minIdx] = min(C, [], 3);
%         mvSub = zeros(rowl, coll, 2);
       
        
        % recover mv from idx
        [r, c] = ind2sub([2*verSearchHalfWinSize+1, 4*verSearchHalfWinSize + 1], minIdx(:));
        mvx = c - 2*verSearchHalfWinSize - 1;
        mvy = r - verSearchHalfWinSize - 1;
       
        mvCurLevel(:,:,1) = reshape(mvx, [rowl, coll]);
%         mvCurLevel(:,:,1) = medfilt2(mvCurLevel(:,:,1), [3 3]);
        mvCurLevel(:,:,2) = reshape(mvy, [rowl, coll]);
%         mvCurLevel(:,:,2) = medfilt2(mvCurLevel(:,:,2), [3 3]);
        mvCurLevel = mvCurLevel + mvPreLevel(1:rowl, 1:coll, :) + mvSub;
       
        
        mvPyd{l} = mvCurLevel;
        if (l > 1)
            %pass to next level, need to upscale mv map size and also the
            %mv magnitude. 
            mvPreLevel = 2*imresize(mvCurLevel, 2, 'nearest');
        end

    end
end

