function [ flow , minC] = ng_sgm( I0, I1 )
% Calculate flow from I0 to I1 using neighbour guided SGM Optical flow Method
% 
%   Detailed explanation goes here
% main parameters
    
    P1 = 6;
    P2 = 32;
    
    % loop pyramidal levels
    row = size(I0, 1);
    col = size(I1, 2);

    tic;
    I1gray = rgb2gray(permute(I0, [2, 1, 3]));
    I2gray = rgb2gray(permute(I1, [2, 1, 3]));
    temporalHints = permute(zeros(row, col), [2, 1, 3]);
    %construct cost volume and SGM

    [minC, flow] = calc_cost_sgm_ng(I1gray, I2gray, temporalHints, 1, 2, 0, P1, P2);
    flow = permute(flow, [2, 1, 3]);
    toc;
        
end

