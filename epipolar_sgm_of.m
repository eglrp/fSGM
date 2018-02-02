function [ flow, flowPP, conf, minC, status ] = epipolar_sgm_of( I0, I1, K, dMax, vMax, post_processing)
%EPIPOLAR_SGM_OF Summary of this function goes here
%   Detailed explanation goes here
debug = 1;
if(nargin < 4)
    dMax = 64;
end
if(nargin < 5)
    vMax = 0.3;
end
if(nargin < 6)
    post_processing = 0;
end
P1 = 6; P2 = 64;
n = dMax+1;

rows = size(I0, 1);
cols = size(I0, 2);


[F, E, H, epi, direction, Pd0, normlizeDirection, O, flowR, status] = epipolar_geometry(I0, I1, K);

if(status)
    flow = zeros(rows, cols, 3);
    flow(:,:,3) = 1;
    flowPP = flow;
    conf = zeros(rows, cols);
    minC = zeros(rows, cols);
    return;
end


I0_ = permute(I0, [2, 1, 3]);
I1_ = permute(I1, [2, 1, 3]);
if(size(I0_, 3) > 1)
    I0_ = rgb2gray(I0_);
    I1_ = rgb2gray(I1_);
end
Pd0_ = permute(Pd0, [2, 1, 3]);
normlizeDirection_ = permute(normlizeDirection, [2, 1, 3]);
O_ = O';
tic;
[bestD, minC, conf, D2] = calc_cost_sgm(I0_, I1_, dMax, vMax, Pd0_, normlizeDirection_, O_, P1, P2);
disparites = double(bestD')/256.0;
toc;
% C = calc_cost(I0, I1, epi, dMax, halfWinSize, vMax, Pd0, normlizeDirection, O);

% [D1, minC, L1, L2, L3, L4] = sgm(C, 6, 64);
% D2 = calc_disp_from_first(D1, Pd0, normlizeDirection,  O, vMax, n );
% D1 = forward_backward_check(D1, D2, Pd0, normlizeDirection, O, vMax, n);
% disparites = vzInd2Disp(D1, O, vMax, n);
flowT = disparites.*normlizeDirection;
flow = flowT + flowR;
% flow(:,:,3) = ~isnan(D1);
flow(:,:,3) = 1;

if(post_processing)
    filterD1 = speckle_filter(D1, 2, 100);
    filterD2 = calc_disp_from_first(filterD1, Pd0, normlizeDirection,  O, vMax, n );
    filterD1 = forward_backward_check(filterD1, filterD2, Pd0, normlizeDirection, O, vMax, n);
    conf = ~isnan(filterD1);

    filterD1 = speckle_filter(filterD1, dMax, rows*cols/10);
    filterD1 = scanline_in_fill(filterD1);
    filterdisparites = vzInd2Disp(filterD1, O, vMax, n);
    flowT2 = filterdisparites.*normlizeDirection;

    flowPP = flowT2 + flowR;
    flowPP(:,:,3) = ~isnan(filterD1);
else
    flowPP = flow;
    conf = true(rows, cols);
end

if(status)
    flow(:,:,1:2) =0;
end

end

