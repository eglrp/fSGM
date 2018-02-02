function [ flow, minC, status ] = epipolar_sgm_of( I0, I1, K, dMax, vMax)
%EPIPOLAR_SGM_OF calll epipolar SGM optical flow to do optical flow
%   Detailed explanation goes here
% Inputs: First/Second image I1/I2, camera intrinsic matrix K
% maximun disparity range dMax, maximun v-z ratio vMax
%
% Outputs: 
% flow: estimated optical flow 
% minC: minimum aggregated path cost from SGM 
% status: 0/1 sucess/fail

if(nargin < 4)
    dMax = 64;
end
if(nargin < 5)
    vMax = 0.3;
end

P1 = 6; P2 = 64;

rows = size(I0, 1);
cols = size(I0, 2);

[Pd0, normlizeDirection, O, flowR, status] = epipolar_geometry(I0, I1, K);

if(status)
    flow = zeros(rows, cols, 3);
    flow(:,:,3) = 1;
    minC = zeros(rows, cols);
    return;
end

I0_ = permute(I0, [2, 1, 3]);
I1_ = permute(I1, [2, 1, 3]);
if(size(I0_, 3) > 1)
    I0_ = rgb2gray(I0_);
    I1_ = rgb2gray(I1_);
end

tic;
Pd0_ = permute(Pd0, [2, 1, 3]);
normlizeDirection_ = permute(normlizeDirection, [2, 1, 3]);
O_ = O';

[bestD, minC] = calc_cost_sgm(I0_, I1_, dMax, vMax, Pd0_, normlizeDirection_, O_, P1, P2);
disparites = double(bestD')/256.0;
toc;

flowT = disparites.*normlizeDirection;
flow = flowT + flowR;
flow(:,:,3) = 1;

% in case that epipolar sgm flow fail, output zero flow map
if(status)
    flow(:,:,1:2) =0;
end

end

