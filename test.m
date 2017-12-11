close all;
clear;

seqIndex = 0;

I0 = imread(sprintf('%06d_10.png', seqIndex));
I1 = imread(sprintf('%06d_11.png', seqIndex));
gtFlow = flow_read_kitti('C:\Users\megamusz\Desktop\KITTI\flow_occ\000000_10.png');
% F1 = estimate_Fmatrix(I0, I1);

P0 = read_calib_file(sprintf('%06d.txt', seqIndex), 1);

% 
% [w, t, E2, e1, e2, F2] = mono_vo(I0, I1, [P0(1,1) P0(2,2)], [P0(1, 3) P0(2, 3)]);

load workspace3.mat

flow = rotation_motion(w, P0(1,1), [P0(1, 3) P0(2, 3)], [size(I0, 2), size(I0, 1)]);

% 
dMax = 128;
vMax = 0.4;
halfWinSize = 2;

[C, normlizeDirection, O] = calc_cost(I0, I1, e2, flow, gtFlow, dMax, halfWinSize, vMax);
% 
% save('workspace.mat', 'C', 'e1', 'e2', 'flow', 'w', 't', 'E2', 'F2', 'normlizeDirection', 'O');
% load workspace.mat;

[bestD, L1, L2, L3, L4] = sgm(C);

n = dMax+1;
bestD = vzInd2Disp(bestD-1, O, vMax, n);
flowT = (bestD-1).*normlizeDirection;

flowFinal = flowT + flow;
flowFinal(:,:,3) = 1;

% gtFlow = flow_read_kitti('000004_10_gtflow.png');


eimage  =  flow_error_image(gtFlow, flowFinal, [3; 0.05]);
figure;
imshow(eimage);
[outlier, aepe] = flow_error(gtFlow, flowFinal, [3; 0.05]);
title(['outlier percentage: ', num2str(outlier), ' AEPE: ', num2str(aepe)]);
figure;
Fc = flow_to_color(flowFinal);
imshow(Fc);


function D = vzInd2Disp(w, O, vMax, n)
    vzRatio = w./n*vMax;
    vzInd = vzRatio ./ (1-vzRatio);
    D = O.*vzInd;
end
