close all;
clear;

seqIndex = 0;

dMax = 64;
vMax = 0.3;
halfWinSize = 2;
postProcessing = 0;

datasetPath = 'C:/Users/megamusz/Desktop/KITTI/';
% I0 = imread(sprintf('%s/training/image_2/%06d_10.png', datasetPath, seqIndex));
% I1 = imread(sprintf('%s/training/image_2/%06d_11.png', datasetPath, seqIndex));
% gtFlow = flow_read_kitti(sprintf('%s/flow_occ/%06d_10.png',datasetPath, seqIndex));
% P0 = read_calib_file(sprintf('%s/training/calib_cam_to_cam/%06d.txt', datasetPath, seqIndex), 1);

I0 = imread(sprintf('%s/training/image_0/%06d_10.png', datasetPath, seqIndex));
I1 = imread(sprintf('%s/training/image_0/%06d_11.png', datasetPath, seqIndex));
gtFlow = flow_read_kitti(sprintf('%s/training/flow_occ/%06d_10.png', datasetPath, seqIndex));
P0 = read_calib_file(sprintf('%s/training/calib/%06d.txt', datasetPath, seqIndex));
K = P0(1:3, 1:3);

rows = size(I0, 1);
cols = size(I0, 2);

[F, E, H, epi] = epipolar_geometry(I0, I1, K);

flow = rotation_motion(H, F, rows, cols);

[C, normlizeDirection, O] = calc_cost(I0, I1, epi, flow, gtFlow, dMax, halfWinSize, vMax);

[bestD, minC, L1, L2, L3, L4] = sgm(C, 6, 64);

n = dMax+1;
bestD = vzInd2Disp(bestD-1, O, vMax, n);
flowT = bestD.*normlizeDirection;

flowFinal = flowT + flow;
flowFinal(:,:,3) = 1;

thr = [3; 0.05];

eimage   =  flow_error_image(gtFlow, flowFinal, thr);
[outlier, aepe] = flow_error(gtFlow, flowFinal, thr);
if postProcessing
    flowFinalPP = vmf(flowFinal);
    eimage2  =  flow_error_image(gtFlow, flowFinalPP, thr);
    figure;
    imshow([eimage;eimage2]);
    [outlier2, aepe2] = flow_error(gtFlow, flowFinalPP, thr);
    title(['outlier percentage: ', num2str(outlier), ' AEPE: ', num2str(aepe), ' After PP: ', 'outlier percentage: ', num2str(outlier2), ' AEPE: ', num2str(aepe2)]);
else 
    figure;
    imshow([eimage]);
    title(['outlier percentage: ', num2str(outlier), ' AEPE: ', num2str(aepe)]);
end

saveas(gcf, sprintf('%06d_error.png', seqIndex));


figure;
Fc = flow_to_color(flowFinal, 1.2*max(gtFlow(:)));
gtFlowColor = flow_to_color(gtFlow, 1.2*max(gtFlow(:)));
imshow([Fc;gtFlowColor]);
saveas(gcf, sprintf('%06d_color.png', seqIndex));


flow_field(flowFinal, I0, gtFlow);
saveas(gcf, sprintf('%06d_flow.png', seqIndex));

function D = vzInd2Disp(w, O, vMax, n)
    vzRatio = w./n*vMax;
    vzInd = vzRatio ./ (1-vzRatio);
    D = O.*vzInd;
end

