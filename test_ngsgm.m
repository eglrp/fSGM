close all;
% clear;

dMax = 64;
vMax = 0.3;
n = dMax+1;
halfWinSize = 2;
postProcessing = 0;

debug = 1;

datasetPath = 'C:/Users/Dong/Documents/dataset/data_scene_flow/';
outlierAll = 0;
epeAll = 0;
seqNum = 0;

for seqIndex = 14:14
I0 = imread(sprintf('%s/training/image_2/%06d_10.png', datasetPath, seqIndex));
I1 = imread(sprintf('%s/training/image_2/%06d_11.png', datasetPath, seqIndex));
gtFlow = flow_read_kitti(sprintf('%s/training/flow_occ/%06d_10.png',datasetPath, seqIndex));
P0 = read_calib_file(sprintf('%s/training/calib_cam_to_cam/%06d.txt', datasetPath, seqIndex), 1);

% I0 = imread(sprintf('%s/training/image_0/%06d_10.png', datasetPath, seqIndex));
% I1 = imread(sprintf('%s/training/image_0/%06d_11.png', datasetPath, seqIndex));
% gtFlow = flow_read_kitti(sprintf('%s/training/flow_occ/%06d_10.png', datasetPath, seqIndex));
% P0 = read_calib_file(sprintf('%s/training/calib/%06d.txt', datasetPath, seqIndex));
K = P0(1:3, 1:3);

rows = size(I0, 1);
cols = size(I0, 2);

[ flow , minC ] = ng_sgm( I0, I1 );
flow(:,:,3) = 1;


thr = [3; 0.05];
[outlier, aepe] = flow_error(gtFlow, flow, thr);
epeAll = epeAll + aepe;
outlierAll = outlierAll + outlier;
eimage   =  flow_error_image(gtFlow, flow, thr);


figure;
imshow(eimage);
title(['outlier percentage: ', num2str(outlier), ' AEPE: ', num2str(aepe)]);
saveas(gcf, sprintf('%06d_error.png', seqIndex));


figure;
maxFlow = 1.2*max(max([abs(gtFlow(:, :, 1)); abs(gtFlow(:, :, 2))]));
Fc = flow_to_color(flow, maxFlow );
gtFlowColor = flow_to_color(gtFlow, maxFlow);
imshow([Fc;gtFlowColor]);
saveas(gcf, sprintf('%06d_color.png', seqIndex));


flow_field(flow, I0, gtFlow);
saveas(gcf, sprintf('%06d_flow.png', seqIndex));


flow_write_kitti(flow, sprintf('%06d_10.png', seqIndex));

end
