close all;
% clear;

dMax = 64;
vMax = 0.3;
n = dMax+1;
halfWinSize = 2;
postProcessing = 0;

debug = 1;

datasetPath = 'C:/Users/megamusz/Desktop/KITTI/';
outlierAll = 0;
epeAll = 0;
seqNum = 0;

for seqIndex = 14:14
I0 = imread(sprintf('%s/training/image_2/%06d_10.png', datasetPath, seqIndex));
I1 = imread(sprintf('%s/training/image_2/%06d_11.png', datasetPath, seqIndex));
gtFlow = flow_read_kitti(sprintf('%s/flow_occ/%06d_10.png',datasetPath, seqIndex));
P0 = read_calib_file(sprintf('%s/training/calib_cam_to_cam/%06d.txt', datasetPath, seqIndex), 1);

% I0 = imread(sprintf('%s/training/image_0/%06d_10.png', datasetPath, seqIndex));
% I1 = imread(sprintf('%s/training/image_0/%06d_11.png', datasetPath, seqIndex));
% gtFlow = flow_read_kitti(sprintf('%s/training/flow_occ/%06d_10.png', datasetPath, seqIndex));
% P0 = read_calib_file(sprintf('%s/training/calib/%06d.txt', datasetPath, seqIndex));
K = P0(1:3, 1:3);

rows = size(I0, 1);
cols = size(I0, 2);


numPyd = 3;
[ mvCurLevel, mvPyd, minCpyd  ] = pyramidal_sgm( I0, I1, numPyd );

mvCurLevel(:,:,3) = 1;


if(debug)
    figure;
    for l = numPyd:-1:1
        subplot(2, 3, l);
        imshow(flow_to_color(mvPyd{numPyd - l + 1}));
        title(['Level: ', num2str(l)]);
    end
    saveas(gcf, sprintf('pyd_mv_%06d.png', seqIndex));
    
    thr = [3; 0.05];
    [outlier, aepe] = flow_error(gtFlow, mvCurLevel, thr);
    eimage   =  flow_error_image(gtFlow, mvCurLevel, thr);

    figure;
    imshow(eimage);
    title(['outlier percentage: ', num2str(outlier), ' AEPE: ', num2str(aepe)]);
    saveas(gcf, sprintf('%06d_error.png', seqIndex));


    figure;
    Fc = flow_to_color(mvCurLevel, 1.2*max(gtFlow(:)));
    gtFlowColor = flow_to_color(gtFlow, 1.2*max(gtFlow(:)));
    imshow([Fc;gtFlowColor]);
    saveas(gcf, sprintf('%06d_color.png', seqIndex));


    flow_field(mvCurLevel, I0, gtFlow);
    saveas(gcf, sprintf('%06d_flow.png', seqIndex));
end

flow_write_kitti(mvCurLevel, sprintf('%06d_10.png', seqIndex));


seqNum = seqNum + 1;
end
