close all;
% clear;

dMax = 64;
vMax = 0.3;
n = dMax+1;
halfWinSize = 2;
postProcessing = 0;

debug = 0;

datasetPath = 'C:/Users/megamusz/Desktop/KITTI/';
outlierAll = 0;
epeAll = 0;
seqNum = 0;

for seqIndex = 0:193
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

[F, E, H, epi, direction, Pd0, normlizeDirection, O, flowR, status] = epipolar_geometry(I0, I1, K);

C = calc_cost(I0, I1, epi, flowR, gtFlow, dMax, halfWinSize, vMax, Pd0, normlizeDirection, O);

[D1, minC, L1, L2, L3, L4] = sgm(C, 6, 64);
% D2 = calc_disp_from_first(D1, Pd0, normlizeDirection,  O, vMax, n );
% D1 = forward_backward_check(D1, D2, Pd0, normlizeDirection, O, vMax, n);
disparites = vzInd2Disp(D1, O, vMax, n);
flowT = disparites.*normlizeDirection;
flow = flowT + flowR;
flow(:,:,3) = ~isnan(D1);


filterD1 = speckle_filter(D1, 2, 100);
filterD2 = calc_disp_from_first(filterD1, Pd0, normlizeDirection,  O, vMax, n );
filterD1 = forward_backward_check(filterD1, filterD2, Pd0, normlizeDirection, O, vMax, n);
filterD1 = speckle_filter(filterD1, dMax, rows*cols/10);
filterD1 = scanline_in_fill(filterD1);
filterdisparites = vzInd2Disp(filterD1, O, vMax, n);
flowT2 = filterdisparites.*normlizeDirection;
flow2 = flowT2 + flowR;
flow2(:,:,3) = ~isnan(filterD1);

if(debug)
    
    D1Image = disp_to_color(D1);
    filterD1(isnan(filterD1)) = 0;
    filterD1Image = disp_to_color(filterD1);
    imwrite([D1Image; filterD1Image], 'sgm_before_vs_after_filter.png');
    
end

if(status)
    flow(:,:,1:2) =0;
end

if(debug)
    thr = [3; 0.05];
    [outlier, aepe] = flow_error(gtFlow, flow, thr)
    epeAll = epeAll + aepe;
    outlierAll = outlierAll + outlier;

    eimage   =  flow_error_image(gtFlow, flow, thr);
    if postProcessing
        flowFinalPP = vmf(flow);
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
    Fc = flow_to_color(flow, 1.2*max(gtFlow(:)));
    gtFlowColor = flow_to_color(gtFlow, 1.2*max(gtFlow(:)));
    imshow([Fc;gtFlowColor]);
    saveas(gcf, sprintf('%06d_color.png', seqIndex));


    flow_field(flow, I0, gtFlow);
    saveas(gcf, sprintf('%06d_flow.png', seqIndex));
end

flow_write_kitti(flow, sprintf('%06d_10.png', seqIndex));
flow_write_kitti(flow2, sprintf('%06d_10_pp.png', seqIndex));

seqNum = seqNum + 1;
end


outlierAll = outlierAll/seqNum;
epeAll = epeAll/seqNum;


