close all;

seqIndex = 4;

I0 = imread(sprintf('%06d_10.png', seqIndex));
I1 = imread(sprintf('%06d_11.png', seqIndex));

% F = estimate_Fmatrix(I0, I1);

% P0 = read_calib_file(sprintf('%06d.txt', seqIndex));
% 
% 
% [w, t, E1, e1, e2, F] = mono_vo(I0, I1, [P0(1,1) P0(2,2)], [P0(1, 3) P0(2, 3)]);
% 
% flow = rotation_motion(w, P0(1,1), [P0(1, 3) P0(2, 3)], [size(I0, 2), size(I0, 1)]);
% 
% % Fc = flow_to_color(flow);
% % imshow(Fc);
% 
% % flow_field(flow, I0);
% 
% 
% C = calc_cost(I0, I1, e2, flow);
% 
% save('workspace.mat', 'C', 'e1', 'e2', 'flow', 'w', 't', 'E1', 'F');

load workspace.mat;
bestD = sgm(C);

rows = size(I0, 1);
cols = size(I0, 2);
P = zeros(rows, cols, 2);
[P(:,:,1), P(:,:,2)] = meshgrid(1:cols, 1:rows);
E2I = zeros(rows, cols, 2);
e2i = [e2(1)/e2(3); e2(2)/e2(3)];
E2I(:,:,1) = repmat(e2i(1), rows, cols);
E2I(:,:,2) = repmat(e2i(2), rows, cols);
normlizeDirection = normlize(P + flow - E2I);
flowT = (bestD-1).*normlizeDirection;

flowFinal = flowT + flow;
flowFinal(:,:,3) = 1;

gtFlow = flow_read_kitti('C:\Users\megamusz\Desktop\KITTI\training\flow_occ\000004_10.png');


eimage  =  flow_error_image(gtFlow, flowFinal, [3; 0.05]);
imshow(eimage);
    
function bestD = sgm(C)

    P1 = 100;
    P2 = 2000;
    
    rows = size(C, 1);
    cols = size(C, 2);
    dMax = size(C, 3);
    
    L1 = zeros(size(C)); %path cost for left -> right direction
    
    for j = 1:rows
        for i = 1:cols
            %initialize the Path cost L at the left bounary
            if i == 1
                L1(j, 1, :) = C(j, 1, :);
                
            else
                for d = 1:dMax
                    minLeft  = min(L1(j, i-1, 1:max(1, d-2)), [], 3);
                    minRight = min(L1(j, i-1, min(dMax, d+2):end), [], 3);
                    
                    L1(j, i, d) = C(j, i, d) + min([L1(j, i-1, d), ...
                                                   L1(j, i-1, min(dMax, d+1)) + P1, ...
                                                   L1(j, i-1, max(1, d-1)) + P1, ...
                                                   min(minLeft, minRight) + P2]);
                end
                
                
            end
            
        end
    end
    
    [~, bestD] = min(L1, [], 3);
    
end



function dnorm = normlize(d)
    if(size(d, 3) > 1)
        dnorm = d./sqrt(sum(d.^2, 3));
    else
        dnorm = d/sqrt(sum(d.^2));
    end
end
