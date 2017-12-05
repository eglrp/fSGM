function C = calc_cost(I1, I2, epipole, Rflow, dMax, halfWinSize)

debug = 1;

if debug
    close all;
end

if nargin < 5
    dMax = 16;
end

if nargin < 6
    halfWinSize = 3;
end

% 
% I1 = double(imread('000004_10.png'));
% I2 = double(imread('000004_11.png'));


assert(isequal(size(I1), size(I2)));

% e1 = [-0.959870232857247;-0.280438303284809;-0.001869257957789];
% e2 = [0.965222503953150;0.261423812450535;0.001762993545078];
% 
% F = [1.160024130294298e-08,1.539326689928344e-04,-0.023099942444453;-1.531504691961629e-04,3.313869025373902e-07,0.078593562935847;0.022703419892627,-0.084325838315150,0.992831445999651];
% P0 = [718.856000000000,0,607.192800000000,0;0,718.856000000000,185.215700000000,0;0,0,1,0];
% w = [-0.002409420424081,-0.046632321275709,-0.002414820939297];
% Rflow = rotation_motion(w, P0(1,1), [P0(1, 3) P0(2, 3)], [size(I1, 2), size(I1, 1)]);

rows = size(I1, 1);
cols = size(I1, 2);

e2i = [epipole(1)/epipole(3); epipole(2)/epipole(3)];
% e2i = [e2(1)/e2(3); e2(2)/e2(3)];

C = 65535 * ones(rows, cols, dMax + 1);

P = zeros(rows, cols, 2);
[P(:,:,1), P(:,:,2)] = meshgrid(1:cols, 1:rows);
E2I = zeros(rows, cols, 2);

E2I(:,:,1) = repmat(e2i(1), rows, cols);
E2I(:,:,2) = repmat(e2i(2), rows, cols);
normlizeDirection = normlize(P + Rflow - E2I);

% flow_field(normlizeDirection, I2);
for j = 1:rows
    for i = 1:cols
        uw = reshape(Rflow(j, i, :), 2, 1);
        p = round([i; j]);
        pw = p + uw;
%         dnorm = normlize(pw - e2i);

        %unit vector pointing from epipole in I2 to current pw (current
        %pixel + rotational offset)
        dnorm = reshape(normlizeDirection(j, i, :), 2, 1);
%         assert(isequal(dnorm, reshape(normlizeDirection(j, i, :), 2, 1)));

        %evaluate disparity from 0 to maximun 
        for d = 0:dMax
            % convert disparity to 2D offset
            offset = d*dnorm;
            
            %target pixel position in I2 needs to evaluate cost
            pt = pw + offset; 
            
            pt = round(pt);
            pt(1) = min(cols, max(1, pt(1)));
            pt(2) = min(rows, max(1, pt(2)));

            cost = 0;
            for dy=-halfWinSize:halfWinSize
                for dx = -halfWinSize:halfWinSize
                    p2x = min(cols, max(1, pt(1) + dx));
                    p2y = min(rows, max(1, pt(2) + dy));
                    p1x = min(cols, max(1, i + dx));
                    p1y = min(rows, max(1, j + dy));
                    cost = cost + abs(I1(p1y, p1x) - I2(p2y, p2x));
                end
            end
            C(j, i, d+1) = cost;
%             diff = abs(patch1 - patch2);
%             C(j, i, d+1) = sum(diff(:));
       
        end
% 
%         if(i==584&&j==46)
%             imshow(uint8(0.5*I2 + 0.5*I1));
%             hold on;
%             plot(e2i(1),e2i(2),'go','MarkerFaceColor','r','MarkerSize',5);
%             hold on;
%             
%             plot(i,j,'bx');
%             hold on;
%             for d = 0:dMax
%                 offset = d*dnorm;
% 
%                 pt = pw + offset;
% 
%                 
%                 plot(pt(1),pt(2),'r+', 'MarkerSize', 2);
%                 hold on;
%             end
%         
%             figure;
%             plot(reshape(C(j, i, :), [], 1));
%         end
    end
end

if (debug)
    [~, Ind] = min(C, [], 3);

    flowT = (Ind-1).*normlizeDirection;

    flow = flowT + Rflow;

    Fc = flowToColor(flow);

    gtFlow = flow_read_kitti('C:\Users\megamusz\Desktop\KITTI\training\flow_occ\000004_10.png');
    gtFlowc = flowToColor(gtFlow(:,:, 1:2));

    flow(:,:,3) = ones(rows, cols);
    rflow = Rflow;
    rflow(:,:,3) = ones(rows, cols);
    eimage  =  flow_error_image(gtFlow, flow, [3; 0.05]);
    eimageR =  flow_error_image(gtFlow, rflow,[3; 0.05]);
    er = flow_error(gtFlow, rflow, [3;0.05]);
    ef = flow_error(gtFlow, flow,  [3;0.05]);
    imshow([Fc;gtFlowc]);
    figure;
    imshow([eimageR;eimage]);
    title(['outlier percentage (>3) before vs after:', num2str([er, ef])]);
    flow_field(flow, I1, gtFlow);
    flow_field(Rflow, I1, gtFlow);
end

end

function dnorm = normlize(d)
    if(size(d, 3) > 1)
        dnorm = d./sqrt(sum(d.^2, 3));
    else
        dnorm = d/sqrt(sum(d.^2));
    end
end
