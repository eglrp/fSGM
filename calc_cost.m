function C = calc_cost(I1, I2, epipole, Rflow, dMax, halfWinSize)

debug = 1;

useSAD = false;

if debug
    close all;
end

if nargin < 5
    dMax = 32;
end

if nargin < 6
    halfWinSize = 2;
end

% 
% I1 = double(imread('000004_10.png'));
% I2 = double(imread('000004_11.png'));


assert(isequal(size(I1), size(I2)));

if ~useSAD
    cen1 = census(I1, 2);
    cen2 = census(I2, 2);
end

% e1 = [-0.959870232857247;-0.280438303284809;-0.001869257957789];
% e2 = [0.965222503953150;0.261423812450535;0.001762993545078];
% 
% F = [1.160024130294298e-08,1.539326689928344e-04,-0.023099942444453;-1.531504691961629e-04,3.313869025373902e-07,0.078593562935847;0.022703419892627,-0.084325838315150,0.992831445999651];
% P0 = [718.856000000000,0,607.192800000000,0;0,718.856000000000,185.215700000000,0;0,0,1,0];
% w = [-0.002409420424081,-0.046632321275709,-0.002414820939297];
% Rflow = rotation_motion(w, P0(1,1), [P0(1, 3) P0(2, 3)], [size(I1, 2), size(I1, 1)]);

rows = size(I1, 1);
cols = size(I1, 2);

% get the inhomogenous coordinate of epipole
% might need to deal with the case when epipole is at infinity
assert(epipole(3) ~= 0);
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
tic;

% validRows = length(1+halfWinSize:rows-halfWinSize);
% validCols = length(1+halfWinSize:cols-halfWinSize);
% numValidPixels = validRows*validCols;

% 
% parfor k = 1:numValidPixels
%     %convert to 2D position
%     j = halfWinSize + 1 + floor((k-1)/validCols);
%     i = halfWinSize + 1 + mod(k-1, validCols);
%     uw = reshape(Rflow(j, i, :), 2, 1);
%     p = round([i; j]);
%     pw = p + uw;
%     
%     dnorm = reshape(normlizeDirection(j, i, :), 2, 1);
%     
%     for d = 0:dMax
%         % convert disparity to 2D offset
%         offset = d*dnorm;
% 
%         %target pixel position in I2 needs to evaluate cost
%         pt = pw + offset; 
% 
%         pt = round(pt);
%         pt(1) = min(cols, max(1, pt(1)));
%         pt(2) = min(rows, max(1, pt(2)));
% 
%         cost = 0;
% 
% 
%         p2xMin = min(cols, max(1, pt(1) - halfWinSize));
%         p2yMin = min(rows, max(1, pt(2) - halfWinSize));
% 
%         p1xMin = i - halfWinSize;
%         p1yMin = j - halfWinSize;
% 
%         p2xMax = min(cols, max(1, pt(1) + halfWinSize));
%         p2yMax = min(rows, max(1, pt(2) + halfWinSize));
% 
% 
%         padR = pt(1) + halfWinSize - p2xMax;
%         padL = p2xMin - (pt(1) - halfWinSize);
%         padD = pt(2) + halfWinSize - p2yMax;
%         padU = p2yMin - (pt(2) - halfWinSize);
% 
%         p1xMax = i + halfWinSize;
%         p1yMax = j + halfWinSize;
% 
%         P1 = cen1(p1yMin+padU:p1yMax-padD, p1xMin+padL:p1xMax-padR, :);
%         P2 = cen2(p2yMin:p2yMax, p2xMin:p2xMax, :);
% 
% %         if(~isequal(size(P1), size(P2)))
% %             disp('size mismatch');
% %         end
%         diff = P1 ~= P2;
%         cost = sum(diff(:));
%                 
%     end
% %             C(j, i, d+1) = cost;
%    
%     
% end

for j = 1+halfWinSize:rows-halfWinSize
    for i = 1+halfWinSize:cols-halfWinSize
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
            if(useSAD)
                for dy=-halfWinSize:halfWinSize
                    for dx = -halfWinSize:halfWinSize
                        p2x = min(cols, max(1, pt(1) + dx));
                        p2y = min(rows, max(1, pt(2) + dy));
%                         p1x = min(cols, max(1, i + dx));
%                         p1y = min(rows, max(1, j + dy));
%                         p2x = pt(1) + dx;
%                         p2y = pt(2) + dy;
                        p1x = i + dx;
                        p1y = j + dy;
                    end
                end
            else
                p2xMin = min(cols, max(1, pt(1) - halfWinSize));
                p2yMin = min(rows, max(1, pt(2) - halfWinSize));
                
                p1xMin = i - halfWinSize;
                p1yMin = j - halfWinSize;
                
                p2xMax = min(cols, max(1, pt(1) + halfWinSize));
                p2yMax = min(rows, max(1, pt(2) + halfWinSize));
  
                padR = pt(1) + halfWinSize - p2xMax;
                padL = p2xMin - (pt(1) - halfWinSize);
                padD = pt(2) + halfWinSize - p2yMax;
                padU = p2yMin - (pt(2) - halfWinSize);
                
                p1xMax = i + halfWinSize;
                p1yMax = j + halfWinSize;
                
                P1 = cen1(p1yMin+padU:p1yMax-padD, p1xMin+padL:p1xMax-padR, :);
                P2 = cen2(p2yMin:p2yMax, p2xMin:p2xMax, :);

                diff = P1 ~= P2;
                cost = sum(diff(:));
                
            end
            
            C(j, i, d+1) = cost;
        end
% 
%         if(i ==halfWinSize+1 && j==halfWinSize+1)
%             imshow(uint8(0.5*I2 + 0.5*I1));
%             hold on;
%             plot(e2i(1),e2i(2),'go','MarkerFaceColor','r','MarkerSize',5);
%             hold on;
%             gtFlow = flow_read_kitti('C:\Users\megamusz\Desktop\KITTI\training\flow_occ\000004_10.png');
%         end
%         
%         if(mod(i, 50) == 0 && mod(j, 20) == 0 && j > 150)
% 
%             gtVec = gtFlow(j, i,:);
%             
%             if(gtVec(3)==0)
%                 continue;
%             end
%             line([p(1), p(1)+gtVec(1)], [p(2), p(2) + gtVec(2)], 'Color', 'g');
%             hold on;
%             
% 
%             plot(i,j,'bx');
%             hold on;
%             line([p(1), pw(1)], [p(2), pw(2)], 'Color', 'y');
%             hold on;
%                
%             line([pw(1), pw(1) + dMax*dnorm(1)], [pw(2), pw(2) + dMax*dnorm(2)], 'Color', 'r');
% %             figure;
% %             plot(reshape(C(j, i, :), [], 1));
%             
%         end
    end
end
toc;

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
    [er, aeper] = flow_error(gtFlow, rflow, [3;0.05]);
    [ef, aepef] = flow_error(gtFlow, flow,  [3;0.05]);
    imshow([Fc;gtFlowc]);
    figure;
    imshow([eimageR;eimage]);
    title(['outlier percentage (>3) before vs after:', num2str([er, ef]), ', EPE: ', num2str([aeper, aepef])]);
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
