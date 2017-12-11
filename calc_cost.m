function C = calc_cost(I1, I2, epipole, Rflow, gtFlow, dMax, halfWinSize)

debug = true;
debug_plot_search_line = true;

if debug
    close all;
end

if nargin < 6
    dMax = 64;
end

if nargin < 7
    halfWinSize = 2;
end

assert(isequal(size(I1), size(I2)));


cen1 = census(I1, halfWinSize);
cen2 = census(I2, halfWinSize);


rows = size(I1, 1);
cols = size(I1, 2);
numPixels = rows*cols;

% get the inhomogenous coordinate of epipole
% might need to deal with the case when epipole is at infinity
assert(epipole(3) ~= 0, 'epipole is at inifinity, not supported');
e2i = [epipole(1)/epipole(3); epipole(2)/epipole(3)];

C = 65535 * ones(rows, cols, dMax + 1);

P = zeros(rows, cols, 2);
[P(:,:,1), P(:,:,2)] = meshgrid(1:cols, 1:rows);
E2I = zeros(rows, cols, 2);

E2I(:,:,1) = repmat(e2i(1), rows, cols);
E2I(:,:,2) = repmat(e2i(2), rows, cols);

normlizeDirection = normlize(P + Rflow - E2I);
normDirectionFlat = zeros(2, numPixels);
normDirectionFlat(1, :) = reshape(normlizeDirection(:,:, 1), 1, []);
normDirectionFlat(2, :) = reshape(normlizeDirection(:,:, 2), 1, []);

cen1Flat = zeros(size(cen1, 3), numPixels);
cen2Flat = zeros(size(cen2, 3), numPixels);
CFlat = 65535*ones(dMax+1, numPixels);

for i = 1:size(cen1, 3)
    cen1Flat(i, :) = reshape(cen1(:,:,i), 1, []);
    cen2Flat(i, :) = reshape(cen2(:,:,i), 1, []);
end
assert(numel(cen1Flat) == numel(cen1));

tic;
d = 0:dMax;
for ind = 1:numPixels
    
    [j, i] = ind2sub([rows, cols], ind);
    uw = reshape(Rflow(j, i, :), 2, 1);
    p = round([i; j]);
    pw = p + uw;

    %unit vector pointing from epipole in I2 to current pw (current
    %pixel + rotational offset)

    dnorm = normDirectionFlat(:, ind);
    %evaluate disparity from 0 to maximun 
    
    offset = d.*repmat(dnorm, 1, dMax+1);
    pt = pw + offset;
    pt = round(pt);
    
    pt(1,:) = min(cols, max(1, pt(1,:)));
    pt(2,:) = min(rows, max(1, pt(2,:)));
    ind2 = sub2ind([rows, cols], pt(2,:), pt(1,:));

    cenCur = cen1Flat(:, ind);
    cenRef = cen2Flat(:, ind2);

    diff = cenCur ~= cenRef;

    cost = sum(diff);
    CFlat(:, ind) = cost'; 
          
    if(debug_plot_search_line)
        if(i ==halfWinSize+1 && j==halfWinSize+1)
            imshow(uint8(0.5*I2 + 0.5*I1));
            hold on;
            plot(e2i(1),e2i(2),'go','MarkerFaceColor','r','MarkerSize',5);
            hold on;
        end

        if(mod(i, 50) == 0 && mod(j, 20) == 0 && j > 150)

            gtVec = gtFlow(j, i,:);

            if(gtVec(3)==0)
                continue;
            end
            line([p(1), p(1)+gtVec(1)], [p(2), p(2) + gtVec(2)], 'Color', 'g');
            hold on;


            plot(i,j,'bx');
            hold on;
            line([p(1), pw(1)], [p(2), pw(2)], 'Color', 'y');
            hold on;

            line([pw(1), pw(1) + dMax*dnorm(1)], [pw(2), pw(2) + dMax*dnorm(2)], 'Color', 'r');

    %             figure;
    %             plot(reshape(C(j, i, :), [], 1));
    %             title(['cost at pixel: ', num2str([i, j])]);

            end
    end
  
end


toc;

for d = 1:dMax+1
    C(:,:, d) = reshape(CFlat(d, :), rows, cols);
end


if (debug)
    [minCost, Ind] = min(C, [], 3);

    flowT = (Ind-1).*normlizeDirection;

    flow = flowT + Rflow;
    flow(:,:,3) = ones(rows, cols);
    Fc = flow_to_color(flow);

    gtFlowc = flow_to_color(gtFlow);
    
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
