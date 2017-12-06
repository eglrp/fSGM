function [] = flow_field(flow, I, gtFlow)


if nargin < 2
    I = zeros(size(flow, 1), size(flow, 2));
end


if ~isinteger(I)
    I = uint8(I);
end

h = size(I, 1);
w = size(I, 2);

step = 20;

[xx, yy] = meshgrid(5:step:w-5, 5:step:h-5);

figure;
imshow(I);

% plot ground truth mv if have
if nargin > 2
    mvx = gtFlow(5:step:end-5, 5:step:end-5, 1);
    mvy = gtFlow(5:step:end-5, 5:step:end-5, 2);
    hold on;
    quiver(xx, yy, mvx, mvy, 'AutoScale', 'off', 'Color', 'r', 'MaxHeadSize', 0.02);
    
end

% plot estimated mv
mvx = flow(5:step:end-5, 5:step:end-5, 1);
mvy = flow(5:step:end-5, 5:step:end-5, 2);
hold on;
quiver(xx, yy, mvx, mvy, 'AutoScale', 'off', 'Color', 'g', 'MaxHeadSize', 0.02);

end