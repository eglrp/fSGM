function F = flow_field(flow, I)


if nargin < 2
    I = zeros(size(flow, 1), size(flow, 2));
end

h = size(I, 1);
w = size(I, 2);
imshow(I);
hold on;

step = 5;

[xx, yy] = meshgrid(1:step:w, 1:step:h);

mvx = flow(1:step:h, 1:step:w, 1);
mvy = flow(1:step:h, 1:step:w, 2);

quiver(xx, yy, mvx, mvy, 'Color', 'g');


end