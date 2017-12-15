function [ mv ] = rotation_motion( H, F, rows, cols )
%ROTATION_MOTION Summary of this function goes here
%   Detailed explanation goes here


%%% general motion
[xx, yy] = meshgrid(1:cols, 1:rows);
mv = zeros(rows, cols, 2);

P0 = zeros(3, cols*rows);
P0(1,:) = reshape(xx-1, 1, []);
P0(2,:) = reshape(yy-1, 1, []);
P0(3,:) = 1;

L2 = computeEpipoleLineI2(F, P0);
% double coefficientToEpipolarLine = -epipolarLine(0)*rotatedPoint(0) - epipolarLine(1)*rotatedPoint(1) - epipolarLine(2);
% 			//coefficientToEpipolarLine = 0;
% 			double noDisparityPointX = rotatedPoint(0) + epipolarLine(0)*coefficientToEpipolarLine;
% 			double noDisparityPointY = rotatedPoint(1) + epipolarLine(1)*coefficientToEpipolarLine;

P1 = H*P0;
P1 = P1./P1(3,:);
OFF = P1(1:2, :) - P0(1:2, :);
for j = 1:rows
    for i = 1:cols
        ind = sub2ind([rows, cols], j, i);
        coefficientToEpipolarLine = -L2(:, ind)'*P1(:, ind);
        OFF(:, ind) = OFF(:, ind) + coefficientToEpipolarLine*L2(1:2, ind);
    end
end


mv(:,:,1) = reshape(OFF(1,:), rows, cols);
mv(:,:,2) = reshape(OFF(2,:), rows, cols);



%%% small motion assumption
% x_ = xx - center(1);
% y_ = yy - center(2);
% xy_ = x_.*y_;
% 
% mv(:,:,1) = -( f*w(2) - w(3)*y_ + w(2)/f*x_.^2 - w(1)/f*xy_);
% mv(:,:,2) = -(f*w(1) + w(3)*x_ + w(2)/f*xy_   - w(1)/f*y_.^2);

end

function l2 = computeEpipoleLineI2(F, p1)
    l2 = F*p1;
    normlizeFactor = sqrt(l2(1,:).^2 + l2(2,:).^2);
    normlizeFactor(normlizeFactor < 1e-6) = 1.0;
    l2 = l2./normlizeFactor;
end
