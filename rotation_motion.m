function [ mv ] = rotation_motion( w, f, center, Isize )
%ROTATION_MOTION Summary of this function goes here
%   Detailed explanation goes here
[xx, yy] = meshgrid(1:Isize(1), 1:Isize(2));
mv = zeros(Isize(2), Isize(1), 2);

x_ = xx - center(1);
y_ = yy - center(2);
xy_ = x_.*y_;

mv(:,:,1) = -( f*w(2) - w(3)*y_ + w(2)/f*x_.^2 - w(1)/f*xy_);
mv(:,:,2) = -(f*w(1) + w(3)*x_ + w(2)/f*xy_   - w(1)/f*y_.^2);

end

