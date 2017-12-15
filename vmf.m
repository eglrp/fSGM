function [ flowMed ] = vmf( flow )
%VMF: Vector median filter for optical flow post-processing
%   Detailed explanation goes here

flowMed = zeros(size(flow));

% placeholder here using sperate median filter for mvx/mvy
for c = 1:size(flow, 3)
    flowMed(:,:,c) = medfilt2(flow(:,:,c), [5 5]);
end


end

