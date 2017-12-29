function [ imageFiltered, labelImage ] = speckle_filter( image, maxDiff, maxSpeckleSize )
%Implemented the speckle filter in OpenCV
%   [ imageFiltered ] = speckle_filter( image, maxDiff, maxSpeckleSize )
%   remove small blobs in image, the maxDiff is the threshold to
%   distinguish pixels into different segment; maxSpeckleSize is the
%   threshold to define small blob

if (nargin < 2)
    maxDiff = 2;
end

if(nargin < 3)
    maxSpeckleSize = 100;
end

height = size(image, 1);
width = size(image, 2);
labels = zeros(height*width, 1);
curLabel = 0;
regionTypes = zeros(1);
tic;

for y = 1:height
    for x =1:width
        ind = sub2ind([width, height], x, y);
        if( ~isnan(image(y, x)) ) 
            if(labels(ind) > 0)
                if(regionTypes(labels(ind)))
                    image(y, x) = nan;
                end
            else
%                  stack=java.util.Stack();
                 queue = Queue(width*height);
%                  stack.push(ind);
                queue.Enqueue(ind);
                 
                 curLabel = curLabel + 1;
                 regionTypes(curLabel) = false;
                 
                 labels(ind) = curLabel;
                 regionPixelNum = 0;
%                  while(~stack.empty())
                while(~queue.IsEmpty())
%                      curPixelIdx = stack.pop();
                     curPixelIdx = queue.Dequeue();
                     
                     [curx, cury] = ind2sub([width, height], curPixelIdx);
                     regionPixelNum = regionPixelNum+1;
                     
                     pixelValue = image(cury, curx);
                     
                     %check right neighbour
                     if(curx < width) 
                         if(labels(curPixelIdx+1) == 0 && ...
                             ~isnan(image(cury, curx+1)) && abs(pixelValue - image(cury, curx+1))< maxDiff)
                             labels(curPixelIdx + 1) = curLabel;
%                              stack.push(curPixelIdx+1);
                            queue.Enqueue(curPixelIdx+1);
                         end
                     end
                     
                     %check left neighbour
                     if(curx > 1) 
                         if(labels(curPixelIdx-1) == 0 && ...
                             ~isnan(image(cury, curx-1)) && abs(pixelValue - image(cury, curx-1))< maxDiff)
                             labels(curPixelIdx - 1) = curLabel;
%                              stack.push(curPixelIdx-1);
                            queue.Enqueue(curPixelIdx-1);
                         end
                     end
                     
                     %check bottom neighbour
                     if(cury < height)
                         if(labels(curPixelIdx + width) == 0 && ...
                             ~isnan(image(cury+1, curx)) && abs(pixelValue - image(cury+1, curx))< maxDiff)
                             labels(curPixelIdx + width) = curLabel;
%                              stack.push(curPixelIdx + width);
                            queue.Enqueue(curPixelIdx + width);
                         end
                     end
                     
                     %check top neighbour
                     if(cury > 1 )
                         if(labels(curPixelIdx - width) == 0 && ...
                             ~isnan(image(cury-1, curx)) && abs(pixelValue - image(cury-1, curx))< maxDiff)
                             labels(curPixelIdx - width) = curLabel;
%                              stack.push(curPixelIdx - width);
                            queue.Enqueue(curPixelIdx - width);
                         end
                     end
                 end
                 
                 if(regionPixelNum < maxSpeckleSize)
                     regionTypes(curLabel) = true;
                     image(y, x) = nan;
                 end
            end
        end
    end
end
labelImage = reshape(labels, [width, height])';
imageFiltered = image;
toc;
end

