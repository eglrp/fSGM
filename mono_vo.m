function [w, t] = mono_vo(I1, I2, f, center)

debug = 0;

if size(I1, 3) > 1
    I1gray = rgb2gray(I1);
else
    I1gray = I1;
end

if size(I2, 3) > 1
    I2gray = rgb2gray(I2);
else
    I2gray = I2;
end

imagePoints1 = detectSURFFeatures(I1gray);
imagePoints2 = detectSURFFeatures(I2gray);

features1 = extractFeatures(I1gray,imagePoints1, 'SURFSize', 128);
features2 = extractFeatures(I2gray,imagePoints2, 'SURFSize', 128);


indexPairs = matchFeatures(features1,features2);
matchedPoints1 = imagePoints1(indexPairs(:,1));
matchedPoints2 = imagePoints2(indexPairs(:,2));
if debug
    figure
    showMatchedFeatures(I1,I2,matchedPoints1,matchedPoints2);
    title('Putative Matches')
end


intrinsics = cameraIntrinsics(f,center,size(I1gray));

[E,inliers, status] = estimateEssentialMatrix(matchedPoints1,matchedPoints2,...
 intrinsics.CameraParameters);

assert(status == 0);

if(debug)
%     show inliers
    figure;
    imshow(I1);
    hold on;
    plot(matchedPoints1.Location(inliers, 1), matchedPoints1.Location(inliers, 2), 'g+');
end


[R, t] = relativeCameraPose(E,intrinsics.CameraParameters,matchedPoints1.Location(inliers, :),matchedPoints2.Location(inliers, :));
w = rotationMatrixToVector(R);