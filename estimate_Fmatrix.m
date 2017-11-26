close all;
I1 = rgb2gray(imread('000000_10.png'));
I2 = rgb2gray(imread('000000_11.png'));


points1 = detectSURFFeatures(I1);
points2 = detectSURFFeatures(I2);
% FAST features have very poor performance in the Fundermental matrix
% points1 = detectFASTFeatures(I1);
% points2 = detectFASTFeatures(I2);

[f1, vpts1] = extractFeatures(I1, points1);
[f2, vpts2] = extractFeatures(I2, points2);

indexPairs = matchFeatures(f1,f2) ;
matchedPoints1 = vpts1(indexPairs(:,1));
matchedPoints2 = vpts2(indexPairs(:,2));


% figure; showMatchedFeatures(I1,I2,matchedPoints1,matchedPoints2);
% legend('matched points 1','matched points 2');


[F, inlierIndx, status] = estimateFundamentalMatrix(matchedPoints1,...
    matchedPoints2,'Method','RANSAC',...
    'NumTrials',2000,'DistanceThreshold',1e-4);


epiLines = epipolarLine(F',matchedPoints2.Location(inlierIndx, :));
points = lineToBorderPoints(epiLines,size(I1));

imshow(imread('000000_10.png'));
hold on;
plot(matchedPoints1.Location(inlierIndx,1), matchedPoints1.Location(inlierIndx, 2), 'go');
hold on;
line(points(:,[1,3])',points(:,[2,4])');

% calculate the epipoles
[isIn, epipole] = isEpipoleInImage(F, size(I1));
epipole

