function [F, status] = estimate_Fmatrix(I1, I2)
% close all;
% I1 = rgb2gray(imread('000000_10.png'));
% I2 = rgb2gray(imread('000000_11.png'));
debug = 1;
if(debug)
    close all;
end

assert(isequal(size(I1), size(I2)));

if(size(I1, 3) > 1)
    I1 = rgb2gray(I1);
end

if(size(I2, 3) > 1)
    I2 = rgb2gray(I2);
end

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

if(debug)
    figure; showMatchedFeatures(I1,I2,matchedPoints1,matchedPoints2);
    legend('matched points 1','matched points 2');
    title('all matches');
    saveas(gcf, 'all_matches.png');
end

[F, inlierIndx, status] = estimateFundamentalMatrix(matchedPoints1,...
    matchedPoints2,'Method','LMedS',...
    'NumTrials',2000,'DistanceThreshold',1e-4);

if(status ~= 0)
    return;
end

epiLines = epipolarLine(F',matchedPoints2.Location(inlierIndx, :));
points = lineToBorderPoints(epiLines,size(I1));

if(debug)
    figure; showMatchedFeatures(I1,I2,matchedPoints1(inlierIndx),matchedPoints2(inlierIndx));
    legend('matched points 1','matched points 2');
    title('inliers');
    saveas(gcf, 'inliers_matches.png');
    
    imshow(I1);
    hold on;
    plot(matchedPoints1.Location(inlierIndx,1), matchedPoints1.Location(inlierIndx, 2), 'g+');
    hold on;
    line(points(:,[1,3])',points(:,[2,4])');
    title('epiloar lines');
    saveas(gcf, 'epipolar_lines.png');
end

% calculate the epipoles
% [isIn, epipole] = isEpipoleInImage(F, size(I1));
% epipole

