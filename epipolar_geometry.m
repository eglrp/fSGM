function [F, E, H, epi, status] = epipolar_geometry(I1, I2, K)
%Calculate the epipolar geometry for two Camera views I1 and I2
% output is Fundamental Matrix F and epipoles in I2. 
%

% if camera intrinsic is avaible then the Esstential Matrix can also be
% calculated.

debug = 0;
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

%%% feature detection
points1 = detectSURFFeatures(I1, 'MetricThreshold', 500);
points2 = detectSURFFeatures(I2, 'MetricThreshold', 500);

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

%%% estimate fundamental matrix F
[F, inlierIndx, status] = estimateFundamentalMatrix(matchedPoints1,...
    matchedPoints2,'Method','LMedS',...
    'NumTrials',10000,'DistanceThreshold',1e-4);

disp(['Inlier percentage: ', num2str(sum(inlierIndx)/length(matchedPoints1))]);
if(status ~= 0)
    return;
end

if(debug)
    epiLines = epipolarLine(F',matchedPoints2.Location(inlierIndx, :));
    points = lineToBorderPoints(epiLines,size(I1));

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

fp = fopen('000000_10_fund.dat', 'rb');

F = fread(fp, 9, 'double');
F = reshape(F, [3, 3])';

%%% compute epipole in I2 %%%
% F'e' = 0, Fe = 0
% e & e' are epipole in I1 & I2 
[~, ~, V] = svd(F');
epiH = V(:, end);
epi = [epiH(1)/epiH(3); epiH(2)/epiH(3); 1];

%%% From F and K, recover E %%%
E = K' * F * K;
[U, ~, V] = svd(E);
W = [0 -1 0; 1 0 0 ; 0 0 1];
R1 = U*W *V';
R2 = U*W'*V';
if(det(R1) < 0) 
    R1 = -R1;
    R2 = -R2;
end

%%% select the correct rotation matrix
if(R1(1, 1) > 0 && R1(2, 2) > 0 && R1(3, 3) > 0)
    R = R1;
else
    R = R2;
end

%%% calculate the Homography to compensate the camera rotation. 
H = K*R/K;



