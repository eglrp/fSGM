function [F, E, H, epi, direction, Pd0, normlizeDirection, O,  Rflow, status] = epipolar_geometry(I1, I2, K)
%Calculate the epipolar geometry for two Camera views I1 and I2
% output is Fundamental Matrix F and esstential Matrix E
% epipole position in I2
% direction
%
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
    cols = size(I1, 2);
    rows = size(I1, 1);

%%% calculate the Fundamental matrix
    [F, matchedPoints1, matchedPoints2, inlierIndx, status] = estimate_fundamental_matrix(I1, I2);
    
    %in case Fundamental matrix estimation fail, set all output to zero
    if(status) 
        F = 0; E = 0; H = 0; epi = 0; direction = 0; Pd0 = 0; normlizeDirection = 0;O = 0; Rflow = 0;
        return;
    end
% fp = fopen('000000_10_fund.dat', 'rb');
% 
% F = fread(fp, 9, 'double');
% F = reshape(F, [3, 3])';

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


%%% detect direction (expansion or contraction)
    expansion = 0;
    inlierNum = sum(inlierIndx);

    for i = 1:length(inlierIndx)
        if(inlierIndx(i)) 
            x1 = matchedPoints1.Location(i,1);
            y1 = matchedPoints1.Location(i,2);
            pr = H*[x1, y1, 1]';
            xr = pr(1)/pr(3);
            yr = pr(2)/pr(3);

            dist1 = sqrt((xr-epi(1))^2 + (yr-epi(2))^2);

            x2 = matchedPoints2.Location(i,1);
            y2 = matchedPoints2.Location(i,2);

            dist2 = sqrt((x2-epi(1))^2 + (y2-epi(2))^2);

            if(dist2 > dist1)
                expansion = expansion + 1;
            end
        end
    end

    if(expansion/inlierNum > 0.5)
        direction = 0;
    else
        direction = 1;
    end


%%% calculate starting search postition in I2
% and the direction along epipolar line for every pixel in I1

    P = zeros(rows, cols, 2);
    [P(:,:,1), P(:,:,2)] = meshgrid(1:cols, 1:rows);
    E2I = zeros(rows, cols, 2);
    E2I(:,:,1) = repmat(epi(1), rows, cols);
    E2I(:,:,2) = repmat(epi(2), rows, cols);

    Rflow = rotation_motion(H, F, rows, cols);

    Pd0 = P + Rflow; % position with zero disparity
    Direct = Pd0 - E2I;
    if(direction)
        Direct = -Direct;
    end

    O = sqrt(sum(Direct.^2, 3)); %Offset from Pd0 to epipole in I2
    normlizeDirection = normlize(Direct); %normlized direction
end

function dnorm = normlize(d)
    if(size(d, 3) > 1)
        dnorm = d./sqrt(sum(d.^2, 3));
    else
        dnorm = d/sqrt(sum(d.^2));
    end
end

function [F, matchedPoints1, matchedPoints2, inlierIndx, status] = estimate_fundamental_matrix(I1, I2)
    debug = 0;
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
        'NumTrials',10000,'DistanceThreshold',1e-4, 'ReportRuntimeError', false);
    
    if(debug)
        disp(['Inlier percentage: ', num2str(sum(inlierIndx)/length(matchedPoints1))]);

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
end