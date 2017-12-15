function [bestD, minC, L1, L2, L3, L4] = sgm(C, P1, P2, e)
%Peform SGM on Cost Volume C

    if nargin < 2
        P1 = 7;  
    end

    if nargin < 3
        P2 = 100;
    end
    
    subpixelRefine = true;
    adpativeP2 = false;
    if(adpativeP2)
        downScaleFactor = 1;
    end
    enableDiagonalPath = false;
    
    P2C = P2;
    
    rows = size(C, 1);
    cols = size(C, 2);
    dMax = size(C, 3);
    
    L1 = 65535*ones(size(C)); %path cost for left -> right direction
    L2 = permute(65535*ones(size(C)), [2,1,3]); %path cost for upper->bottom direction
    L3 = permute(65535*ones(size(C)), [2,1,3]); %path cost for bottom -> upper direction 
    L4 = 65535*ones(size(C)); %path cost for right -> left direction
    if(enableDiagonalPath)
        L5 = 65535*ones(size(C)); %path cost for 45 degree top-left to bottom-right
        L6 = 65535*ones(size(C)); %path cost for 45 degree top-right to bottom-left
        L7 = permute(65535*ones(size(C)), [2, 1, 3]); %path cost for 45 degree bottom-left to top-right
        L8 = permute(65535*ones(size(C)), [2, 1, 3]); %path cost for 45 degree bottom-right to top-left
    end
    tic;
    Ct = permute(C, [2,1,3]);
    for i=1:cols
        ir = cols+1-i;
        
        if(i == 1)
            L1(:, i, :) = C(:, i, :);
            L4(:, ir,:) = C(:, ir,:);
            if(enableDiagonalPath)
                L5(:, i, :) = C(:, i, :);
                L6(:, ir,:) = C(:, ir,:);
            end
        else
            if(enableDiagonalPath)
                L5(1, i, :) = C(1, i, :);
                minLpre = min(L5(1:rows-1, i-1, :), [], 3);
                if(adpativeP2)
                    P2 = P2C ./ (1 + downScaleFactor*e(2:rows, i));
                end
                L5(2:rows, i, 1) = C(2:rows, i, 1) + min([L5(1:rows-1, i-1,  1), ...
                                                (L5(1:rows-1, i-1, 2) + P1), ...
                                                (L5(1:rows-1, i-1, 1) + P1), ...
                                                minLpre + P2], [], 2) - minLpre;
                d = 2:dMax-1;
                minLpre2 = repmat(minLpre, 1, 1, length(d));
                L5(2:rows, i, d) = C(2:rows, i, d) + min([L5(1:rows-1, i-1,  d), ...
                                                (L5(1:rows-1, i-1, d+1) + P1), ...
                                                (L5(1:rows-1, i-1, d-1) + P1), ...
                                                minLpre2 + P2], [], 2) - minLpre2;

                L5(2:rows, i, dMax) = C(2:rows, i, dMax) + min([L5(1:rows-1, i-1,  dMax), ...
                                                (L5(1:rows-1, i-1, dMax) + P1), ...
                                                (L5(1:rows-1, i-1, dMax-1) + P1), ...
                                                minLpre + P2], [], 2) - minLpre;

                if(adpativeP2)
                    P2 = P2C ./ (1 + downScaleFactor*e(2:rows, ir));
                end

                L6(1, ir, :) = C(1, ir, :);
                minLpre = min(L6(1:rows-1, ir+1, :), [], 3);

                L6(2:rows, ir, 1) = C(2:rows, ir, 1) + min([L6(1:rows-1, ir+1,  1), ...
                                                (L6(1:rows-1, ir+1, 2) + P1), ...
                                                (L6(1:rows-1, ir+1, 1) + P1), ...
                                                minLpre + P2], [], 2) - minLpre;
                d = 2:dMax-1;
                minLpre2 = repmat(minLpre, 1, 1, length(d));
                L6(2:rows, ir, d) = C(2:rows, ir, d) + min([L6(1:rows-1, ir+1,  d), ...
                                                (L6(1:rows-1, ir+1, d+1) + P1), ...
                                                (L6(1:rows-1, ir+1, d-1) + P1), ...
                                                minLpre2 + P2], [], 2) - minLpre2;

                L6(2:rows, ir, dMax) = C(2:rows, ir, dMax) + min([L6(1:rows-1, ir+1,  dMax), ...
                                                (L6(1:rows-1, ir+1, dMax) + P1), ...
                                                (L6(1:rows-1, ir+1, dMax-1) + P1), ...
                                                minLpre + P2], [], 2) - minLpre;
            end
% readable, but slow code
%             for d = 1:dMax
% 
%                 minLpre = min(L1(:, i-1, :), [], 3);
%                 L1(:, i, d) = C(:, i, d) + min([L1(:, i-1, d), ...
%                                                 (L1(:, i-1, min(dMax, d+1)) + P1), ...
%                                                 (L1(:, i-1, max(1, d-1)) + P1), ...
%                                                  minLpre + P2], [], 2) - minLpre;
%                               
% 
%                 minLpre = min(L4(:, ir+1, :), [], 3);
%                 L4(:, ir, d) = C(:, ir, d) + min([L4(:, ir+1, d), ...
%                                                  (L4(:, ir+1, min(dMax, d+1)) + P1), ...
%                                                  (L4(:, ir+1, max(1, d-1)) + P1), ...
%                                                  minLpre + P2], [], 2) - minLpre;
%             end
% fast, but not that reable code
            %%%%%%%%% L1
            minLpre = min(L1(:, i-1, :), [], 3);
            if(adpativeP2)
                P2 = P2C ./ (1 + downScaleFactor*e(:, i));
            end
            % d = 1
            L1(:, i, 1) = C(:, i, 1) + min([L1(:, i-1,  1), ...
                                            (L1(:, i-1, 2) + P1), ...
                                            (L1(:, i-1, 1) + P1), ...
                                             minLpre + P2], [], 2) - minLpre;
            % d = 2~ dMax-1
            d = 2:dMax-1;
            
            minLpre2 = repmat(minLpre, 1, 1, length(d));
            L1(:, i, d) = C(:, i, d) + min([L1(:, i-1, d), ...
                                            (L1(:, i-1, d+1) + P1), ...
                                            (L1(:, i-1, d-1) + P1), ...
                                             minLpre2 + P2], [], 2) - minLpre2;
                                         
            % d = dMax
            L1(:, i, dMax) = C(:, i, dMax) + min([L1(:, i-1,  dMax), ...
                                            (L1(:, i-1, dMax) + P1), ...
                                            (L1(:, i-1, dMax-1) + P1), ...
                                             minLpre + P2], [], 2) - minLpre;
                                         

            
            %%%%%%%%% L4       
            if(adpativeP2)
                P2 = P2C ./ (1 + downScaleFactor*e(:, ir));
            end
            minLpre = min(L4(:, ir+1, :), [], 3);
            % d = 1
            L4(:, ir, 1) = C(:, ir, 1) + min([L4(:, ir+1, 1), ...
                                             (L4(:, ir+1, 2) + P1), ...
                                             (L4(:, ir+1, 1) + P1), ...
                                             minLpre + P2], [], 2) - minLpre;
            % d = 2~ dMax-1                            
            minLpre2 = repmat(minLpre, 1, 1, length(d));
            L4(:, ir, d) = C(:, ir, d) + min([L4(:, ir+1, d), ...
                                             (L4(:, ir+1, d+1) + P1), ...
                                             (L4(:, ir+1, d-1) + P1), ...
                                             minLpre2 + P2], [], 2) - minLpre2;

            % d = dMax
            L4(:, ir, dMax) = C(:, ir, dMax) + min([L4(:, ir+1, dMax), ...
                                             (L4(:, ir+1, dMax) + P1), ...
                                             (L4(:, ir+1, dMax-1) + P1), ...
                                             minLpre + P2], [], 2) - minLpre;
                                         
        end
    end
    
    
    for i=1:rows
        ir = rows+1-i;
        
        if(i == 1)
            L2(:, i, :) = Ct(:, i, :);
            L3(:, ir,:) = Ct(:, ir,:);
            if(enableDiagonalPath)
                L7(:, i, :) = Ct(:, i, :);
                L8(:, ir,:) = Ct(:, ir,:);
            end
        else
            if(enableDiagonalPath)
                L7(1, i, :) = Ct(1, i, :);
                minLpre = min(L7(1:cols-1, i-1, :), [], 3);
                if(adpativeP2)
                    P2 = P2C ./ (1 + downScaleFactor*e(2:cols, i));
                end
                L7(2:cols, i, 1) = Ct(2:cols, i, 1) + min([L7(1:cols-1, i-1,  1), ...
                                                (L7(1:cols-1, i-1, 2) + P1), ...
                                                (L7(1:cols-1, i-1, 1) + P1), ...
                                                minLpre + P2], [], 2) - minLpre;
                d = 2:dMax-1;
                minLpre2 = repmat(minLpre, 1, 1, length(d));
                L7(2:cols, i, d) = Ct(2:cols, i, d) + min([L7(1:cols-1, i-1,  d), ...
                                                (L7(1:cols-1, i-1, d+1) + P1), ...
                                                (L7(1:cols-1, i-1, d-1) + P1), ...
                                                minLpre2 + P2], [], 2) - minLpre2;

                L7(2:cols, i, dMax) = Ct(2:cols, i, dMax) + min([L7(1:cols-1, i-1,  dMax), ...
                                                (L7(1:cols-1, i-1, dMax) + P1), ...
                                                (L7(1:cols-1, i-1, dMax-1) + P1), ...
                                                minLpre + P2], [], 2) - minLpre;

                if(adpativeP2)
                    P2 = P2C ./ (1 + downScaleFactor*e(2:cols, ir));
                end

                L8(1, ir, :) = Ct(1, ir, :);
                minLpre = min(L8(1:cols-1, ir+1, :), [], 3);

                L8(2:cols, ir, 1) = Ct(2:cols, ir, 1) + min([L8(1:cols-1, ir+1,  1), ...
                                                (L8(1:cols-1, ir+1, 2) + P1), ...
                                                (L8(1:cols-1, ir+1, 1) + P1), ...
                                                minLpre + P2], [], 2) - minLpre;
                d = 2:dMax-1;
                minLpre2 = repmat(minLpre, 1, 1, length(d));
                L8(2:cols, ir, d) = Ct(2:cols, ir, d) + min([L8(1:cols-1, ir+1,  d), ...
                                                (L8(1:cols-1, ir+1, d+1) + P1), ...
                                                (L8(1:cols-1, ir+1, d-1) + P1), ...
                                                minLpre2 + P2], [], 2) - minLpre2;

                L8(2:cols, ir, dMax) = Ct(2:cols, ir, dMax) + min([L8(1:cols-1, ir+1,  dMax), ...
                                                (L8(1:cols-1, ir+1, dMax) + P1), ...
                                                (L8(1:cols-1, ir+1, dMax-1) + P1), ...
                                                minLpre + P2], [], 2) - minLpre;
            end
% readable, but slow code
%             for d = 1:dMax
% 
%                 minLpre = min(L2(:, i-1, :), [], 3);
%                 L2(:, i, d) = Ct(:, i, d) + min([L2(:, i-1, d), ...
%                                                 (L2(:, i-1, min(dMax, d+1)) + P1), ...
%                                                 (L2(:, i-1, max(1, d-1)) + P1), ...
%                                                  minLpre + P2], [], 2) - minLpre;
%                               
%      
%                 minLpre = min(L3(:, ir+1, :), [], 3);
%                 L3(:, ir, d) = Ct(:, ir, d) + min([L3(:, ir+1, d), ...
%                                                  (L3(:, ir+1, min(dMax, d+1)) + P1), ...
%                                                  (L3(:, ir+1, max(1, d-1)) + P1), ...
%                                                  minLpre + P2], [], 2) - minLpre;
%             end
% fast, but not that reable code
            if(adpativeP2)
                P2 = P2C ./ (1 + downScaleFactor*e(i, :)');
            end
            minLpre = min(L2(:, i-1, :), [], 3);
            % d = 1
            L2(:, i, 1) = Ct(:, i, 1) + min([L2(:, i-1, 1), ...
                                            (L2(:, i-1, 2) + P1), ...
                                            (L2(:, i-1, 1) + P1), ...
                                             minLpre + P2], [], 2) - minLpre;

            % d = 2~dMax-1                          
            d = 2:dMax-1;
            minLpre2 = repmat(minLpre, 1, 1, length(d));
            L2(:, i, d) = Ct(:, i, d) + min([L2(:, i-1, d), ...
                                            (L2(:, i-1, d+1) + P1), ...
                                            (L2(:, i-1, d-1) + P1), ...
                                             minLpre2 + P2], [], 2) - minLpre2;

            % d = dMax
            L2(:, i, dMax) = Ct(:, i, dMax) + min([L2(:, i-1, dMax), ...
                                            (L2(:, i-1, dMax) + P1), ...
                                            (L2(:, i-1, dMax-1) + P1), ...
                                             minLpre + P2], [], 2) - minLpre;
            
            
            if(adpativeP2)
                P2 = P2C ./ (1 + downScaleFactor*e(ir, :)');
            end
            minLpre = min(L3(:, ir+1, :), [], 3);
            % d = 1
            L3(:, ir, 1) = Ct(:, ir, 1) + min([L3(:, ir+1, 1), ...
                                             (L3(:, ir+1,  2) + P1), ...
                                             (L3(:, ir+1,  1) + P1), ...
                                             minLpre + P2], [], 2) - minLpre;
            % d = 2~dMax-1
            minLpre2 = repmat(minLpre, 1, 1, length(d));
            L3(:, ir, d) = Ct(:, ir, d) + min([L3(:, ir+1, d), ...
                                             (L3(:, ir+1, d+1) + P1), ...
                                             (L3(:, ir+1, d-1) + P1), ...
                                             minLpre2 + P2], [], 2) - minLpre2;
            % d = dMax
            L3(:, ir, dMax) = Ct(:, ir, dMax) + min([L3(:, ir+1, dMax), ...
                                             (L3(:, ir+1,  dMax) + P1), ...
                                             (L3(:, ir+1,  dMax-1) + P1), ...
                                             minLpre + P2], [], 2) - minLpre;
        end
    end

    toc;
    L2 = permute(L2, [2, 1, 3]);
    L3 = permute(L3, [2, 1, 3]);
    L = L1 + L2 + L3 + L4 ;
    if(enableDiagonalPath)
        L7 = permute(L7, [2, 1, 3]);
        L8 = permute(L8, [2, 1, 3]);
        L = L + L5 + L6 + L7 + L8;
    end
%     L = L5 + L6;
    
    [minC, bestD] = min(L, [], 3);
    
    
if(subpixelRefine)
    % do subpixel quadratic interpolation:
    % fit parabola into (x1=d-1, y1=C[d-1]), (x2=d, y2=C[d]), (x3=d+1, y3=C[d+1])
    % then find minimum of the parabola.
    for j = 1:rows
        for i =1:cols
            
            if (bestD(j, i) > 1 && bestD(j, i) < dMax)
                c_1 = L(j, i, bestD(j, i)-1);
                c = L(j, i, bestD(j, i));
                c1 = L(j, i, bestD(j, i)+1);
                denorm2 = max(c_1+c1-2*c, 1);
                bestD(j, i) = bestD(j, i) + ((c_1-c1) + denorm2)/(denorm2*2);
            end
        end
    end
end

end
