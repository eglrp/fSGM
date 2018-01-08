function [bestD, minC, mvSub, L1, L2, L3, L4] = sgm2d(C, m, n, P1, P2, subpixelRefine)
%Peform SGM on Cost Volume C

    if nargin < 4
        P1 = 7;  
    end

    if nargin < 5
        P2 = 100;
    end
      
    if nargin < 6
        subpixelRefine = 0;
    end

   
    rows = size(C, 1);
    cols = size(C, 2);
    dMax = size(C, 3);
    
    L1 = uint8(ones(size(C))); %path cost for left -> right direction
    L2 = uint8(permute(ones(size(C)), [2,1,3])); %path cost for upper->bottom direction
    L3 = uint8(permute(ones(size(C)), [2,1,3])); %path cost for bottom -> upper direction 
    L4 = uint8(ones(size(C))); %path cost for right -> left direction

    tic;
    Ct = permute(C, [2,1,3]);

    [lutR, lutC] = ind2sub([m, n], 1:dMax);
    
    for i=1:cols
        ir = cols+1-i;
        
        if(i == 1)
            L1(:, i, :) = C(:, i, :);
            L4(:, ir,:) = C(:, ir,:);

        else
% readable, but slow code
            minLpre = min(L1(:, i-1, :), [], 3);
            minLpre4 = min(L4(:, ir+1, :), [], 3);
     
            for d = 1:dMax
                r = lutR(d); c = lutC(d);
                left  = min(n, max(1, c - 1));
                right = min(n, max(1, c + 1));
                up    = min(m, max(1, r - 1));
                down  = min(m, max(1, r + 1));
                
                neighborIdx = sub2ind([m, n], [r, r, up, down],[left, right, c, c]);
                L1(:, i, d) = C(:, i, d) + min([L1(:, i-1, d), ...
                                                (min(L1(:, i-1, neighborIdx), [], 3) + P1), ...
                                                 minLpre + P2], [], 2) - minLpre;
                              

                
                L4(:, ir, d) = C(:, ir, d) + min([L4(:, ir+1, d), ...
                                                 (min(L4(:, ir+1, neighborIdx), [], 3) + P1) ...
                                                 minLpre4 + P2], [], 2) - minLpre4;
            end            
        end
    end
    
    
    for i=1:rows
        ir = rows+1-i;
        
        if(i == 1)
            L2(:, i, :) = Ct(:, i, :);
            L3(:, ir,:) = Ct(:, ir,:);

        else
            
% readable, but slow code

            minLpre2 = min(L2(:, i-1, :), [], 3);
            minLpre3 = min(L3(:, ir+1, :), [], 3);
            for d = 1:dMax

                r = lutR(d); c = lutC(d);
                left  = min(n, max(1, c - 1));
                right = min(n, max(1, c + 1));
                up    = min(m, max(1, r - 1));
                down  = min(m, max(1, r + 1));
                neighborIdx = sub2ind([m, n], [r, r, up, down],[left, right, c, c]);
                
                L2(:, i, d) = Ct(:, i, d) + min([L2(:, i-1, d), ...
                                                (min(L2(:, i-1, neighborIdx), [], 3) +  P1), ...
                                                 minLpre2 + P2], [], 2) - minLpre2;
                              
     
                L3(:, ir, d) = Ct(:, ir, d) + min([L3(:, ir+1, d), ...
                                                 (min(L3(:, ir+1, neighborIdx), [], 3) + P1), ...
                                                 minLpre3 + P2], [], 2) - minLpre3;
            end
        end
    end
% 
    toc;
    L2 = permute(L2, [2, 1, 3]);
    L3 = permute(L3, [2, 1, 3]);
    L = L1 + L2 + L3 + L4 ;
%     L = L1 + L4;
    [minC, bestD] = min(L, [], 3);
    
    mvSub = zeros(rows, cols, 2);
    if(subpixelRefine)
        % do subpixel quadratic interpolation:
        % fit parabola into (x1=d-1, y1=C[d-1]), (x2=d, y2=C[d]), (x3=d+1, y3=C[d+1])
        % then find minimum of the parabola.
        for j = 1:rows
            for i =1:cols

                [r, c] = ind2sub([m, n], bestD(j, i));
                c0 = double(L(j, i, bestD(j, i)));
                
                if(r > 1 && r < m)
                    cLeft = double(L(j, i, bestD(j, i)-1));
                    cRight  = double(L(j, i, bestD(j, i)+1));
                  
                    
                    if (cRight < cLeft)
                        mvSub(j, i, 2) = (cRight-cLeft)/(c0 - cLeft)/2.0;
                    else
                        mvSub(j, i, 2) = (cRight-cLeft)/(c0 - cRight)/2.0;
                    end
                end
                
                if(c>1 && c < n)
                    cLeft = double(L(j, i, bestD(j, i)-m));
                    cRight  = double(L(j, i, bestD(j, i)+m));
                    
                    if (cRight < cLeft)
                        mvSub(j, i, 1) = (cRight-cLeft)/(c0 - cLeft)/2.0;
                    else
                        mvSub(j, i, 1) = (cRight-cLeft)/(c0 - cRight)/2.0;
                    end
                end

            end
        end
    end
    
    

%     bestD = bestD - 1;

end
