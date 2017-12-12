function [bestD, L1, L2, L3, L4] = sgm(C)
%Peform SGM on Cost Volume C

    P1 = 7;  
    P2 = 100;  
    
    rows = size(C, 1);
    cols = size(C, 2);
    dMax = size(C, 3);
    
    L1 = zeros(size(C)); %path cost for left -> right direction
    L2 = permute(zeros(size(C)), [2,1,3]); %path cost for upper->bottom direction
    L3 = permute(zeros(size(C)), [2,1,3]); %path cost for bottom -> upper direction 
    L4 = zeros(size(C)); %path cost for right -> left direction
    
    tic;
    Ct = permute(C, [2,1,3]);
    for i=1:cols
        ir = cols+1-i;
        
        if(i == 1)
            L1(:, i, :) = C(:, i, :);
            L4(:, ir,:) = C(:, ir,:);
        else
%% readable, but slow code
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
%% fast, but not that reable code
            %%%%%%%%% L1
            minLpre = min(L1(:, i-1, :), [], 3);
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
        else
%% readable, but slow code
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
%% fast, but not that reable code
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
    [~, bestD] = min(L1 + L2 + L3 + L4, [], 3);

end
