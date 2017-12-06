function bestD = sgm(C)
%Peform SGM on Cost Volume C

%     P1 = 10;
%     P2 = 100;
    P1 = 200;  
    P2 = 3000;  
    
    rows = size(C, 1);
    cols = size(C, 2);
    dMax = size(C, 3);
    
    L1 = zeros(size(C)); %path cost for left -> right direction
    L2 = permute(zeros(size(C)), [2,1,3]); %path cost for upper->bottom direction
    L3 = permute(zeros(size(C)), [2,1,3]); %path cost for bottom -> upper direction
%     L2 = zeros(size(C)); %path cost for upper->bottom direction
%     L3 = zeros(size(C)); %path cost for bottom -> upper direction
    
    L4 = zeros(size(C)); %path cost for right -> left direction
    tic;
    Ct = permute(C, [2,1,3]);
   
    
    
    for i=1:cols
        ir = cols+1-i;
        
        if(i == 1)
            L1(:, i, :) = C(:, i, :);
            L4(:, ir,:) = C(:, ir,:);
        else
            for d = 1:dMax
                minLeft  = min(L1(:, i-1, 1:max(1, d-2)), [], 3);
                minRight = min(L1(:, i-1, min(dMax, d+2):end), [], 3);


                L1(:, i, d) = C(:, i, d) + min([L1(:, i-1, d), ...
                                                (L1(:, i-1, min(dMax, d+1)) + P1), ...
                                                (L1(:, i-1, max(1, d-1)) + P1), ...
                                                (min(minLeft, minRight) + P2)], [], 2) - min(L1(:, i-1, :), [], 3);
                              
                                            
                minLeft  = min(L4(:, ir+1, 1:max(1, d-2)), [], 3);
                minRight = min(L4(:, ir+1, min(dMax, d+2):end), [], 3);


                L4(:, ir, d) = C(:, ir, d) + min([L4(:, ir+1, d), ...
                                                 (L4(:, ir+1, min(dMax, d+1)) + P1), ...
                                                 (L4(:, ir+1, max(1, d-1)) + P1), ...
                                                 (min(minLeft, minRight) + P2)], [], 2) - min(L4(:, ir+1, :), [], 3);
            end

        end
    end
    
    
    for i=1:rows
        ir = rows+1-i;
        
        if(i == 1)
            L2(:, i, :) = Ct(:, i, :);
            L3(:, ir,:) = Ct(:, ir,:);
        else
            for d = 1:dMax
                minLeft  = min(L2(:, i-1, 1:max(1, d-2)), [], 3);
                minRight = min(L2(:, i-1, min(dMax, d+2):end), [], 3);


                L2(:, i, d) = Ct(:, i, d) + min([L2(:, i-1, d), ...
                                                (L2(:, i-1, min(dMax, d+1)) + P1), ...
                                                (L2(:, i-1, max(1, d-1)) + P1), ...
                                                (min(minLeft, minRight) + P2)], [], 2) - min(L2(:, i-1, :), [], 3);
                              
                                            
                minLeft  = min(L3(:, ir+1, 1:max(1, d-2)), [], 3);
                minRight = min(L3(:, ir+1, min(dMax, d+2):end), [], 3);


                L3(:, ir, d) = Ct(:, ir, d) + min([L3(:, ir+1, d), ...
                                                 (L3(:, ir+1, min(dMax, d+1)) + P1), ...
                                                 (L3(:, ir+1, max(1, d-1)) + P1), ...
                                                 (min(minLeft, minRight) + P2)], [], 2) - min(L3(:, ir+1, :), [], 3);
            end

        end
    end
%     
    
%     
%     for j=1:rows
%         
%         if(j == 1)
%             L2(1, :, :) = C(1, :, :);
%         else
%             for d = 1:dMax
%                 minLeft  = min(L2(j-1, :, 1:max(1, d-2)), [], 3);
%                 minRight = min(L2(j-1, :, min(dMax, d+2):end), [], 3);
% 
% 
%                 L2(j, :, d) = C(j, :, d) + min([L2(j-1, :, d); ...
%                                                 L2(j-1, :, min(dMax, d+1)) + P1; ...
%                                                 L2(j-1, :, max(1, d-1)) + P1; ...
%                                                 min([minLeft; minRight]) + P2], [], 1) - min(L2(j-1, :, :), [], 3);
%             end
%         end
%     end
    
%     for j = 1:rows
%         for i = 1:cols
% %             initialize the Path cost L2 at the top bounary
%             if(j == 1)
%                 L2(1, i, :) = C(1, i, :);
%             else
%                 for d = 1:dMax
%                     minLeft  = min(L2(j-1, i, 1:max(1, d-2)), [], 3);
%                     minRight = min(L2(j-1, i, min(dMax, d+2):end), [], 3);
%                     
%                     
%                     L2(j, i, d) = C(j, i, d) + min([L2(j-1, i, d), ...
%                                                     L2(j-1, i, min(dMax, d+1)) + P1, ...
%                                                     L2(j-1, i, max(1, d-1)) + P1, ...
%                                                     min(minLeft, minRight) + P2]) - min(L2(j-1, i, :), [], 3);
%                 end
%             end
%    
% %             initialize the Path cost L1 at the left bounary
%             if(i == 1)
%                 L1(j, 1, :) = C(j, 1, :);
%                 
%             else
%                 for d = 1:dMax
%                     minLeft  = min(L1(j, i-1, 1:max(1, d-2)), [], 3);
%                     minRight = min(L1(j, i-1, min(dMax, d+2):end), [], 3);
%                     
%                     
%                     L1(j, i, d) = C(j, i, d) + min([L1(j, i-1, d), ...
%                                                     (L1(j, i-1, min(dMax, d+1)) + P1), ...
%                                                     (L1(j, i-1, max(1, d-1)) + P1), ...
%                                                     (min(minLeft, minRight) + P2)]) - min(L1(j, i-1, :), [], 3);
%                 end
%                 
%                 
%             end
%             
%         end
%     end
%     
%     %reverse path to calculate the L3/L4
%     for j = rows:-1:1
%         for i = cols:-1:1
%             %initialize the Path cost L3 at the bottom bounary
%             if(j == rows)
%                 L3(j, i, :) = C(j, i, :);
%             else
%                 for d = 1:dMax
%                     minLeft  = min(L3(j+1, i, 1:max(1, d-2)), [], 3);
%                     minRight = min(L3(j+1, i, min(dMax, d+2):end), [], 3);
%                     
%                     
%                     L3(j, i, d) = C(j, i, d) + min([L3(j+1, i, d), ...
%                                                     L3(j+1, i, min(dMax, d+1)) + P1, ...
%                                                     L3(j+1, i, max(1, d-1)) + P1, ...
%                                                     min(minLeft, minRight) + P2]) - min(L3(j+1, i, :), [], 3);
%                 end
%             end
%    
%             %initialize the Path cost L1 at the left bounary
%             if(i == cols)
%                 L4(j, 1, :) = C(j, 1, :);
%                 
%             else
%                 for d = 1:dMax
%                     minLeft  = min(L4(j, i+1, 1:max(1, d-2)), [], 3);
%                     minRight = min(L4(j, i+1, min(dMax, d+2):end), [], 3);
%                     
%                     
%                     L4(j, i, d) = C(j, i, d) + min([L4(j, i+1, d), ...
%                                                     L4(j, i+1, min(dMax, d+1)) + P1, ...
%                                                     L4(j, i+1, max(1, d-1)) + P1, ...
%                                                     min(minLeft, minRight) + P2]) - min(L4(j, i+1, :), [], 3);
%                 end
%                 
%                 
%             end
%             
%         end
%     end
    
    toc;
    L2 = permute(L2, [2, 1, 3]);
    L3 = permute(L3, [2, 1, 3]);
%     [~, bestD] = min(L1 + L2 + L3 + L4, [], 3);
    [~, bestD] = min(L1+L4+L2+L3, [], 3);
end
