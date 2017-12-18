function [ D1 ] = forward_backward_check( D1, D2, Pd0, normDirect, O, vMax, n )
%FORWARD_BACKWARD_CHECK Summary of this function goes here
%   Detailed explanation goes here
rows = size(D1, 1);
cols = size(D1, 2);
thr = 2.0;

for j = 1:rows
    for i = 1:cols
        vzInd = D1(j, i);
        
        if(isnan(vzInd))
            continue;
        end
        disp = vzInd2Disp(vzInd, O(j, i), vMax, n);
        
        off = reshape(normDirect(j, i, :), 2, 1);
        p2 = reshape(Pd0(j, i,:), 2, 1) + disp.*off;
        
        p2 = round(p2);
        
        if(p2(1) < 1 || p2(1) > cols || p2(2)<1 || p2(2)>rows)
            D1(j, i) = nan;
            continue;
        end
        
        if(D2(p2(2), p2(1)) == -1)
            D1(j, i) = nan;
            continue;
        end
        
        if(abs(D1(j, i) - D2(p2(2), p2(1))) > thr)
            D1(j, i) = nan;
        end

    end

end

