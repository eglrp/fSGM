function [ D2 ] = calc_disp_from_first( D1, Pd0, normDirect, O, vMax, n )
%CALC_DISP_FROM_FIRST Summary of this function goes here
%   Detailed explanation goes here
rows = size(D1, 1);
cols = size(D1, 2);
D2 = -1*ones(size(D1));

for j = 1:rows
    for i = 1:cols
        vzInd = D1(j, i);
        disp = vzInd2Disp(vzInd, O(j, i), vMax, n);
        
        off = reshape(normDirect(j, i, :), 2, 1);
        p2 = reshape(Pd0(j, i,:), 2, 1) + disp.*off;
        
        s0 = floor(p2);
        s1 = s0 + 1;
        
        sx0 = s0(1);
        sx1 = s1(1);
        sy0 = s0(2);
        sy1 = s1(2);
        
        if(sx0 >= 1&& sx0 <= cols && sy0 >= 1 && sy0 <= rows) 
            if(D2(sy0, sx0) == 0 || D2(sy0, sx0) < D1(j, i))
                D2(sy0, sx0) = D1(j, i);
            end
        end
        
        if(sx1 >= 1&& sx1 <= cols && sy0 >= 1 && sy0 <= rows) 
            if(D2(sy0, sx1) == 0 || D2(sy0, sx1) < D1(j, i))
                D2(sy0, sx1) = D1(j, i);
            end
        end
        
        if(sx0 >= 1&& sx0 <= cols && sy1 >= 1 && sy1 <= rows) 
            if(D2(sy1, sx0) == 0 || D2(sy1, sx0) < D1(j, i))
                D2(sy1, sx0) = D1(j, i);
            end
        end
        
        if(sx1 >= 1&& sx1 <= cols && sy1 >= 1 && sy1 <= rows) 
            if(D2(sy1, sx1) == 0 || D2(sy1, sx1) < D1(j, i))
                D2(sy1, sx1) = D1(j, i);
            end
        end
    end
end


end

