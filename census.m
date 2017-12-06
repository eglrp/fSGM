function cen = census(img, halfWin)
%Census transform of img
% halfWin the half window size of the census transform

% cen is the returned result, in MxNxD logical array, where D =
% (2*halfWin+1)*(2*halfWin + 1) - 1
% -1 means excluding the self compare result (always 0)
% logical "true" means the center is pixel larger than neighouring  value.

if(nargin < 2)
    halfWin = 2;
end

if(size(img, 3) > 1)
    img = rgb2gray(img);
end

rows = size(img, 1);
cols = size(img, 2);

bins = (2*halfWin + 1)*(2*halfWin+1);

cen = false(rows, cols,  bins- 1); % exclude the center comparision 

for j = 1+halfWin:rows-halfWin
    for i = 1+halfWin:cols-halfWin
        ind = 1;
        for wj=-halfWin:halfWin
            for wi = -halfWin:halfWin
                if((wj ~= 0 || wi ~= 0)) %exclude the center
                    if(img(j, i) > img(j+wj, i+wi))
                        cen(j, i, ind) = true;
                    end
                    ind = ind + 1;
                end
            end
        end
        assert(ind == bins);
    end
end
