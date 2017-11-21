function C = calc_cost(I0, I1, R)

I0 = double(I0);
I1 = double(I1);

% calculate Cost Volumn for I0/I1
if nargin < 3
    R = 10; % the half search range 
end

rows = size(I0, 1);
cols = size(I0, 2);
C = nan(rows, cols, (R*2 + 1)*(R*2+1));

for j = -R:R
    for i = -R:R
        It = imtranslate(I1,[-i, -j]);
        diff = abs(It-I0);
        C(:,:, (j+R)*(2*R+1) + i+R+1) = sum(diff, 3);
    end 
end

end