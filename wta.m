function flow = wta(C)
%winner takes all on the Cost volumn to get the flow vector
R = (sqrt(size(C, 3)) - 1)/2;

[M, I] = min(C, [],  3);


[y, x] = ind2sub([2*R + 1, 2*R + 1],I);

flow(:,:,1) = x-R-1;
flow(:,:,2) = y-R-1;


end