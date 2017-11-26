function [P0, P1, P2, P3] = read_calib_file(filename)

fp = fopen(filename);

for i = 1:4
    fscanf(fp, "%s:");
    P{i} = fscanf(fp, "%f %f %f %f %f %f %f %f %f %f %f %f\n");
end
P0 = reshape(P{1}', [4, 3]);
P1 = reshape(P{2}', [4, 3]);
P2 = reshape(P{3}', [4, 3]);
P3 = reshape(P{4}', [4, 3]);

P0 = P0';
P1 = P1';
P2 = P2';
P3 = P3';

end