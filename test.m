seqIndex = 4;

I0 = imread(sprintf('%06d_10.png', seqIndex));
I1 = imread(sprintf('%06d_11.png', seqIndex));

% F = estimate_Fmatrix(I0, I1);

P0 = read_calib_file(sprintf('%06d.txt', seqIndex));


[w, t, E1] = mono_vo(I0, I1, [P0(1,1) P0(2,2)], [P0(1, 3) P0(2, 3)]);

flow = rotation_motion(w, P0(1,1), [P0(1, 3) P0(2, 3)], [size(I0, 2), size(I0, 1)]);

% Fc = flow_to_color(flow);
% imshow(Fc);

flow_field(flow, I0);



