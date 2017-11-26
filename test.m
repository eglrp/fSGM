I0 = imread('000004_10.png');
I1 = imread('000004_11.png');
P0 = read_calib_file('000004.txt');


[w, t] = mono_vo(I0, I1, [P0(1,1) P0(2,2)], [P0(1, 3) P0(2, 3)]);
flow = rotation_motion(w, P0(1,1), [P0(1, 3) P0(2, 3)], [size(I0, 2), size(I0, 1)]);

Fc = flow_to_color(flow);
imshow(Fc);

flow_field(flow, I0);
% C = calc_cost(I0, I1, 20);
% 
% flow = wta(C);
% 
% fColor = flowToColor(flow, 17.6);
% 
% 
% 
% gtFlow = flow_read('flow10.flo');
% gtColor = flowToColor(gtFlow);
% imshow([fColor gtColor]);
% 


