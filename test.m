I0 = imread('000004_10.png');
I1 = imread('000004_11.png');
P0 = read_calib_file('000004.txt');


[w, t] = mono_vo(I0, I1, [P0(1,1) P0(2,2)], [P0(1, 3) P0(2, 3)]);
w
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


