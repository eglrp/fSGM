I0 = imread('frame10.png');
I1 = imread('frame11.png');

C = calc_cost(I0, I1, 20);

flow = wta(C);

fColor = flowToColor(flow, 17.6);



gtFlow = flow_read('flow10.flo');
gtColor = flowToColor(gtFlow);
imshow([fColor gtColor]);



