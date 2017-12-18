% the hole filling algorithm from KITTI benchmark
function output = scanline_in_fill(input)
    height = size(input, 1);
    width = size(input, 2);
    
    for v=1:height
        count = 0;
        
        for u=1:width
            if(~isnan(input(v, u)))
                if(count >=1)
                    u1 = u-count;
                    u2 = u-1;
                    if(u1 > 1 && u2 <width)                
                        fu_ipol = min(input(v, u1-1, 1), input(v, u2 + 1, 1));
%                         fv_ipol = min(input(v, u1-1, 2), input(v, u2 + 1, 2));
                        for u_curr = u1:u2
                            input(v, u_curr, 1) = fu_ipol;
%                             input(v, u_curr, 2) = fv_ipol;
                        end
                    end
                end
                
                count = 0;
            else
                count = count + 1;
            end
        end
        %extrapolate to the left
        for u=1:width
            if ~isnan(input(v, u))
                for u2=1:u-1
                    input(v, u2, :) = input(v, u, :);
                end
                break;
            end
        end
        %extrapolate to the right
        for u=width:-1:1
            if ~isnan(input(v, u))
                for u2=u+1:width
                    input(v, u2, :) = input(v, u, :);
                end
                break;
            end
        end
    end
    
    
    for u= 1:width
        %extrapolate to the top
        for v=1:height
            if ~isnan(input(v, u))
                for v2=1:v-1
                    input(v2, u, :) = input(v, u, :);
                end
                break;
            end
        end
        %extrapolate to the bottom
        for v=height:-1:1
            if ~isnan(input(v, u))
                for v2=v+1:height
                    input(v2, u, :) = input(v, u, :);
                end
                break;
            end
        end
    end
    output = input;
end