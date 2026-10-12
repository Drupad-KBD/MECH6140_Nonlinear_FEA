function [x] = Coordinate_cal(ShapeFun_selection,node_location,x_cap)
% calling the shape function %
N = ShapeFun_selection(node_location);

[m,n]=size(x_cap);
if n == 2
    converted_X_cap = X_2_to_X_1converter(x_cap);

% calculating the xy coordinates
x = N * converted_X_cap;

elseif n== 1
 x = N * x_cap;  
end

end