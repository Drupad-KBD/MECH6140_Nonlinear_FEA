function [u] = Displacement_cal(ShapeFun_8Noded,Node_location,q)
% calling and calculating the shape function %
N = ShapeFun_8Noded(Node_location);

% calculating the displacments
u = N * q;

end