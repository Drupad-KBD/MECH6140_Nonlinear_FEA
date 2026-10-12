function [H] = Displacement_Gradient_cal(ShapeFun_derivative,xi_cordi,x_cap,q)
dN_cap = ShapeFun_derivative(xi_cordi);
J_inv = Jacobian_inv_cal(ShapeFun_derivative,xi_cordi,x_cap);

G = J_inv * dN_cap;
h = G * q;

H = [h(1),h(2);h(3),h(4)];

end