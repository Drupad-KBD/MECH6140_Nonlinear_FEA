function [B] = Nonlinear_B_Matrix_cal( ...
    ShapeFun_derivative,xi_cordi,x_cap,F)

%% ================================================================
% Nonlinear Total Lagrangian B Matrix
% ================================================================

%% Shape function derivative

dN_cap = ...
    ShapeFun_derivative(xi_cordi);


%% Reference Jacobian inverse

J_inv = ...
    Jacobian_inv_cal( ...
    ShapeFun_derivative, ...
    xi_cordi, ...
    x_cap);


%% Finite deformation A matrix

A_finite = [ ...
    F(1,1)   0        F(2,1)   0;
    0        F(1,2)   0        F(2,2);
    F(1,2)   F(1,1)   F(2,2)   F(2,1)];


%% Nonlinear B matrix

B = ...
    A_finite ...
    * J_inv ...
    * dN_cap;

end