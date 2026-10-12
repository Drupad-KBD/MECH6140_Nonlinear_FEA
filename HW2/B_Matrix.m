function [B] = B_Matrix(ShapeFun_selection,xi_cordi_4N, x_cap)

%% function to calculate B_matrix %%
A = [1 0 0 0;0 0 0 1;0 1 1 0];

% calling function to compute N_cap matrix
N_cap = ShapeFun_selection(xi_cordi_4N);

% calling function to compute J inverse matrix 
J_inv = Jacobian_inv_cal(ShapeFun_selection,xi_cordi_4N, x_cap);

% Caluculating B matrix
B = A * J_inv * N_cap;

end