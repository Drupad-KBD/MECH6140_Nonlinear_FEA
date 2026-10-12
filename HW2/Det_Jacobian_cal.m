function [J_det] = Det_Jacobian_cal(ShapeFun_selection,xi_cordi_4N, x_cap)
%% function to compute det of jacobian (to plot folding element)  

dN_cap = ShapeFun_selection(xi_cordi_4N);

J_dum = dN_cap * x_cap;

J = [J_dum(1,1), J_dum(3,1);J_dum(2,1),J_dum(4,1)];

J_det = det(J);

end