function [J] = Jacobian_cal(ShapeFun_selection,xi_cordi_4N, x_cap)
%% function to compute jacobian 
[m,n]=size(x_cap);
% if n == 2
%     converted_X_cap = X_2_to_X_1converter(x_cap);
%     
% dN_cap = ShapeFun_selection(xi_cordi_4N);
% 
% J_dum = dN_cap * converted_X_cap;
% 
% J = [J_dum(1,1), J_dum(3,1);J_dum(2,1),J_dum(4,1)];
% 
% elseif n==1
     
dN_cap = ShapeFun_selection(xi_cordi_4N);

J_dum = dN_cap * x_cap;

J = [J_dum(1,1), J_dum(3,1);J_dum(2,1),J_dum(4,1)];
        

end