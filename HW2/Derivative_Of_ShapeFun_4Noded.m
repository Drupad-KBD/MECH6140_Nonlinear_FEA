
function [dN_cap] = Derivative_Of_ShapeFun_4Noded(xi_cordi_4N)

%% function to compute Ncap matrix using derivatives of shape functions %%

% accessing the xi and eta from the input matrix (xi,eta) % 
xi  = xi_cordi_4N(1,1); % accessing xi (xi,eta) so (1,1)
eta = xi_cordi_4N(1,2); % accessing eta (xi,eta) so (1,2)
 
% eqns to calculate shape functions %
% xi derivatives
dN1dxi = eta/4 - 1/4;
dN2dxi = 1/4 - eta/4;
dN3dxi = eta/4 + 1/4;
dN4dxi = -eta/4 - 1/4;

% deta derivatives
dN1deta = xi/4 - 1/4;
dN2deta = -xi/4-1/4;
dN3deta = xi/4 + 1/4;
dN4deta = 1/4 - xi/4;

% storing every thing in 1 matrix %
dN_cap = [dN1dxi,0,dN2dxi,0,dN3dxi,0,dN4dxi,0;... 
         dN1deta,0,dN2deta,0,dN3deta,0,dN4deta,0;...
         0,dN1dxi,0,dN2dxi,0,dN3dxi,0,dN4dxi;...
         0,dN1deta,0,dN2deta,0,dN3deta,0,dN4deta];
end