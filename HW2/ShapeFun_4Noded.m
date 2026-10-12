function [N] = ShapeFun_4Noded(xi_cordi_4N)
%% shape function 4 node calculator %%

% accessing the xi and eta from the input matrix (xi,eta) % 
xi  = xi_cordi_4N(1,1); % accessing xi (xi,eta) so (1,1)
eta = xi_cordi_4N(1,2); % accessing eta (xi,eta) so (1,2)

% eqns to calculate shape functions %

N1 = (xi-1)*(eta-1)*0.25;
N2 = (xi+1)*(eta-1)*(-0.25);
N3 = (xi+1)*(eta+1)*0.25;
N4 = (xi-1)*(eta+1)*(-0.25);

% storing every thing in 1 matrix %
N = [N1,0,N2,0,N3,0,N4,0; 0,N1,0,N2,0,N3,0,N4];

end