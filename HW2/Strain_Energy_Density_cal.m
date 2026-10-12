function [W] = Strain_Energy_Density_cal(E,lambda_L,mu)

%% St. Venant-Kirchhoff strain energy density

trace_E = trace(E);

E_double_E = sum(sum(E.*E));

W = (lambda_L/2)*(trace_E^2) + mu*E_double_E;

end