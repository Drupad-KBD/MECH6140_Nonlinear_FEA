function [S] = Second_Piola_Kirchhoff_Stress_cal(E,lambda_L,mu)

%% Second Piola-Kirchhoff stress

S = lambda_L*trace(E)*eye(2) + 2*mu*E;

end