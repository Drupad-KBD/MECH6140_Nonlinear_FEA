function [sigma] = Cauchy_Stress_cal(F,S,J_deformation)

%% Cauchy stress

sigma = (1/J_deformation)*F*S*F.';

end