function [F] = Deformation_Gradient_cal(H)
F = eye(2) + H;
end