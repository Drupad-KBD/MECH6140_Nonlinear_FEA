function [epsilon,epsilon_vec] = Infinitesimal_Strain_cal(H)
epsilon = 0.5 * (H + H.');
epsilon_vec = [epsilon(1,1); epsilon(2,2); 2*epsilon(1,2)];
end