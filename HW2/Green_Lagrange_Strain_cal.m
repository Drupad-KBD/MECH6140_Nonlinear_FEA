function [E,E_vec] = Green_Lagrange_Strain_cal(H)
E = 0.5 * (H + H.' + H.' * H);

F_check = eye(2) + H;
E_check = 0.5 * (F_check.' * F_check - eye(2));

if norm(E - E_check,'fro') > 100 * eps(max(1,norm(H,'fro')^2))
    error('Green-Lagrange strain calculations do not agree.');
end

E_vec = [E(1,1); E(2,2); 2*E(1,2)];
end