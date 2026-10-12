function [D_pl_strain] = D_matrix_PlaneStrain(E,v)
% plane Strain condition 
% Calculate D matrix
D_pl_strain=(E/((1+v)*(1-2*v)))*[(1-v),v,0; v,(1-v),0; 0,0,(0.5*(1-2*v))];

end