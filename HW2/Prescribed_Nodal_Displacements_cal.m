function [q] = Prescribed_Nodal_Displacements_cal(nodalLocations,u_fun,v_fun)
nNode = size(nodalLocations,1);
q = zeros(2*nNode,1);

for i = 1:nNode
    X = nodalLocations(i,1);
    Y = nodalLocations(i,2);
    q(2*i-1) = u_fun(X,Y);
    q(2*i) = v_fun(X,Y);
end

end