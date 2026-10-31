%% Q4 and Q8 element verification - same homogeneous deformation
clearvars;clc;
E=70000;nu=0.3;t=1;
mat=fea.materials.LinearElastic(E,nu);
X4=[0 0;1 0;1 1;0 1];
X8=[X4;0.5 0;1 0.5;0.5 1;0 0.5];
Q4=fea.elements.Quad4(X4,mat,t);
Q8=fea.elements.Quad8(X8,mat,t);
assert(norm(Q4.ShapeFunctions([0 0]))>0);
assert(norm(Q8.ShapeFunctions([0 0]))>0);
% Affine uniform stretch u_x=0.1 X, u_y=-0.02 Y
for el={Q4,Q8}
    elem=el{1};X=elem.Nodal_Locations;
    q=reshape([0.1*X(:,1)';-0.02*X(:,2)'],[],1);
    [f,KL,KNL,details]=elem.Element_Response_cal(q);
    assert(all(isfinite(f))&&all(isfinite(KL(:)))&&all(isfinite(KNL(:))));
    for j=1:numel(details)
        assert(norm(details(j).F-diag([1.1 0.98]),'fro')<1e-11);
    end
    fprintf('%s: %d DOFs, %d Gauss points, net internal force = %.3e\n', ...
        class(elem),numel(q),numel(details),norm([sum(f(1:2:end));sum(f(2:2:end))]));
    assert(norm([sum(f(1:2:end));sum(f(2:2:end))])<1e-8);
end
fprintf('PASS: Q4 and Q8 homogeneous deformation tests.\n');
