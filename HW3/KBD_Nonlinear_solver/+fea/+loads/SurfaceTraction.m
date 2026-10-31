classdef SurfaceTraction
    % Uniform dead traction on a reference Q4/Q8 edge.
    % Local edges for Q4: [1 2], [2 3], [3 4], [4 1].
    % Local edges for Q8: [1 5 2], [2 6 3], [3 7 4], [4 8 1].
    methods (Static)
        function f=EdgeForce_cal(element,edge,traction)
            assert(numel(traction)==2,'Traction must be [tx ty]');
            n=element.nNodes;f=zeros(2*n,1);
            [g,w]=fea.elements.ElementCalculations.gauss_location_weight( ...
                2+(n==8));
            X=element.Nodal_Locations;
            for k=1:numel(g)
                a=g(k);
                switch edge
                    case 1,xi=[a,-1]; direction=[1;0];
                    case 2,xi=[1,a];direction=[0;1];
                    case 3,xi=[-a,1];direction=[-1;0];
                    case 4,xi=[-1,-a];direction=[0;-1];
                    otherwise,error('Edge must be 1,2,3,4');
                end
                [N,dN]=element.ShapeFunctions(xi);
                t=dN'*direction; dL=norm(X'*t);
                for i=1:n
                    f(2*i-1:2*i)=f(2*i-1:2*i)+ ...
                        N(i)*traction(:)*element.thickness*dL*w(k);
                end
            end
        end
    end
end
