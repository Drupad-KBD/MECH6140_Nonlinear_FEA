classdef ElementBase
    % Shared public element interface, calculations delegated to ElementCalculations.
    properties (SetAccess=protected)
        Nodal_Locations
        thickness
        material
        nNodes
        gaussOrder
    end
    methods
        function obj=ElementBase(X,material,thickness,nNodes,gaussOrder)
            assert(isequal(size(X),[nNodes,2]),'Incorrect element coordinates');
            assert(isscalar(thickness)&&thickness>0,'Invalid thickness');
            obj.Nodal_Locations=X; obj.material=material; obj.thickness=thickness;
            obj.nNodes=nNodes; obj.gaussOrder=gaussOrder;
        end
        function [J,detJ,dNdX]=Jacobian_cal(obj,xi)
            [J,detJ,dNdX]=fea.elements.ElementCalculations.Jacobian_cal(obj,xi);
        end
        function H=Displacement_Gradient_cal(obj,xi,q)
            H=fea.elements.ElementCalculations.Displacement_Gradient_cal(obj,xi,q);
        end
        function B=B_Matrix(obj,xi)
            B=fea.elements.ElementCalculations.B_Matrix(obj,xi);
        end
        function B=Nonlinear_B_Matrix_cal(obj,xi,F)
            B=fea.elements.ElementCalculations.Nonlinear_B_Matrix_cal(obj,xi,F);
        end
        function [W,F,E,S,H,eps]=evaluatePoint(obj,xi,q)
            [W,F,E,S,H,eps]=fea.elements.ElementCalculations.evaluatePoint(obj,xi,q);
        end
        function q=Prescribed_Nodal_Displacements_cal(obj,u_fun,v_fun)
            q=zeros(2*obj.nNodes,1);
            for k=1:obj.nNodes
                x=obj.Nodal_Locations(k,1);y=obj.Nodal_Locations(k,2);
                q(2*k-1)=u_fun(x,y);q(2*k)=v_fun(x,y);
            end
        end
        function U=Total_Strain_Energy_cal(obj,q)
            U=fea.elements.ElementCalculations.Total_Strain_Energy_cal(obj,q);
        end
        function fint=Internal_Force_GP_cal(obj,q)
            fint=fea.elements.ElementCalculations.Internal_Force_GP_cal(obj,q);
        end
        function K=Gauss_Element_stiffness_cal(obj)
            K=fea.elements.ElementCalculations.Gauss_Element_stiffness_cal(obj);
        end
        function J=Deformation_Jacobian_cal(obj,xi,q)
            H=obj.Displacement_Gradient_cal(xi,q);J=det(eye(2)+H);
        end
        function sigma=Cauchy_Stress_cal(obj,xi,q)
            [~,F,~,S]=obj.evaluatePoint(xi,q); J=det(F);
            assert(J>0,'Invalid deformation');sigma=F*S*F'/J;
        end
        function out=calculate(obj,q)
            out.q=q(:);out.Fint=obj.Internal_Force_GP_cal(q);
            out.Energy=obj.Total_Strain_Energy_cal(q);
            out.GaussPoints=struct('xi',{},'F',{},'E',{},'S',{},'W',{},'epsilon',{});
            [gp,~]=fea.elements.ElementCalculations.gauss_location_weight(obj.gaussOrder);
            for i=1:numel(gp)
                for j=1:numel(gp)
                    xi=[gp(i),gp(j)];
                    [W,F,E,S,~,epsilon]=obj.evaluatePoint(xi,q);
                    out.GaussPoints(end+1)=struct('xi',xi,'F',F,'E',E,'S',S,'W',W,'epsilon',epsilon);
                end
            end
        end
        function [fint,KL,KNL,detail]=Element_Response_cal(obj,q)
            [fint,KL,KNL,detail]=fea.elements.ElementCalculations.Element_Response_cal(obj,obj.material,q);
        end
    end
end
