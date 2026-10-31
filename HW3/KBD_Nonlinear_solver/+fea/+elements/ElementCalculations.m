classdef ElementCalculations
    % Shared 2D plane-strain element kinematics, quadrature and TL response.
    % Voigt convention: [xx yy engineering_xy] for strain, [Sxx Syy Sxy] for stress.
    methods (Static)
        function [gp,w]=gauss_location_weight(n)
            switch n
                case 1, gp=0; w=2;
                case 2, gp=[-1,1]/sqrt(3);w=[1,1];
                case 3, gp=[-sqrt(3/5),0,sqrt(3/5)];w=[5/9,8/9,5/9];
                otherwise,error('Unsupported quadrature order');
            end
        end
        function F=Deformation_Gradient_cal(H)
            F=eye(2)+H;
        end
        function [E,E_vec]=Green_Lagrange_Strain_cal(H)
            E=.5*(H+H'+H'*H);
            E_vec=[E(1,1);E(2,2);2*E(1,2)];
        end
        function [e,e_vec]=Infinitesimal_Strain_cal(H)
            e=.5*(H+H');e_vec=[e(1,1);e(2,2);2*e(1,2)];
        end
        function [J,detJ,dNdX]=Jacobian_cal(element,xi)
            [~,dN]=element.ShapeFunctions(xi);
            J=dN*element.Nodal_Locations;detJ=det(J);
            assert(detJ>0,'Inverted or degenerate element');
            dNdX=J\dN;
        end
        function H=Displacement_Gradient_cal(element,xi,q)
            assert(numel(q)==2*element.nNodes,'Wrong element DOFs');
            [~,~,d]=fea.elements.ElementCalculations.Jacobian_cal(element,xi);
            H=reshape(q(:),2,element.nNodes)*d';
        end
        function B=B_Matrix(element,xi)
            [~,~,d]=fea.elements.ElementCalculations.Jacobian_cal(element,xi);
            B=zeros(3,2*element.nNodes);
            for a=1:element.nNodes
                k=2*a-1;B(1,k)=d(1,a);B(2,k+1)=d(2,a);
                B(3,k)=d(2,a);B(3,k+1)=d(1,a);
            end
        end
        function B=Nonlinear_B_Matrix_cal(element,xi,F)
            [~,~,d]=fea.elements.ElementCalculations.Jacobian_cal(element,xi);
            B=zeros(3,2*element.nNodes);
            for a=1:element.nNodes
                dx=d(1,a);dy=d(2,a);k=2*a-1;
                B(:,k)=[F(1,1)*dx;F(1,2)*dy;F(1,1)*dy+F(1,2)*dx];
                B(:,k+1)=[F(2,1)*dx;F(2,2)*dy;F(2,1)*dy+F(2,2)*dx];
            end
        end
        function [AL,BL,BNL,F,detJ]=Strain_Operators_cal(element,xi,q)
            [~,detJ,d]=fea.elements.ElementCalculations.Jacobian_cal(element,xi);
            H=fea.elements.ElementCalculations.Displacement_Gradient_cal(element,xi,q);
            F=eye(2)+H;BNL=zeros(4,2*element.nNodes);
            for a=1:element.nNodes
                k=2*a-1;BNL(1:2,k)=d(:,a);BNL(3:4,k+1)=d(:,a);
            end
            AL=[F(1,1),0,F(2,1),0;0,F(1,2),0,F(2,2); ...
                F(1,2),F(1,1),F(2,2),F(2,1)];
            BL=AL*BNL;
        end
        function [W,F,E,S,H,eps]=evaluatePoint(element,xi,q)
            H=fea.elements.ElementCalculations.Displacement_Gradient_cal(element,xi,q);F=eye(2)+H;
            [E,~]=fea.elements.ElementCalculations.Green_Lagrange_Strain_cal(H);
            [eps,~]=fea.elements.ElementCalculations.Infinitesimal_Strain_cal(H);
            S=element.material.Second_Piola_Kirchhoff_Stress_cal(E);
            if ismethod(element.material,'Strain_Energy_Density_cal')
                W=element.material.Strain_Energy_Density_cal(E);
            else
                W=NaN; % No energy potential provided by constitutive model.
            end
        end
        function U=Total_Strain_Energy_cal(element,q)
            [gp,w]=fea.elements.ElementCalculations.gauss_location_weight(element.gaussOrder);U=0;
            for i=1:numel(gp)
                for j=1:numel(gp)
                    xi=[gp(i),gp(j)];[W,~,~,~,~,~]=fea.elements.ElementCalculations.evaluatePoint(element,xi,q);
                    [~,detJ]=fea.elements.ElementCalculations.Jacobian_cal(element,xi);
                    U=U+W*detJ*w(i)*w(j)*element.thickness;
                end
            end
        end
        function fint=Internal_Force_GP_cal(element,q)
            [fint,~,~]=fea.elements.ElementCalculations.Element_Response_cal(element,element.material,q);
        end
        function K=Gauss_Element_stiffness_cal(element)
            D=element.material.D_matrix_PlaneStrain();ndof=2*element.nNodes;
            K=zeros(ndof);[gp,w]=fea.elements.ElementCalculations.gauss_location_weight(element.gaussOrder);
            for i=1:numel(gp)
                for j=1:numel(gp)
                    xi=[gp(i),gp(j)];B=fea.elements.ElementCalculations.B_Matrix(element,xi);
                    [~,detJ]=fea.elements.ElementCalculations.Jacobian_cal(element,xi);
                    K=K+B'*D*B*detJ*w(i)*w(j)*element.thickness;
                end
            end
        end
        function [fint,KL,KNL,detail]=Element_Response_cal(element,material,q)
            C=fea.elements.ElementCalculations;ndof=2*element.nNodes;
            fint=zeros(ndof,1);KL=zeros(ndof);KNL=zeros(ndof);
            [gp,w]=fea.elements.ElementCalculations.gauss_location_weight(element.gaussOrder);
            detail=struct('xi',{},'F',{},'E',{},'S',{},'BL',{},'BNL',{},'dV0',{});
            for i=1:numel(gp)
                for j=1:numel(gp)
                    xi=[gp(i),gp(j)];
                    [~,BL,BNL,F,detJ]=fea.elements.ElementCalculations.Strain_Operators_cal(element,xi,q);
                    [E,~]=fea.elements.ElementCalculations.Green_Lagrange_Strain_cal(F-eye(2));
                    if ismethod(material,'Material_Tangent_cal')
                        [S,D]=material.Material_Tangent_cal(E);
                    else
                        S=material.Second_Piola_Kirchhoff_Stress_cal(E);
                        D=material.D_matrix_PlaneStrain();
                    end
                    sv=[S(1,1);S(2,2);S(1,2)];Stilde=kron(eye(2),S);
                    dV0=detJ*w(i)*w(j)*element.thickness;
                    fint=fint+BL'*sv*dV0;
                    KL=KL+BL'*D*BL*dV0;
                    KNL=KNL+BNL'*Stilde*BNL*dV0;
                    detail(end+1)=struct('xi',xi,'F',F,'E',E,'S',S, ...
                        'BL',BL,'BNL',BNL,'dV0',dV0);
                end
            end
        end
    end
end
