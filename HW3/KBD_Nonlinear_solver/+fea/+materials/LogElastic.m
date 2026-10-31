
classdef LogElastic
    %LOGELASTIC Abaqus-like finite-strain elasticity, plane strain.
    %
    % Uses logarithmic strain and an isotropic Cauchy stress law.
    % Returns second Piola-Kirchhoff stress and dS/dE for TL FEA.
    %
    % This is an Abaqus-like comparison model, not an exact
    % reproduction of every Abaqus material integration path.

    properties (SetAccess = private)
        E
        nu
        lambda_L
        mu
    end

    methods
        function obj = LogElastic(E, nu)

            validateattributes(E, {'numeric'}, ...
                {'scalar','real','positive','finite'});

            assert(isscalar(nu) && isfinite(nu) && ...
                nu > -1 && nu < 0.5, ...
                'Invalid Poisson ratio.');

            obj.E = E;
            obj.nu = nu;

            obj.lambda_L = E*nu/((1+nu)*(1-2*nu));
            obj.mu = E/(2*(1+nu));
        end

        function S = Second_Piola_Kirchhoff_Stress_cal(obj, EGL)
            % Convert Green-Lagrange strain into logarithmic strain,
            % calculate spatial stress, then return PK2 stress.

            EGL = 0.5*(EGL + EGL');
            Cmat = eye(2) + 2*EGL;

            [V, D] = eig(Cmat);
            principalC = diag(D);

            if any(principalC <= 0) || any(~isfinite(principalC))
                error('Invalid deformation: C must be positive definite.');
            end

            % Principal logarithmic strains
            eLog = 0.5*log(principalC);

            % Plane strain: out-of-plane stretch = 1
            J = sqrt(prod(principalC));

            % Principal Cauchy stresses
            sigmaPrincipal = obj.lambda_L*sum(eLog) ...
                + 2*obj.mu*eLog;

            % Second Piola-Kirchhoff stress
            SPrincipal = J*(sigmaPrincipal./principalC);

            S = V*diag(SPrincipal)*V';
            S = 0.5*(S + S');
        end

        function [S, C] = Material_Tangent_cal(obj, EGL)
            % Numerical consistent constitutive tangent:
            %
            % C = d[S11 S22 S12] / d[E11 E22 2*E12]

            S = obj.Second_Piola_Kirchhoff_Stress_cal(EGL);
            C = zeros(3,3);

            basis = zeros(2,2,3);
            basis(:,:,1) = [1 0; 0 0];
            basis(:,:,2) = [0 0; 0 1];
            basis(:,:,3) = [0 0.5; 0.5 0];

            h = 1e-6*max(1,norm(EGL,'fro'));

            for j = 1:3

                Ep = EGL + h*basis(:,:,j);
                Em = EGL - h*basis(:,:,j);

                Sp = obj.Second_Piola_Kirchhoff_Stress_cal(Ep);
                Sm = obj.Second_Piola_Kirchhoff_Stress_cal(Em);

                SpVec = [Sp(1,1); Sp(2,2); Sp(1,2)];
                SmVec = [Sm(1,1); Sm(2,2); Sm(1,2)];

                C(:,j) = (SpVec-SmVec)/(2*h);
            end
        end

        function D = D_matrix_PlaneStrain(obj)
            % Small-strain elastic stiffness for comparison.

            l = obj.lambda_L;
            m = obj.mu;

            D = [l+2*m, l, 0;
                 l, l+2*m, 0;
                 0, 0, m];
        end
    end
end
