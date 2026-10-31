classdef LinearElastic
    %LINEARELASTIC St. Venant-Kirchhoff plane-strain material.
    % The out-of-plane Green-Lagrange strain is zero (F33=1).
    properties (SetAccess = private)
        E; nu; lambda_L; mu;
    end
    methods
        function obj = LinearElastic(E,nu)
            % Constructor: Young's modulus E and Poisson's ratio nu.
            assert(isscalar(E) && E>0 && nu>-1 && nu<0.5,'Invalid material');
            obj.E=E; obj.nu=nu;
            obj.lambda_L=E*nu/((1+nu)*(1-2*nu));
            obj.mu=E/(2*(1+nu));
        end
        function S = Second_Piola_Kirchhoff_Stress_cal(obj,EGL)
            % Return 2x2 second Piola-Kirchhoff stress, in reference frame.
            S=obj.lambda_L*trace(EGL)*eye(2)+2*obj.mu*EGL;
        end
        function W = Strain_Energy_Density_cal(obj,EGL)
            % Return linear elastic strain-energy density per reference volume.
            W=0.5*obj.lambda_L*trace(EGL)^2+obj.mu*sum(EGL(:).^2);
        end
        function D=D_matrix_PlaneStrain(obj)
            % Return small-strain plane-strain D in [xx,yy,engineering xy].
            l=obj.lambda_L; m=obj.mu;
            D=[l+2*m,l,0;l,l+2*m,0;0,0,m];
        end
    end
end
