% classdef PenaltyMethod
%     %PENALTYMETHOD Linear displacement-constraint penalty, shared by solvers.
%     % beta = gamma*max(abs(diag(Kreference))) is fixed for the analysis.
%     % The constraint stiffness is Kp and the penalty internal force is
%     % fp(q)=Kp*(q-uBar). No stress calculation occurs in this class.
%     methods (Static)
%         function penalty = create(Kreference,fixedDOF,prescribedValues,gamma)
%             if nargin<4 || isempty(gamma), gamma=1e6; end
%             if nargin<3 || isempty(prescribedValues), prescribedValues=0; end
%             n=size(Kreference,1);
%             assert(size(Kreference,2)==n && n>0,'Reference stiffness must be square.');
%             assert(isscalar(gamma) && isfinite(gamma) && gamma>0, ...
%                 'Penalty multiplier must be positive and finite.');
%             fixedDOF=fixedDOF(:);
%             assert(numel(unique(fixedDOF))==numel(fixedDOF) && ...
%                 all(fixedDOF>=1 & fixedDOF<=n & fixedDOF==fix(fixedDOF)), ...
%                 'Invalid or repeated constrained DOF.');
%             if isscalar(prescribedValues)
%                 prescribedValues=repmat(prescribedValues,numel(fixedDOF),1);
%             end
%             prescribedValues=prescribedValues(:);
%             assert(numel(prescribedValues)==numel(fixedDOF) && ...
%                 all(isfinite(prescribedValues)), 'Incorrect prescribed displacements.');
%             scale=max(abs(diag(Kreference)));
%             assert(isfinite(scale) && scale>0, ...
%                 'Reference stiffness must have a positive diagonal scale.');
%             alpha=gamma*scale;
%             penalty.Kp=sparse(fixedDOF,fixedDOF,alpha,n,n);
%             penalty.uBar=zeros(n,1);
%             penalty.uBar(fixedDOF)=prescribedValues;
%             penalty.alpha=alpha;
%             penalty.gamma=gamma;
%             penalty.fixedDOF=fixedDOF;
%         end
%         function [Kmodified,Fmodified]=applyLinear(K,F,penalty)
%             Kmodified=K+penalty.Kp;
%             Fmodified=F+penalty.Kp*penalty.uBar;
%         end
%         function [fModified,Kmodified]=applyNonlinear(fint,Kt,q,penalty)
%             fModified=fint+penalty.Kp*(q-penalty.uBar);
%             Kmodified=Kt+penalty.Kp;
%         end
%     end
% end

classdef PenaltyMethod
    %PENALTYMETHOD Linear displacement-constraint penalty, shared by solvers.
    %
    % beta = gamma*max(abs(diag(Kreference))) is fixed for the analysis.
    %
    % The constraint stiffness is Kp and the penalty internal force is:
    %
    % fp(q) = Kp*(q-uBar)
    %
    % No stress calculation occurs in this class.

    methods (Static)

        function penalty = create(Kreference,fixedDOF,prescribedValues,gamma)

            if nargin < 4 || isempty(gamma)
                gamma = 1e6;
            end

            if nargin < 3 || isempty(prescribedValues)
                prescribedValues = 0;
            end

            n = size(Kreference,1);

            assert(size(Kreference,2)==n && n>0, ...
                'Reference stiffness must be square.');

            assert(isscalar(gamma) && isfinite(gamma) && gamma>0, ...
                'Penalty multiplier must be positive and finite.');

            fixedDOF = fixedDOF(:);

            assert(numel(unique(fixedDOF))==numel(fixedDOF) && ...
                all(fixedDOF>=1 & fixedDOF<=n & ...
                fixedDOF==fix(fixedDOF)), ...
                'Invalid or repeated constrained DOF.');

            if isscalar(prescribedValues)
                prescribedValues = repmat( ...
                    prescribedValues,numel(fixedDOF),1);
            end

            prescribedValues = prescribedValues(:);

            assert(numel(prescribedValues)==numel(fixedDOF) && ...
                all(isfinite(prescribedValues)), ...
                'Incorrect prescribed displacements.');

            %% Penalty stiffness

            scale = max(abs(diag(Kreference)));

            assert(isfinite(scale) && scale>0, ...
                'Reference stiffness must have a positive diagonal scale.');

            alpha = gamma*scale;

            penalty.Kp = sparse(fixedDOF,fixedDOF,alpha,n,n);

            %% Prescribed displacement vector

            penalty.uBar = zeros(n,1);
            penalty.uBar(fixedDOF) = prescribedValues;

            %% Store parameters

            penalty.alpha = alpha;
            penalty.gamma = gamma;
            penalty.fixedDOF = fixedDOF;
        end

        %% Linear equilibrium

        function [Kmodified,Fmodified] = applyLinear(K,F,penalty)

            Kmodified = K + penalty.Kp;

            Fmodified = F + penalty.Kp*penalty.uBar;
        end

        %% Nonlinear equilibrium

        function [fModified,Kmodified] = ...
                applyNonlinear(fint,Kt,q,penalty)

            fModified = fint + ...
                penalty.Kp*(q-penalty.uBar);

            Kmodified = Kt + penalty.Kp;
        end

    end
end
