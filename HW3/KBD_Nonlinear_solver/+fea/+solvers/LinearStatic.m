
classdef LinearStatic
    %LINEARSTATIC Linear equilibrium with elimination or penalty BCs.
    %
    % Original call remains supported:
    %
    % solve(K,loadHistory,fixedDOF)
    %
    % Default boundary condition method: elimination.

    methods (Static)

        function result = solve(K,loadHistory,fixedDOF,options)

            if nargin < 4 || isempty(options)
                options = struct();
            end

            %% Default settings

            if ~isfield(options,'bcMethod')
                options.bcMethod = 'elimination';
            end

            if ~isfield(options,'penaltyFactor')
                options.penaltyFactor = 1e6;
            end

            if ~isfield(options,'prescribedValues')
                options.prescribedValues = 0;
            end

            %% Check dimensions

            nDOF = size(K,1);

            assert(ismatrix(K) && ...
                size(K,2)==nDOF && ...
                size(loadHistory,1)==nDOF, ...
                'Incompatible stiffness/load sizes.');

            fixedDOF = unique(fixedDOF(:)');

            assert(all(fixedDOF>=1 & ...
                fixedDOF<=nDOF & ...
                fixedDOF==fix(fixedDOF)), ...
                'Invalid fixed DOF indices.');

            freeDOF = setdiff(1:nDOF,fixedDOF);

            method = lower(char(options.bcMethod));

            result.penalty = [];

            %% Boundary condition enforcement

            switch method

                case 'elimination'

                    assert(all(options.prescribedValues==0), ...
                        ['Nonzero prescribedValues are supported ' ...
                         'with bcMethod=penalty.']);

                    qHistory = zeros( ...
                        nDOF,size(loadHistory,2));

                    qHistory(freeDOF,:) = ...
                        K(freeDOF,freeDOF) \ ...
                        loadHistory(freeDOF,:);

                case 'penalty'

                    penalty = ...
                        fea.assembly.PenaltyMethod.create( ...
                        K, ...
                        fixedDOF, ...
                        options.prescribedValues, ...
                        options.penaltyFactor);

                    [Kmod,Fmod] = ...
                        fea.assembly.PenaltyMethod.applyLinear( ...
                        K,loadHistory,penalty);

                    qHistory = Kmod \ Fmod;

                    result.penalty = penalty;

                otherwise

                    error( ...
                        'Unknown boundary condition method: %s', ...
                        method);
            end

            %% Return results

            result.q = qHistory(:,end);
            result.qHistory = qHistory;

            result.fixedDOF = fixedDOF;
            result.freeDOF = freeDOF;

            result.bcMethod = method;

        end
    end
end
