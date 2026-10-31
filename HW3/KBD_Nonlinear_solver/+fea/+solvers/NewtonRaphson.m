
classdef NewtonRaphson
    %NEWTONRAPHSON Incremental nonlinear equilibrium solver.
    %
    % response(q) returns [fint,Kt].
    %
    % Methods:
    %   full      - exact tangent every iteration
    %   modified  - tangent fixed within each load increment
    %   secant    - Broyden rank-one tangent update
    %   broyden   - alias for secant
    %   numerical - centered finite-difference tangent
    %
    % Boundary conditions:
    %   elimination (default)
    %   penalty
    %
    % Optional load-step bisection and iteration output.

    methods (Static)

        function result = solve( ...
                response,loadHistory,fixedDOF,options)

            if nargin < 4 || isempty(options)
                options = struct();
            end

            %% Default settings

            defaults = struct( ...
                'tolResidual',1e-10, ...
                'maxIter',25, ...
                'verbose',false, ...
                'printEvery',1, ...
                'method','full', ...
                'cutback',false, ...
                'maxCutbacks',8, ...
                'fdStep',1e-6, ...
                'bcMethod','elimination', ...
                'penaltyFactor',1e6, ...
                'prescribedValues',0, ...
                'checkpointFile','', ...
                'restartFrom','', ...
                'saveEvery',1, ...
                'restartKey','');

            names = fieldnames(defaults);

            for i = 1:numel(names)

                if ~isfield(options,names{i})
                    options.(names{i}) = defaults.(names{i});
                end
            end

            %% Solver method

            method = lower(char(options.method));

            assert(ismember(method, ...
                {'full','modified','secant','broyden','numerical'}), ...
                'Unknown nonlinear solver method.');

            assert(isscalar(options.verbose) && ...
                (islogical(options.verbose) || ...
                isnumeric(options.verbose)), ...
                'options.verbose must be true or false.');

            assert(isscalar(options.printEvery) && ...
                options.printEvery>=1 && ...
                options.printEvery==fix(options.printEvery), ...
                'printEvery must be a positive integer.');

            assert(isscalar(options.maxIter) && ...
                options.maxIter>=1 && ...
                options.maxIter==fix(options.maxIter), ...
                'Invalid maxIter.');

            assert(isscalar(options.maxCutbacks) && ...
                options.maxCutbacks>=0 && ...
                options.maxCutbacks==fix(options.maxCutbacks), ...
                'Invalid maxCutbacks.');

            assert(options.fdStep>0 && ...
                isfinite(options.fdStep), ...
                'fdStep must be positive.');

            assert(isscalar(options.saveEvery) && ...
                options.saveEvery>=1 && options.saveEvery==fix(options.saveEvery), ...
                'saveEvery must be a positive integer.');

            %% Degrees of freedom

            nDOF = size(loadHistory,1);
            nSteps = size(loadHistory,2)-1;

            assert(nSteps>=1 && ...
                all(isfinite(loadHistory(:))), ...
                'Load history must have at least two finite columns.');

            fixedDOF = unique(fixedDOF(:)');

            assert(all(fixedDOF>=1 & ...
                fixedDOF<=nDOF & ...
                fixedDOF==fix(fixedDOF)), ...
                'Invalid fixed DOF indices.');

            freeDOF = setdiff(1:nDOF,fixedDOF);

            assert(~isempty(freeDOF), ...
                'No unconstrained DOFs.');

            %% Boundary condition method

            bcMethod = lower(char(options.bcMethod));

            penalty = [];

            responseEffective = response;

            activeDOF = freeDOF;
            fixedForIteration = fixedDOF;

            switch bcMethod

                case 'elimination'

                    assert(all(options.prescribedValues==0), ...
                        ['Nonzero prescribedValues require ' ...
                         'bcMethod=penalty.']);

                case 'penalty'

                    % Initial tangent is used as the reference
                    % stiffness for the penalty parameter.

                    [~,Kreference] = response(zeros(nDOF,1));

                    penalty = ...
                        fea.assembly.PenaltyMethod.create( ...
                        Kreference, ...
                        fixedDOF, ...
                        options.prescribedValues, ...
                        options.penaltyFactor);

                    % Incorporate penalty force and tangent
                    % into the existing response callback.

                    responseEffective = @(qTrial) ...
                        fea.solvers.NewtonRaphson.penaltyResponse( ...
                        response,penalty,qTrial);

                    % Penalty solves the complete DOF system.
                    activeDOF = 1:nDOF;
                    fixedForIteration = [];

                otherwise

                    error( ...
                        'Unknown boundary condition method: %s', ...
                        bcMethod);
            end

            %% Initialize solution

            q = zeros(nDOF,1);

            if strcmp(bcMethod,'penalty')
                q = penalty.uBar;
            end

            qHistory = zeros(nDOF,nSteps+1);
            qHistory(:,1) = q;

            corrections = zeros(nSteps,1);
            residualHistory = cell(nSteps,1);

            acceptedSubsteps = zeros(nSteps,1);
            totalCutbacks = zeros(nSteps,1);

            %% Optional restart from completed load step

            startStep = 1;
            if ~isempty(options.restartFrom)
                saved = load(options.restartFrom,'checkpoint');
                assert(isfield(saved,'checkpoint'), ...
                    'Restart file does not contain a checkpoint.');
                cp = saved.checkpoint;
                assert(isequal(cp.loadHistory,loadHistory) && ...
                    isequal(cp.fixedDOF,fixedDOF) && ...
                    strcmp(cp.bcMethod,bcMethod) && ...
                    strcmp(cp.method,method) && ...
                    isequal(cp.restartKey,options.restartKey) && ...
                    isequal(cp.penaltyFactor,options.penaltyFactor) && ...
                    isequal(cp.prescribedValues,options.prescribedValues), ...
                    'Restart settings, loads or model key do not match.');
                assert(cp.completedStep>=0 && cp.completedStep<=nSteps, ...
                    'Invalid restart step.');
                startStep = cp.completedStep+1;
                q = cp.q;
                qHistory = cp.qHistory;
                corrections = cp.corrections;
                residualHistory = cp.residualHistory;
                acceptedSubsteps = cp.acceptedSubsteps;
                totalCutbacks = cp.totalCutbacks;
                fprintf('RESTART: resuming from completed step %d\n',cp.completedStep);
            end

            %% Iteration output header

            if options.verbose

                fprintf( ...
                    'Step   Iteration    Relative Residual\n');

                fprintf( ...
                    '-------------------------------------\n');
            end

            %% ========================================================
            % LOAD STEPS
            % =========================================================

            for step = startStep:nSteps

                previousLoad = loadHistory(:,step);

                pending = {loadHistory(:,step+1)};
                currentLoad = previousLoad;

                rrAll = [];

                while ~isempty(pending)

                    target = pending{1};

                    [ok,qTrial,rr,nCorr] = ...
                        fea.solvers.NewtonRaphson.iterate( ...
                        responseEffective, ...
                        q, ...
                        target, ...
                        activeDOF, ...
                        fixedForIteration, ...
                        options, ...
                        method, ...
                        step);

                    if ok

                        %% Accept converged state

                        q = qTrial;
                        currentLoad = target;

                        pending(1) = [];

                        rrAll = [rrAll;rr(:)]; %#ok<AGROW>

                        corrections(step) = ...
                            corrections(step) + nCorr;

                        acceptedSubsteps(step) = ...
                            acceptedSubsteps(step) + 1;

                        if options.verbose

                            fprintf( ...
                                '%3d    CONVERGED    %d iterations\n', ...
                                step,nCorr);
                        end

                    else

                        %% Failed increment: optional bisection

                        if ~options.cutback || ...
                                totalCutbacks(step)>=options.maxCutbacks

                            error('FEA:Newton:NoConvergence', ...
                                ['Load step %d failed with %s ' ...
                                 'after %d iterations ' ...
                                 '(residual %.3e).'], ...
                                step, ...
                                method, ...
                                options.maxIter, ...
                                rr(end));
                        end

                        midpoint = 0.5*(currentLoad+target);

                        if norm(midpoint-currentLoad)==0

                            error('FEA:Newton:CutbackLimit', ...
                                'Load increment too small to bisect.');
                        end

                        % Retry midpoint first, then original target.
                        pending = [{midpoint},pending]; %#ok<AGROW>

                        totalCutbacks(step) = ...
                            totalCutbacks(step) + 1;

                        if options.verbose

                            fprintf( ...
                                '%3d    CUTBACK      load increment bisected\n', ...
                                step);
                        end
                    end
                end

                qHistory(:,step+1) = q;
                residualHistory{step} = rrAll;

                % Checkpoint only after full load-step convergence.
                if ~isempty(options.checkpointFile) && ...
                        (mod(step,options.saveEvery)==0 || step==nSteps)
                    checkpoint = struct( ...
                        'completedStep',step, 'q',q, ...
                        'qHistory',qHistory, ...
                        'corrections',corrections, ...
                        'residualHistory',{residualHistory}, ...
                        'acceptedSubsteps',acceptedSubsteps, ...
                        'totalCutbacks',totalCutbacks, ...
                        'loadHistory',loadHistory, ...
                        'fixedDOF',fixedDOF, ...
                        'bcMethod',bcMethod, ...
                        'method',method, ...
                        'penaltyFactor',options.penaltyFactor, ...
                        'prescribedValues',options.prescribedValues, ...
                        'restartKey',options.restartKey);
                    save(options.checkpointFile,'checkpoint');
                end
            end

            %% Return results

            result.q = q;
            result.qHistory = qHistory;

            result.corrections = corrections;
            result.residualHistory = residualHistory;

            result.fixedDOF = fixedDOF;
            result.freeDOF = freeDOF;

            result.method = method;

            result.acceptedSubsteps = acceptedSubsteps;
            result.cutbacks = totalCutbacks;

            result.bcMethod = bcMethod;
            result.penalty = penalty;

        end
    end

    methods (Static,Access=private)

        %% ============================================================
        % PENALTY RESPONSE
        % =============================================================

        function [fint,Kt] = ...
                penaltyResponse(response,penalty,q)

            [fint,Kt] = response(q);

            [fint,Kt] = ...
                fea.assembly.PenaltyMethod.applyNonlinear( ...
                fint,Kt,q,penalty);
        end

        %% ============================================================
        % NEWTON ITERATIONS
        % =============================================================

        function [ok,q,rr,nCorr] = iterate( ...
                response,qStart,Fext,free,fixed,opts,method,step)

            q = qStart;

            nf = numel(free);

            rr = zeros(opts.maxIter+1,1);

            ok = false;
            nCorr = 0;

            Kapprox = [];

            prevQ = [];
            prevF = [];

            %% Iterative equilibrium solution

            for iter = 0:opts.maxIter

                [fint,Kexact] = response(q);

                R = Fext-fint;

                relR = norm(R(free)) / ...
                    max(1,norm(Fext(free)));

                rr(iter+1) = relR;

                %% Print iteration

                if opts.verbose && ...
                    (mod(iter,opts.printEvery)==0 || ...
                    relR<opts.tolResidual)

                    fprintf( ...
                        '%3d       %3d        %.6e\n', ...
                        step,iter,relR);
                end

                %% Convergence check

                if relR<opts.tolResidual

                    ok = true;
                    nCorr = iter;

                    rr = rr(1:iter+1);
                    return;
                end

                if iter==opts.maxIter
                    break;
                end

                %% Select nonlinear tangent method

                switch method

                    case 'full'

                        Kuse = Kexact(free,free);

                    case 'modified'

                        if isempty(Kapprox)
                            Kapprox = Kexact(free,free);
                        end

                        Kuse = Kapprox;

                    case 'numerical'

                        Kuse = zeros(nf,nf);

                        for j = 1:nf

                            qplus = q;
                            qminus = q;

                            h = opts.fdStep * ...
                                max(1,abs(q(free(j))));

                            qplus(free(j)) = ...
                                qplus(free(j)) + h;

                            qminus(free(j)) = ...
                                qminus(free(j)) - h;

                            fplus = response(qplus);
                            fminus = response(qminus);

                            Kuse(:,j) = ...
                                (fplus(free)-fminus(free))/(2*h);
                        end

                    case {'secant','broyden'}

                        if isempty(Kapprox)

                            Kapprox = Kexact(free,free);

                        elseif ~isempty(prevQ)

                            s = q(free)-prevQ;
                            y = fint(free)-prevF;

                            denom = s'*s;

                            if denom>eps(max(1,norm(s)^2))

                                Kapprox = Kapprox + ...
                                    ((y-Kapprox*s)*s')/denom;
                            end
                        end

                        Kuse = Kapprox;

                end

                %% Solve displacement correction

                if any(~isfinite(Kuse(:))) || ...
                        rcond(full(Kuse))<1e-14
                    break;
                end

                dq = Kuse \ R(free);

                if any(~isfinite(dq))
                    break;
                end

                %% Update configuration

                prevQ = q(free);
                prevF = fint(free);

                q(free) = q(free)+dq;
                q(fixed) = 0;
            end

            rr = rr(1:iter+1);

            % Do not commit a failed trial configuration.
            q = qStart;
        end
    end
end


% classdef NewtonRaphson
%     %NEWTONRAPHSON Incremental nonlinear equilibrium solver.
%     %
%     % response(q) returns [fint,Kt].
%     %
%     % Methods:
%     %   full      - exact tangent every iteration
%     %   modified  - tangent fixed within each load increment
%     %   secant    - Broyden rank-one tangent update
%     %   broyden   - alias for secant
%     %   numerical - centered finite-difference tangent
%     %
%     % Boundary conditions:
%     %   elimination (default)
%     %   penalty
%     %
%     % Optional load-step bisection and iteration output.
% 
%     methods (Static)
% 
%         function result = solve( ...
%                 response,loadHistory,fixedDOF,options)
% 
%             if nargin < 4 || isempty(options)
%                 options = struct();
%             end
% 
%             %% Default settings
% 
%             defaults = struct( ...
%                 'tolResidual',1e-10, ...
%                 'maxIter',25, ...
%                 'verbose',false, ...
%                 'printEvery',1, ...
%                 'method','full', ...
%                 'cutback',false, ...
%                 'maxCutbacks',8, ...
%                 'fdStep',1e-6, ...
%                 'bcMethod','elimination', ...
%                 'penaltyFactor',1e6, ...
%                 'prescribedValues',0);
% 
%             names = fieldnames(defaults);
% 
%             for i = 1:numel(names)
% 
%                 if ~isfield(options,names{i})
%                     options.(names{i}) = defaults.(names{i});
%                 end
%             end
% 
%             %% Solver method
% 
%             method = lower(char(options.method));
% 
%             assert(ismember(method, ...
%                 {'full','modified','secant','broyden','numerical'}), ...
%                 'Unknown nonlinear solver method.');
% 
%             assert(isscalar(options.verbose) && ...
%                 (islogical(options.verbose) || ...
%                 isnumeric(options.verbose)), ...
%                 'options.verbose must be true or false.');
% 
%             assert(isscalar(options.printEvery) && ...
%                 options.printEvery>=1 && ...
%                 options.printEvery==fix(options.printEvery), ...
%                 'printEvery must be a positive integer.');
% 
%             assert(isscalar(options.maxIter) && ...
%                 options.maxIter>=1 && ...
%                 options.maxIter==fix(options.maxIter), ...
%                 'Invalid maxIter.');
% 
%             assert(isscalar(options.maxCutbacks) && ...
%                 options.maxCutbacks>=0 && ...
%                 options.maxCutbacks==fix(options.maxCutbacks), ...
%                 'Invalid maxCutbacks.');
% 
%             assert(options.fdStep>0 && ...
%                 isfinite(options.fdStep), ...
%                 'fdStep must be positive.');
% 
%             %% Degrees of freedom
% 
%             nDOF = size(loadHistory,1);
%             nSteps = size(loadHistory,2)-1;
% 
%             assert(nSteps>=1 && ...
%                 all(isfinite(loadHistory(:))), ...
%                 'Load history must have at least two finite columns.');
% 
%             fixedDOF = unique(fixedDOF(:)');
% 
%             assert(all(fixedDOF>=1 & ...
%                 fixedDOF<=nDOF & ...
%                 fixedDOF==fix(fixedDOF)), ...
%                 'Invalid fixed DOF indices.');
% 
%             freeDOF = setdiff(1:nDOF,fixedDOF);
% 
%             assert(~isempty(freeDOF), ...
%                 'No unconstrained DOFs.');
% 
%             %% Boundary condition method
% 
%             bcMethod = lower(char(options.bcMethod));
% 
%             penalty = [];
% 
%             responseEffective = response;
% 
%             activeDOF = freeDOF;
%             fixedForIteration = fixedDOF;
% 
%             switch bcMethod
% 
%                 case 'elimination'
% 
%                     assert(all(options.prescribedValues==0), ...
%                         ['Nonzero prescribedValues require ' ...
%                          'bcMethod=penalty.']);
% 
%                 case 'penalty'
% 
%                     % Initial tangent is used as the reference
%                     % stiffness for the penalty parameter.
% 
%                     [~,Kreference] = response(zeros(nDOF,1));
% 
%                     penalty = ...
%                         fea.solvers.PenaltyMethod.create( ...
%                         Kreference, ...
%                         fixedDOF, ...
%                         options.prescribedValues, ...
%                         options.penaltyFactor);
% 
%                     % Incorporate penalty force and tangent
%                     % into the existing response callback.
% 
%                     responseEffective = @(qTrial) ...
%                         fea.solvers.NewtonRaphson.penaltyResponse( ...
%                         response,penalty,qTrial);
% 
%                     % Penalty solves the complete DOF system.
%                     activeDOF = 1:nDOF;
%                     fixedForIteration = [];
% 
%                 otherwise
% 
%                     error( ...
%                         'Unknown boundary condition method: %s', ...
%                         bcMethod);
%             end
% 
%             %% Initialize solution
% 
%             q = zeros(nDOF,1);
% 
%             if strcmp(bcMethod,'penalty')
%                 q = penalty.uBar;
%             end
% 
%             qHistory = zeros(nDOF,nSteps+1);
%             qHistory(:,1) = q;
% 
%             corrections = zeros(nSteps,1);
%             residualHistory = cell(nSteps,1);
% 
%             acceptedSubsteps = zeros(nSteps,1);
%             totalCutbacks = zeros(nSteps,1);
% 
%             %% Iteration output header
% 
%             if options.verbose
% 
%                 fprintf( ...
%                     'Step   Iteration    Relative Residual\n');
% 
%                 fprintf( ...
%                     '-------------------------------------\n');
%             end
% 
%             %% ========================================================
%             % LOAD STEPS
%             % =========================================================
% 
%             for step = 1:nSteps
% 
%                 previousLoad = loadHistory(:,step);
% 
%                 pending = {loadHistory(:,step+1)};
%                 currentLoad = previousLoad;
% 
%                 rrAll = [];
% 
%                 while ~isempty(pending)
% 
%                     target = pending{1};
% 
%                     [ok,qTrial,rr,nCorr] = ...
%                         fea.solvers.NewtonRaphson.iterate( ...
%                         responseEffective, ...
%                         q, ...
%                         target, ...
%                         activeDOF, ...
%                         fixedForIteration, ...
%                         options, ...
%                         method, ...
%                         step);
% 
%                     if ok
% 
%                         %% Accept converged state
% 
%                         q = qTrial;
%                         currentLoad = target;
% 
%                         pending(1) = [];
% 
%                         rrAll = [rrAll;rr(:)]; %#ok<AGROW>
% 
%                         corrections(step) = ...
%                             corrections(step) + nCorr;
% 
%                         acceptedSubsteps(step) = ...
%                             acceptedSubsteps(step) + 1;
% 
%                         if options.verbose
% 
%                             fprintf( ...
%                                 '%3d    CONVERGED    %d iterations\n', ...
%                                 step,nCorr);
%                         end
% 
%                     else
% 
%                         %% Failed increment: optional bisection
% 
%                         if ~options.cutback || ...
%                                 totalCutbacks(step)>=options.maxCutbacks
% 
%                             error('FEA:Newton:NoConvergence', ...
%                                 ['Load step %d failed with %s ' ...
%                                  'after %d iterations ' ...
%                                  '(residual %.3e).'], ...
%                                 step, ...
%                                 method, ...
%                                 options.maxIter, ...
%                                 rr(end));
%                         end
% 
%                         midpoint = 0.5*(currentLoad+target);
% 
%                         if norm(midpoint-currentLoad)==0
% 
%                             error('FEA:Newton:CutbackLimit', ...
%                                 'Load increment too small to bisect.');
%                         end
% 
%                         % Retry midpoint first, then original target.
%                         pending = [{midpoint},pending]; %#ok<AGROW>
% 
%                         totalCutbacks(step) = ...
%                             totalCutbacks(step) + 1;
% 
%                         if options.verbose
% 
%                             fprintf( ...
%                                 '%3d    CUTBACK      load increment bisected\n', ...
%                                 step);
%                         end
%                     end
%                 end
% 
%                 qHistory(:,step+1) = q;
%                 residualHistory{step} = rrAll;
%             end
% 
%             %% Return results
% 
%             result.q = q;
%             result.qHistory = qHistory;
% 
%             result.corrections = corrections;
%             result.residualHistory = residualHistory;
% 
%             result.fixedDOF = fixedDOF;
%             result.freeDOF = freeDOF;
% 
%             result.method = method;
% 
%             result.acceptedSubsteps = acceptedSubsteps;
%             result.cutbacks = totalCutbacks;
% 
%             result.bcMethod = bcMethod;
%             result.penalty = penalty;
% 
%         end
%     end
% 
%     methods (Static,Access=private)
% 
%         %% ============================================================
%         % PENALTY RESPONSE
%         % =============================================================
% 
%         function [fint,Kt] = ...
%                 penaltyResponse(response,penalty,q)
% 
%             [fint,Kt] = response(q);
% 
%             [fint,Kt] = ...
%                 fea.solvers.PenaltyMethod.applyNonlinear( ...
%                 fint,Kt,q,penalty);
%         end
% 
%         %% ============================================================
%         % NEWTON ITERATIONS
%         % =============================================================
% 
%         function [ok,q,rr,nCorr] = iterate( ...
%                 response,qStart,Fext,free,fixed,opts,method,step)
% 
%             q = qStart;
% 
%             nf = numel(free);
% 
%             rr = zeros(opts.maxIter+1,1);
% 
%             ok = false;
%             nCorr = 0;
% 
%             Kapprox = [];
% 
%             prevQ = [];
%             prevF = [];
% 
%             %% Iterative equilibrium solution
% 
%             for iter = 0:opts.maxIter
% 
%                 [fint,Kexact] = response(q);
% 
%                 R = Fext-fint;
% 
%                 relR = norm(R(free)) / ...
%                     max(1,norm(Fext(free)));
% 
%                 rr(iter+1) = relR;
% 
%                 %% Print iteration
% 
%                 if opts.verbose && ...
%                     (mod(iter,opts.printEvery)==0 || ...
%                     relR<opts.tolResidual)
% 
%                     fprintf( ...
%                         '%3d       %3d        %.6e\n', ...
%                         step,iter,relR);
%                 end
% 
%                 %% Convergence check
% 
%                 if relR<opts.tolResidual
% 
%                     ok = true;
%                     nCorr = iter;
% 
%                     rr = rr(1:iter+1);
%                     return;
%                 end
% 
%                 if iter==opts.maxIter
%                     break;
%                 end
% 
%                 %% Select nonlinear tangent method
% 
%                 switch method
% 
%                     case 'full'
% 
%                         Kuse = Kexact(free,free);
% 
%                     case 'modified'
% 
%                         if isempty(Kapprox)
%                             Kapprox = Kexact(free,free);
%                         end
% 
%                         Kuse = Kapprox;
% 
%                     case 'numerical'
% 
%                         Kuse = zeros(nf,nf);
% 
%                         for j = 1:nf
% 
%                             qplus = q;
%                             qminus = q;
% 
%                             h = opts.fdStep * ...
%                                 max(1,abs(q(free(j))));
% 
%                             qplus(free(j)) = ...
%                                 qplus(free(j)) + h;
% 
%                             qminus(free(j)) = ...
%                                 qminus(free(j)) - h;
% 
%                             fplus = response(qplus);
%                             fminus = response(qminus);
% 
%                             Kuse(:,j) = ...
%                                 (fplus(free)-fminus(free))/(2*h);
%                         end
% 
%                     case {'secant','broyden'}
% 
%                         if isempty(Kapprox)
% 
%                             Kapprox = Kexact(free,free);
% 
%                         elseif ~isempty(prevQ)
% 
%                             s = q(free)-prevQ;
%                             y = fint(free)-prevF;
% 
%                             denom = s'*s;
% 
%                             if denom>eps(max(1,norm(s)^2))
% 
%                                 Kapprox = Kapprox + ...
%                                     ((y-Kapprox*s)*s')/denom;
%                             end
%                         end
% 
%                         Kuse = Kapprox;
% 
%                 end
% 
%                 %% Solve displacement correction
% 
%                 if any(~isfinite(Kuse(:))) || ...
%                         rcond(full(Kuse))<1e-14
%                     break;
%                 end
% 
%                 dq = Kuse \ R(free);
% 
%                 if any(~isfinite(dq))
%                     break;
%                 end
% 
%                 %% Update configuration
% 
%                 prevQ = q(free);
%                 prevF = fint(free);
% 
%                 q(free) = q(free)+dq;
%                 q(fixed) = 0;
%             end
% 
%             rr = rr(1:iter+1);
% 
%             % Do not commit a failed trial configuration.
%             q = qStart;
%         end
%     end
% end
