
%% ================================================================
% MECH 6140 - HW3 SINGLE Q4 ELEMENT VERIFICATION
%
% Linear FEA vs Total Lagrangian Nonlinear FEA
%
% Supports:
%   - StVK and LogElastic materials
%   - Full, Modified, Secant, Numerical Newton methods
%   - Elimination and Penalty boundary conditions
%   - Optional load-step bisection
%   - Analytical reference solution
%   - Tangent stiffness verification
%   - Unified field contour plotting
%
% Q4 node ordering:
%   4 -------- 3
%   |          |
%   |          |
%   1 -------- 2
%
% DOF order:
% [ux1 uy1 ux2 uy2 ux3 uy3 ux4 uy4]
% ================================================================

clearvars;
clc;
close all;

format long g;

%% ================================================================
% 1. USER INPUTS
% ================================================================

% Geometry
L = 1;
H = 1;
thickness = 1;

% Material
Eyoung = 70000;
nu = 0.30;

materialModel = 'StVK';
% materialModel = 'LogElastic';

% Loading
Ptotal = 25000;
nSteps = 10;

% Boundary conditions
boundaryCase = 'freeY';
% boundaryCase = 'fixedY';

bcMethod = 'elimination';
% bcMethod = 'penalty';

penaltyFactor = 1e6;

% Nonlinear solution method
nonlinearMethod = 'full';
% nonlinearMethod = 'modified';
% nonlinearMethod = 'secant';
% nonlinearMethod = 'numerical';

% Newton settings
maxIter = 50;
tolResidual = 1e-10;

showNewtonIterations = true;
printEvery = 1;

% Automatic load-step bisection
useLoadBisection = false;
maxCutbacks = 8;

% Verification
tolVerification = 1e-8;
verifyTangent = true;

% Plot controls
plotForceDisplacement = true;
plotConvergence = true;
plotShapeComparison = true;

plotLinearContours = true;
plotNonlinearContours = true;

contourPalette = 'jet';
contourScale = 1;

showNodes = false;

%% ================================================================
% 2. SINGLE Q4 ELEMENT
% ================================================================

X = [ ...
    0 0;
    L 0;
    L H;
    0 H];

conn = [1 2 3 4];

nNodes = size(X,1);
nDOF = 2*nNodes;

%% ================================================================
% 3. MATERIAL SELECTION
% ================================================================

switch lower(materialModel)

    case 'stvk'

        material = ...
            fea.materials.LinearElastic( ...
            Eyoung,nu);

    case 'logelastic'

        material = ...
            fea.materials.LogElastic( ...
            Eyoung,nu);

    otherwise

        error('Unknown material model.');
end

%% ================================================================
% 4. Q4 ELEMENT AND TL FORMULATION
% ================================================================

element = fea.elements.Quad4( ...
    X,material,thickness);

% Total Lagrangian response is now provided by the element directly.

%% ================================================================
% 5. BOUNDARY CONDITIONS
% ================================================================

switch lower(boundaryCase)

    case 'freey'

        % Fix UX on left edge.
        % Fix UY at node 1 to prevent rigid translation.
        %
        % Transverse contraction is allowed.

        fixedDOF = [1 2 7];

        caseDescription = ...
            'Transverse contraction allowed';

    case 'fixedy'

        % Fix UX on left edge.
        % Fix all vertical DOFs.

        fixedDOF = [1 2 4 6 7 8];

        caseDescription = ...
            'Transverse contraction constrained';

    otherwise

        error('Unknown boundary-condition case.');
end

freeDOF = setdiff(1:nDOF,fixedDOF);

%% ================================================================
% 6. EXTERNAL FORCE AND LOAD HISTORY
% ================================================================

% Apply equal horizontal loads at nodes 2 and 3.

rightDOF = [3 5];

Ftotal = zeros(nDOF,1);

Ftotal(rightDOF) = ...
    Ptotal/numel(rightDOF);

loadFactors = linspace(0,1,nSteps+1);

Psteps = Ptotal*loadFactors;

loadHistory = Ftotal*loadFactors;

%% ================================================================
% 7. LINEAR SOLUTION
% ================================================================

Klinear = element.Gauss_Element_stiffness_cal();

linearOpts = struct( ...
    'bcMethod',bcMethod, ...
    'penaltyFactor',penaltyFactor);

LIN = fea.solvers.LinearStatic.solve( ...
    Klinear,loadHistory,fixedDOF,linearOpts);

qLinear = LIN.q;
qLinearHistory = LIN.qHistory;

uLinear = mean( ...
    qLinearHistory(rightDOF,:),1);

%% ================================================================
% 8. NONLINEAR SOLUTION
% ================================================================

response = @(qTrial) q4Response( ...
    element,qTrial);

nonlinearOpts = struct( ...
    'method',nonlinearMethod, ...
    'maxIter',maxIter, ...
    'tolResidual',tolResidual, ...
    'verbose',showNewtonIterations, ...
    'printEvery',printEvery, ...
    'cutback',useLoadBisection, ...
    'maxCutbacks',maxCutbacks, ...
    'bcMethod',bcMethod, ...
    'penaltyFactor',penaltyFactor);

NL = fea.solvers.NewtonRaphson.solve( ...
    response,loadHistory,fixedDOF,nonlinearOpts);

qNonlinear = NL.q;
qHistory = NL.qHistory;

uNonlinear = mean( ...
    qHistory(rightDOF,:),1);

%% ================================================================
% 9. ANALYTICAL REFERENCE SOLUTION
% ================================================================

% Homogeneous deformation:
%
% F = diag(lambdaX,lambdaY)
%
% This analytical solution assumes ideal zero-displacement
% constraints. Penalty BCs approximate these constraints.

lambdaL = material.lambda_L;
mu = material.mu;

C11 = lambdaL+2*mu;

switch lower(materialModel)

    case 'stvk'

        E11fun = @(lx) 0.5*(lx.^2-1);

        if strcmpi(boundaryCase,'freeY')

            E22fun = @(lx) ...
                -(lambdaL/C11)*E11fun(lx);

            lambdaYfun = @(lx) ...
                sqrt(1+2*E22fun(lx));

        else

            E22fun = @(lx) zeros(size(lx));

            lambdaYfun = @(lx) ones(size(lx));

        end

        Pfun = @(lx) ...
            lx.*(C11*E11fun(lx) + ...
            lambdaL*E22fun(lx))*H*thickness;

    case 'logelastic'

        if strcmpi(boundaryCase,'freeY')

            ratio = lambdaL/C11;

            lambdaYfun = @(lx) lx.^(-ratio);

            effectiveModulus = ...
                C11-lambdaL^2/C11;

        else

            lambdaYfun = @(lx) ones(size(lx));

            effectiveModulus = C11;

        end

        Pfun = @(lx) ...
            lambdaYfun(lx).* ...
            effectiveModulus.* ...
            log(lx)*H*thickness;
end

%% Solve analytical equilibrium at each load step

lambdaGrid = linspace(1,5,4001);
Pgrid = Pfun(lambdaGrid);

uAnalytical = zeros(1,nSteps+1);
vAnalytical = zeros(1,nSteps+1);

for step = 2:nSteps+1

    targetForce = Psteps(step);

    crossing = find( ...
        Pgrid>=targetForce,1,'first');

    if isempty(crossing) || crossing<2

        error(['Analytical equilibrium not found ' ...
               'at load step %d.'],step-1);
    end

    lambdaX = fzero( ...
        @(lx) Pfun(lx)-targetForce, ...
        [lambdaGrid(crossing-1), ...
         lambdaGrid(crossing)]);

    lambdaY = lambdaYfun(lambdaX);

    uAnalytical(step) = (lambdaX-1)*L;

    vAnalytical(step) = (lambdaY-1)*H;
end

%% ================================================================
% 10. FINAL RESULTS
% ================================================================

uNewton = mean(qNonlinear(rightDOF));

vNewton = qNonlinear(6)-qNonlinear(4);

uExact = uAnalytical(end);
vExact = vAnalytical(end);

% Physical internal forces and total tangent
[fint,KL,KNL] = ...
    element.Element_Response_cal(qNonlinear);

Ktotal = KL+KNL;

Rphysical = Ftotal-fint;

% Penalty equilibrium includes the penalty internal force.
if strcmpi(bcMethod,'penalty')

    Kp = NL.penalty.Kp;
    uBar = NL.penalty.uBar;

    Rsolver = Rphysical - Kp*(qNonlinear-uBar);

    relativeResidual = ...
        norm(Rsolver)/max(1,norm(Ftotal));

    constraintError = norm( ...
        qNonlinear(fixedDOF)-uBar(fixedDOF));

else

    relativeResidual = ...
        norm(Rphysical(freeDOF))/ ...
        max(1,norm(Ftotal(freeDOF)));

    constraintError = norm(qNonlinear(fixedDOF));
end

uError = abs(uNewton-uExact)/max(1,abs(uExact));

vError = abs(vNewton-vExact)/max(1,abs(vExact));

%% Print results

fprintf('\n========================================\n');
fprintf(' SINGLE Q4 VERIFICATION\n');
fprintf('========================================\n');

fprintf('Material          : %s\n',materialModel);
fprintf('BC configuration  : %s\n',boundaryCase);
fprintf('BC method         : %s\n',bcMethod);
fprintf('Nonlinear method  : %s\n',nonlinearMethod);
fprintf('Load steps        : %d\n',nSteps);
fprintf('Applied force     : %.6f\n',Ptotal);

fprintf('\n--- DISPLACEMENT COMPARISON ---\n');

fprintf('Analytical UX     : %.12f\n',uExact);
fprintf('Linear UX         : %.12f\n',uLinear(end));
fprintf('Nonlinear UX      : %.12f\n',uNewton);

fprintf('Analytical dY     : %.12f\n',vExact);
fprintf('Nonlinear dY      : %.12f\n',vNewton);

fprintf('\n--- VERIFICATION ---\n');

fprintf('Relative UX error : %.6e\n',uError);
fprintf('Relative dY error : %.6e\n',vError);
fprintf('Equilibrium error : %.6e\n',relativeResidual);
fprintf('Constraint error  : %.6e\n',constraintError);

fprintf('\n--- NEWTON ITERATIONS ---\n');

for step = 1:nSteps

    fprintf('Step %2d : %d iterations\n', ...
        step,NL.corrections(step));
end

%% Verification assessment

if strcmpi(bcMethod,'elimination')

    if uError<tolVerification && ...
       vError<tolVerification && ...
       relativeResidual<tolResidual

        fprintf('\nPASS: Analytical verification.\n');

    else

        warning(['Analytical or equilibrium verification ' ...
                 'exceeded tolerance.']);
    end

else

    fprintf(['\nPenalty BC active: analytical differences ' ...
             'include finite penalty-constraint error.\n']);

    if relativeResidual>=tolResidual
        warning('Penalty equilibrium residual exceeds tolerance.');
    end
end

%% ================================================================
% 11. NUMERICAL TANGENT VERIFICATION
% ================================================================

if verifyTangent

    hFD = 1e-6;

    Kfd = zeros(nDOF,nDOF);

    for j = 1:nDOF

        dq = zeros(nDOF,1);
        dq(j) = hFD;

        fPlus = ...
            element.Element_Response_cal(qNonlinear+dq);

        fMinus = ...
            element.Element_Response_cal(qNonlinear-dq);

        Kfd(:,j) = ...
            (fPlus-fMinus)/(2*hFD);
    end

    tangentError = ...
        norm(Ktotal-Kfd,'fro')/ ...
        max(1,norm(Kfd,'fro'));

    fprintf('\n--- ELEMENT TANGENT CHECK ---\n');

    fprintf('Relative FD error : %.6e\n', ...
        tangentError);

    if tangentError<1e-5

        fprintf('PASS: Element tangent check.\n');

    else

        warning('Element tangent verification failed.');
    end
end

%% ================================================================
% 12. FORCE-DISPLACEMENT COMPARISON
% ================================================================

if plotForceDisplacement

    figure('Name','Single Q4 Force Displacement', ...
           'Color','w');

    hA = plot(uAnalytical,Psteps, ...
        'k-','LineWidth',1.5);
    hold on;

    hN = plot(uNonlinear,Psteps, ...
        'b-o','LineWidth',1.3, ...
        'MarkerSize',5);

    hL = plot(uLinear,Psteps, ...
        'r--s','LineWidth',1.3, ...
        'MarkerSize',4);

    xlabel('Right-edge displacement u_x');
    ylabel('Applied force P');

    title('Single Q4 Force-Displacement', ...
        'FontWeight','normal');

    legend([hA hN hL], ...
        {'Analytical','Nonlinear TL','Linear'}, ...
        'Location','best');

    grid on;
    box on;
end

%% ================================================================
% 13. NEWTON CONVERGENCE
% ================================================================

if plotConvergence

    figure('Name','Newton Convergence','Color','w');
    hold on;

    for step = 1:nSteps

        rr = NL.residualHistory{step};

        semilogy( ...
            0:numel(rr)-1, ...
            max(rr,eps), ...
            'o-', ...
            'LineWidth',1.1, ...
            'MarkerSize',3, ...
            'DisplayName',sprintf('Load Step %d',step));
    end

    xlabel('Newton iteration');
    ylabel('Relative equilibrium residual');

    title('Newton Convergence', ...
        'FontWeight','normal');

    legend('Location','bestoutside');

    grid on;
    box on;

    if any(NL.cutbacks>0)
        warning(['Cutbacks occurred. Residual history ' ...
                 'contains multiple Newton attempts.']);
    end
end

%% ================================================================
% 14. DEFORMED SHAPE COMPARISON
% ================================================================

uLinearXY = reshape(qLinear,2,nNodes).';

uNonlinearXY = reshape(qNonlinear,2,nNodes).';

Xlinear = X + contourScale*uLinearXY;

Xnonlinear = X + contourScale*uNonlinearXY;

nodeLoop = [1 2 3 4 1];

if plotShapeComparison

    figure('Name','Single Q4 Shape Comparison', ...
           'Color','w');

    hRef = plot( ...
        X(nodeLoop,1), ...
        X(nodeLoop,2), ...
        'k--o','LineWidth',1.2);

    hold on;

    hLin = plot( ...
        Xlinear(nodeLoop,1), ...
        Xlinear(nodeLoop,2), ...
        'r--s','LineWidth',1.2);

    hNL = plot( ...
        Xnonlinear(nodeLoop,1), ...
        Xnonlinear(nodeLoop,2), ...
        'b-o','LineWidth',1.5);

    xlabel('X');
    ylabel('Y');

    title('Reference, Linear and Nonlinear Shapes', ...
        'FontWeight','normal');

    legend([hRef hLin hNL], ...
        {'Reference','Linear','Nonlinear TL'}, ...
        'Location','best');

    axis equal;
    grid on;
    box on;
end

%% ================================================================
% 15. UNIFIED FIELD CONTOUR PLOTTING
% ================================================================

% The generic FieldPlotter receives arbitrary scalar fields.
% Here we supply displacement components and magnitude.

fieldsLinear = [ ...
    uLinearXY(:,1), ...
    uLinearXY(:,2), ...
    hypot(uLinearXY(:,1),uLinearXY(:,2))];

fieldsNonlinear = [ ...
    uNonlinearXY(:,1), ...
    uNonlinearXY(:,2), ...
    hypot(uNonlinearXY(:,1),uNonlinearXY(:,2))];

plotOptions = struct();

plotOptions.location = 'nodal';
plotOptions.colormap = contourPalette;
plotOptions.scale = contourScale;

plotOptions.showNodes = showNodes;

plotOptions.fieldNames = ...
    {'u_x','u_y','|u|'};

%% Linear contour fields

if plotLinearContours

    plotOptions.q = qLinear;

    plotOptions.figureName = ...
        'Linear Displacement Fields';

    fea.results.FieldPlotter.plot( ...
        X,conn,fieldsLinear,plotOptions);
end

%% Nonlinear contour fields

if plotNonlinearContours

    plotOptions.q = qNonlinear;

    plotOptions.figureName = ...
        'Nonlinear Displacement Fields';

    fea.results.FieldPlotter.plot( ...
        X,conn,fieldsNonlinear,plotOptions);
end

%% ================================================================
% LOCAL ELEMENT RESPONSE ADAPTER
% ================================================================

function [fint,Kt] = q4Response( ...
    element,q)

    [fint,KL,KNL] = ...
        element.Element_Response_cal(q);

    Kt = KL+KNL;
end



% %% HW3_Verify_Solvers_OOP.m
% % MECH 6140 - Single Q4 element: TL Newton-Raphson vs HW2 linear response
% % Run from the project root containing the +fea_hw2 package.
% %
% % Q4 node ordering: 1=(0,0), 2=(L,0), 3=(L,H), 4=(0,H).
% % DOF ordering: [ux1 uy1 ux2 uy2 ux3 uy3 ux4 uy4].
% % Plane strain; St Venant-Kirchhoff material; 2x2 Gaussian integration.
% % Equal horizontal nodal loads on right edge (nodes 2 and 3).
% % No penalty or global assembly. Solver algorithms now live in +solvers.
% 
% clearvars; clc; close all;
% format long g;
% 
% %% 1. INPUTS
% L = 1; H = 1; thickness = 1;
% Eyoung = 70000; nu = 0.30;
% Ptotal = 25000;               % Total applied horizontal force
% nSteps = 10;
% maxIter = 25;                 % Maximum Newton corrections PER step
% tolResidual = 1e-10;
% tolVerification = 1e-8;
% 
% % Newton Command Window iteration log (true / false)
% showNewtonIterations = true;
% nonlinearMethod = 'full'; % 'full', 'modified', 'secant', 'numerical'
% useLoadBisection = false; % true: halve failed load increments
% maxCutbacks = 8;
% 
% 
% % Plot either solution independently (or both in separate figures).
% plotLinearContours = true;
% plotNonlinearContours = true;
% % Contour styling applies to any field (displacement, stress, strain).
% contourPalette = 'jet';
% plotOptions.showMax = true;   % Show maximum
% plotOptions.showNodes = false;   % Show node numbers
% % plotOptions.showNodes = false; % Hide node numbers
% % plotOptions.showMax = false; % Hide maximum
% contourScale = 1;
% 
% % Choose ONE boundary-condition case:
% %   'freeY'  : transverse contraction allowed; bottom-left uy fixed
% %   'fixedY' : all uy fixed; no transverse contraction
% boundaryCase = 'freeY';
% %materialModel = 'StVK';  % 'StVK' (professor) or 'LogElastic' (Abaqus-like)
% materialModel = 'LogElastic';
% X = [0 0; L 0; L H; 0 H];
% switch lower(materialModel)
%     case 'stvk'
%         material = fea.materials.LinearElastic(Eyoung,nu);
%     case 'logelastic'
%         material = fea.materials.LogElastic(Eyoung,nu);
%     otherwise
%         error('HW3:Material','Use materialModel = ''StVK'' or ''LogElastic''.');
% end
% element = fea.elements.Quad4(X,material,thickness);
% % Total Lagrangian response is now provided by the element directly.
% 
% switch lower(boundaryCase)
%     case 'freey'
%         fixedDOF = [1 2 7];         % Left edge ux=0; node 1 uy=0
%         caseDescription = 'Transverse contraction allowed';
%     case 'fixedy'
%         fixedDOF = [1 2 4 6 7 8]; % Left edge ux=0; all uy=0
%         caseDescription = 'Transverse contraction constrained';
%     otherwise
%         error('HW3:BC','Select boundaryCase = ''freeY'' or ''fixedY''.');
% end
% freeDOF  = setdiff(1:8,fixedDOF);
% rightDOF = [3 5];                 % Horizontal loading at nodes 2 and 3
% 
% loadFactor = (0:nSteps)'/nSteps;
% Psteps = Ptotal * loadFactor;
% q = zeros(8,1);
% qHistory = zeros(8,nSteps+1);
% uNL = zeros(nSteps+1,1);
% iters = zeros(nSteps,1);
% residualHistory = cell(nSteps,1);
% 
% %% 2. NEWTON-RAPHSON: REUSABLE OOP SOLVER
% % Keep this verification script simple; Q4 kinematics stay in TL.
% loadHistory = zeros(8,nSteps+1);
% loadHistory(rightDOF,:) = repmat(Psteps'/numel(rightDOF),numel(rightDOF),1);
% response = @(qTrial) q4Response(TL,element,material,qTrial);
% opts = struct('tolResidual',tolResidual,'maxIter',maxIter, ...
%     'verbose',showNewtonIterations,'printEvery',1, ...
%     'method',nonlinearMethod,'cutback',useLoadBisection, ...
%     'maxCutbacks',maxCutbacks);
% NL = fea.solvers.NewtonRaphson.solve(response,loadHistory,fixedDOF,opts);
% qHistory = NL.qHistory;
% q = NL.q;
% iters = NL.corrections;
% residualHistory = NL.residualHistory;
% uNL = mean(qHistory(rightDOF,:),1)';
% 
% fprintf('\nSingle Q4: %s | %s | P = %.6g\n',materialModel,caseDescription,Ptotal);
% fprintf('Step  Factor  Corrections    ux(right)    Rel. residual\n');
% for step=1:nSteps
%     rr=residualHistory{step};
%     fprintf('%3d  %6.3f     %3d        %11.8f      %.3e\n', ...
%         step,loadFactor(step+1),iters(step),uNL(step+1),rr(end));
% end
% 
% %% 3. LINEAR SOLUTION: EXISTING HW2 Q4 STIFFNESS
% Klinear = element.Gauss_Element_stiffness_cal();
% LIN = fea.solvers.LinearStatic.solve(Klinear,loadHistory,fixedDOF);
% qLinearHistory=LIN.qHistory;
% uLinear=mean(qLinearHistory(rightDOF,:),1)';
% 
% %% 4. HOMOGENEOUS ANALYTICAL REFERENCE FOR SELECTED MATERIAL
% % Both laws use the same plane-strain boundary conditions and reference
% % traction P11 = J*sigma11/lambdaX = lambdaY*sigma11 (LogElastic),
% % or P11 = lambdaX*S11 (StVK).
% lambdaL=material.lambda_L;
% mu=material.mu;
% C11=lambdaL+2*mu;
% if strcmpi(materialModel,'StVK')
%     E11fun=@(lx) .5*(lx.^2-1);
%     if strcmpi(boundaryCase,'freeY')
%         E22fun=@(lx) -(lambdaL/C11)*E11fun(lx);
%         lambdaYfun=@(lx) sqrt(1+2*E22fun(lx));
%     else
%         E22fun=@(lx) zeros(size(lx));
%         lambdaYfun=@(lx) ones(size(lx));
%     end
%     Pfun=@(lx) lx.*(C11*E11fun(lx)+lambdaL*E22fun(lx))*H*thickness;
% else
%     % sigma = lambda*tr(logV)*I+2*mu*logV, with epsilon_zz=0.
%     if strcmpi(boundaryCase,'freeY')
%         ratio=lambdaL/C11;
%         lambdaYfun=@(lx) lx.^(-ratio);     % sigma22=0
%         effectiveModulus=C11-lambdaL^2/C11;
%     else
%         lambdaYfun=@(lx) ones(size(lx));
%         effectiveModulus=C11;
%     end
%     Pfun=@(lx) lambdaYfun(lx).*effectiveModulus.*log(lx)*H*thickness;
% end
% 
% % StVK and log-elastic freeY may develop a turning point in nominal P(lx).
% % Find the FIRST positive-stretch root by scanning, rather than assuming
% % the force remains monotonic out to an arbitrary upper stretch.
% lambdaGrid=linspace(1,5,4001);
% Pgrid=Pfun(lambdaGrid);
% uExactSteps=zeros(nSteps+1,1);
% vExactSteps=zeros(nSteps+1,1);
% for k=2:nSteps+1
%     crossing=find(Pgrid>=Psteps(k),1,'first');
%     if isempty(crossing) || crossing<2 || ~isfinite(Pgrid(crossing))
%         error('HW3:AnalyticalRange', ...
%             'No first tensile equilibrium root at load %.6g; check limits.',Psteps(k));
%     end
%     lambdaXk=fzero(@(lx)Pfun(lx)-Psteps(k), ...
%         [lambdaGrid(crossing-1),lambdaGrid(crossing)]);
%     uExactSteps(k)=(lambdaXk-1)*L;
%     vExactSteps(k)=(lambdaYfun(lambdaXk)-1)*H;
% end
% uExact=uExactSteps(end);
% vExact=vExactSteps(end);
% 
% %% 5. FINAL REACTION FORCES AND VERIFICATION
% [fint,KL,KNL] = TL.Element_Response_cal(element,material,q);
% Fext = zeros(8,1);
% Fext(rightDOF) = Ptotal/numel(rightDOF);
% R = Fext-fint;
% reactions = fint-Fext;
% rightForceComputed = sum(fint(rightDOF));
% leftReactionX = sum(reactions([1 7])); % Signed, negative in tension
% uNewton = mean(q(rightDOF));
% vNewton = q(6)-q(4);               % Change in vertical edge length
% relAxialError = abs(uNewton-uExact)/max(1,abs(uExact));
% relTransverseError = abs(vNewton-vExact)/max(1,abs(vExact));
% relForceError = abs(rightForceComputed-Ptotal)/max(1,abs(Ptotal));
% relEquilibriumError = norm(R(freeDOF))/max(1,norm(Fext(freeDOF)));
% rightSymmetryError = abs(q(3)-q(5));
% 
% fprintf('\n--- FINAL FORCE-CONTROLLED COMPARISON ---\n');
% fprintf('Applied total force                    : %.12f\n',Ptotal);
% fprintf('Newton right-edge ux                   : %.12f\n',uNewton);
% fprintf('Analytical right-edge ux               : %.12f\n',uExact);
% fprintf('HW2 linear right-edge ux               : %.12f\n',uLinear(end));
% fprintf('Newton vertical change (uy3-uy2)      : %.12f\n',vNewton);
% fprintf('Analytical vertical change             : %.12f\n',vExact);
% fprintf('Computed right-edge internal force     : %.12f\n',rightForceComputed);
% fprintf('Signed left-edge X reaction            : %.12f\n',leftReactionX);
% fprintf('Axial displacement error               : %.6e\n',relAxialError);
% fprintf('Transverse displacement error          : %.6e\n',relTransverseError);
% fprintf('Force error                            : %.6e\n',relForceError);
% fprintf('Free-DOF equilibrium residual          : %.6e\n',relEquilibriumError);
% fprintf('Right-edge ux symmetry error           : %.6e\n',rightSymmetryError);
% fprintf('Linear vs nonlinear ux difference      : %.4f %%\n', ...
%     100*abs(uLinear(end)-uNewton)/max(eps,abs(uNewton)));
% 
% pass = relAxialError<tolVerification && ...
%        relTransverseError<tolVerification && ...
%        relForceError<tolVerification && ...
%        relEquilibriumError<tolResidual && ...
%        rightSymmetryError<tolVerification;
% if pass
%     fprintf('PASS - %s Newton result matches the chosen analytical reference.\n',materialModel);
% else
%     warning('HW3:VerificationFailed', ...
%         'One or more verification checks failed; inspect the reported errors.');
% end
% 
% %% 6. PLOT: FORCE-DISPLACEMENT (ONE FIGURE, THREE CURVES)
% figure('Name','Single Q4 - Force displacement','Color','w');
% hA = plot(uExactSteps,Psteps,'k-','LineWidth',1.6); hold on;
% hN = plot(uNL,Psteps,'bo-','LineWidth',1.2,'MarkerSize',5);
% hL = plot(uLinear,Psteps,'r--s','LineWidth',1.2,'MarkerSize',4);
% xlabel('Average right-edge horizontal displacement, u_x');
% ylabel('Total applied horizontal force, P');
% title(['Single Q4: ',caseDescription],'FontWeight','normal');
% legend([hA hN hL], ...
%     {'Selected analytical law','TL Newton-Raphson','HW2 linear'}, ...
%     'Location','best');
% grid on; box on;
% 
% %% 7. PLOT: NEWTON RESIDUALS (ONE FIGURE)
% figure('Name','Single Q4 - Newton convergence','Color','w');
% hold on;
% for step = 1:nSteps
%     rr = residualHistory{step};
%     semilogy(0:numel(rr)-1,max(rr,eps),'o-', ...
%         'LineWidth',1,'MarkerSize',3, ...
%         'DisplayName',sprintf('Step %d',step));
% end
% xlabel('Newton correction number');
% ylabel('Relative equilibrium residual');
% title('Newton convergence by load step','FontWeight','normal');
% legend('Location','bestoutside'); grid on; box on;
% 
% %% 8. PLOT: REFERENCE, LINEAR AND NONLINEAR DEFORMED SHAPES
% qXY = reshape(q,2,4).';
% Xdef = X + qXY;
% qLinearXY = reshape(qLinearHistory(:,end),2,4).';
% Xlin = X + qLinearXY;
% nodeLoop = [1 2 3 4 1];
% 
% figure('Name','Single Q4 - Shape comparison','Color','w');
% hRef = plot(X(nodeLoop,1),X(nodeLoop,2),'k--o','LineWidth',1.2); hold on;
% hDef = plot(Xdef(nodeLoop,1),Xdef(nodeLoop,2),'b-o','LineWidth',1.5);
% hLin = plot(Xlin(nodeLoop,1),Xlin(nodeLoop,2),'r--s','LineWidth',1.2);
% for node = 1:4
%     text(Xdef(node,1)+0.02*L,Xdef(node,2)+0.02*H, ...
%         sprintf('%d',node),'Color','b','FontSize',10);
% end
% axis equal; grid on; box on; xlabel('X'); ylabel('Y');
% title('Reference, linear and nonlinear shapes','FontWeight','normal');
% legend([hRef hDef hLin],{'Reference','TL Newton','HW2 linear'}, ...
%     'Location','best');
% 
% %% 8B. FIELD CONTOURS - EACH SOLUTION PLOTTED INDEPENDENTLY
% % One generic plotter handles displacement, strain, stress, etc.
% % Change the true/false flags below to select which result to display.
% fieldsLinear = [qLinearXY(:,1), qLinearXY(:,2), ...
%     hypot(qLinearXY(:,1),qLinearXY(:,2))];
% fieldsNonlinear = [qXY(:,1), qXY(:,2), ...
%     hypot(qXY(:,1),qXY(:,2))];
% 
% contourOptions = struct( ...
%     'location','nodal', ...
%     'fieldNames',{{'u_x','u_y','|u|'}}, ...
%     'colormap',contourPalette, ...
%     'scale',contourScale, ...
%     'showNodes',false);
% 
% if plotLinearContours
%     contourOptions.q = qLinearHistory(:,end);
%     contourOptions.figureName = 'Q4 - Linear Displacement Fields';
%     fea.results.FieldPlotter.plot( ...
%         X,[1 2 3 4],fieldsLinear,contourOptions);
% end
% 
% if plotNonlinearContours
%     contourOptions.q = q;
%     contourOptions.figureName = 'Q4 - Nonlinear Displacement Fields';
%     fea.results.FieldPlotter.plot( ...
%         X,[1 2 3 4],fieldsNonlinear,contourOptions);
% end
% 
% %% 9. INDEPENDENT NUMERICAL TANGENT CHECK AT CONVERGED DEFORMATION
% % Central difference of the ENTIRE element internal force vs the analytical
% % total element tangent. This validates stress/tangent integration together.
% [fcheck,Kmat,Kgeo]=TL.Element_Response_cal(element,material,q);
% Ktan=Kmat+Kgeo;
% Kfd=zeros(8,8);
% hFD=1e-6;
% for j=1:8
%     e=zeros(8,1); e(j)=hFD;
%     fp=TL.Element_Response_cal(element,material,q+e);
%     fm=TL.Element_Response_cal(element,material,q-e);
%     Kfd(:,j)=(fp-fm)/(2*hFD);
% end
% errorKT=norm(Ktan-Kfd,'fro')/max(1,norm(Kfd,'fro'));
% fprintf('Material %-10s: element tangent FD relative error = %.6e\n', ...
%     materialModel,errorKT);
% if errorKT>1e-5
%     warning('HW3:TangentCheck','Tangent FD check exceeded 1e-5.');
% end
% 
% %% Local adapter: future global assembly will replace only this callback.
% function [fint,Kt] = q4Response(TL,element,material,q)
%     [fint,KL,KNL] = TL.Element_Response_cal(element,material,q);
%     Kt=KL+KNL;
% end
