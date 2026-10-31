%% Cantilever_Cantilever_From_Abaqus.m
% Imported Abaqus Q4 mesh -> existing TL assembler -> existing OOP solvers.
% Run from the project root containing +fea_hw2 and the .inp file.
clearvars; clc; close all;

%% USER SETTINGS
inpFile = 'Job-1.inp';             % Abaqus mesh only
Eyoung = 1;                         % Defined in MATLAB
nu = 0.3;                           % Defined in MATLAB
Ptotal = -0.001;                    % Total force, negative is downward
materialModel = 'StVK';           % 'StVK' or 'LogElastic'
thickness = 1;                   % 2D reference thickness; confirm Abaqus section thickness
nSteps = 10;
bcMethod = 'elimination';        % 'elimination' or 'penalty'
penaltyFactor = 1e6;
nonlinearMethod = 'full';        % full, modified, secant, numerical
verbose = true;
useLoadBisection = false;
maxIter = 30;
tolResidual = 1e-9;
plotLinearContours = true;
plotNonlinearContours = true;
showNodes = false;
contourPalette = 'jet';
contourScale = 1;

%% 1. READ ABAQUS FILE
mesh = fea.mesh.AbaqusINP.read(inpFile);
X = mesh.X;
connectivity = mesh.connectivity;
nNodes = mesh.nNodes;
nDOF = 2*nNodes;
% Define supports and loading entirely in MATLAB.
span = max([max(X(:,1))-min(X(:,1)), max(X(:,2))-min(X(:,2)), 1]);
tol = 1e-8*span;
leftNodes = find(abs(X(:,1)-min(X(:,1)))<tol);
rightNodes = find(abs(X(:,1)-max(X(:,1)))<tol);
[~,topIdx] = max(X(rightNodes,2));
tipNode = rightNodes(topIdx);
fixedDOF = sort([2*leftNodes-1; 2*leftNodes]);
Ftotal = zeros(nDOF,1);
Ftotal(2*tipNode) = Ptotal;

fprintf('\nABAQUS INPUT: %s\n',inpFile);
fprintf('Element: %s | nodes: %d | Q4 elements: %d | DOFs: %d\n', ...
    mesh.elementType,nNodes,mesh.nElements,nDOF);
fprintf('X = [%g,%g], Y = [%g,%g]\n', ...
    min(X(:,1)),max(X(:,1)),min(X(:,2)),max(X(:,2)));
fprintf('E = %g, nu = %g | fixed DOFs = %d\n', ...
    Eyoung,nu,numel(fixedDOF));
fprintf('MATLAB applied force = [%g, %g]\n', ...
    sum(Ftotal(1:2:end)),sum(Ftotal(2:2:end)));

if endsWith(mesh.elementType,'R')
    warning('Cantilever:IntegrationMismatch', ...
       ['Abaqus file uses %s (reduced integration), but MATLAB ' ...
        'Q4 is full 2x2 Gauss integration. Use CPE4 in Abaqus ' ...
        'for a direct element comparison.'],mesh.elementType);
end
assert(~isempty(leftNodes), 'Could not identify left support nodes.');

%% 2. MODEL AND LOAD HISTORY
switch lower(materialModel)
    case 'stvk'
        material = fea.materials.LinearElastic(Eyoung,nu);
    case 'logelastic'
        material = fea.materials.LogElastic(Eyoung,nu);
    otherwise
        error('Invalid materialModel.');
end
TL = fea.formulations.TotalLagrangian();
loadFactors = linspace(0,1,nSteps+1);
loadHistory = Ftotal*loadFactors;

% Callback used by the existing Newton solver.
response = @(q) globalResponse(X,connectivity,material,thickness,q,TL);

%% 3. GLOBAL LINEAR STIFFNESS: INITIAL TL MATERIAL TANGENT
fprintf('\nAssembling initial global stiffness...\n');
[~,KL0,KNL0] = fea.assembly.Assembler.Assemble_TL_Q4_cal( ...
    X,connectivity,material,thickness,zeros(nDOF,1),TL);
Klinear = KL0+KNL0;

linOpts = struct('bcMethod',bcMethod,'penaltyFactor',penaltyFactor);
LIN = fea.solvers.LinearStatic.solve( ...
    Klinear,loadHistory,fixedDOF,linOpts);

%% 4. NONLINEAR SOLUTION
nlOpts = struct('method',nonlinearMethod,'maxIter',maxIter, ...
    'tolResidual',tolResidual,'verbose',verbose,'printEvery',1, ...
    'cutback',useLoadBisection,'maxCutbacks',8, ...
    'bcMethod',bcMethod,'penaltyFactor',penaltyFactor);
NL = fea.solvers.NewtonRaphson.solve( ...
    response,loadHistory,fixedDOF,nlOpts);

%% 5. REFERENCE TIP NODE AND RESULTS
assert(isscalar(tipNode),'Could not identify unique upper-right tip node.');
tipX=2*tipNode-1; tipY=2*tipNode;
qLin=LIN.q; qNL=NL.q;
fprintf('\nUPPER-RIGHT TIP (Abaqus node %d)\n',mesh.nodeLabels(tipNode));
fprintf('             u_x             u_y\n');
fprintf('Linear     % .9e   % .9e\n',qLin(tipX),qLin(tipY));
fprintf('Nonlinear  % .9e   % .9e\n',qNL(tipX),qNL(tipY));

%% 6. FORCE-DISPLACEMENT
% MATLAB applies a vertical corner load; plot vertical tip displacement.
figure('Name','Cantilever Cantilever - Tip Force-Displacement','Color','w');
plot(LIN.qHistory(tipY,:),loadFactors*Ftotal(tipY),'r--s', ...
    'LineWidth',1.3,'MarkerSize',4); hold on;
plot(NL.qHistory(tipY,:),loadFactors*Ftotal(tipY),'b-o', ...
    'LineWidth',1.3,'MarkerSize',4);
xlabel('Upper-right tip displacement u_y');
ylabel('Applied tip force F_y');
legend({'Linear','Nonlinear TL'},'Location','best');
title('Cantilever tip response','FontWeight','normal'); grid on; box on;

%% 7. DEFORMED SHAPE: ORDINARY MATLAB FIGURE
figure('Name','Cantilever Cantilever - Deformation','Color','w'); hold on;
plotQ4Edges(X,connectivity,[0.4 0.4 0.4],':');
plotQ4Edges(X + reshape(qLin,2,[])'*contourScale, ...
    connectivity,[0.8 0.1 0.1],'--');
plotQ4Edges(X + reshape(qNL,2,[])'*contourScale, ...
    connectivity,[0.1 0.3 0.9],'-');
axis equal; box on; grid on;
xlabel('X'); ylabel('Y');
title('Reference / Linear / Nonlinear mesh','FontWeight','normal');

%% 8. GENERAL FIELD CONTOURS (ONE SOLUTION PER CALL)
contourOptions = struct( ...
    'location','nodal','colormap',contourPalette, ...
    'scale',contourScale,'showNodes',showNodes, ...
    'fieldNames',{{'u_x','u_y','|u|'}});
if plotLinearContours
    uv = reshape(qLin,2,[])';
    vals = [uv(:,1),uv(:,2),hypot(uv(:,1),uv(:,2))];
    contourOptions.q = qLin;
    contourOptions.figureName = 'Cantilever - Linear Displacement';
    fea.results.FieldPlotter.plot(X,connectivity,vals,contourOptions);
end
if plotNonlinearContours
    uv = reshape(qNL,2,[])';
    vals = [uv(:,1),uv(:,2),hypot(uv(:,1),uv(:,2))];
    contourOptions.q = qNL;
    contourOptions.figureName = 'Cantilever - Nonlinear Displacement';
    fea.results.FieldPlotter.plot(X,connectivity,vals,contourOptions);
end

%% 9. CONVERGENCE HISTORY (NO CUTBACK-TRACES IF CUTBACKS USED)
figure('Name','Cantilever Cantilever - Newton Convergence','Color','w'); hold on;
for step=1:nSteps
    rr=NL.residualHistory{step};
    semilogy(0:numel(rr)-1,max(rr,eps),'o-', ...
        'LineWidth',1,'MarkerSize',3, ...
        'DisplayName',sprintf('Load step %d',step));
end
xlabel('Newton iteration'); ylabel('Relative equilibrium residual');
title('Newton convergence','FontWeight','normal');
legend('Location','bestoutside'); grid on; box on;
if any(NL.cutbacks>0)
    warning('Cantilever:CutbackPlot', ...
        'Residual plot concatenates substep attempts at steps with cutbacks.');
end

%% LOCAL FUNCTIONS
function [fint,Kt] = globalResponse(X,conn,material,t,q,TL)
    [fint,~,~,Kt] = fea.assembly.Assembler.Assemble_TL_Q4_cal( ...
        X,conn,material,t,q,TL);
end
function plotQ4Edges(X,conn,color,linestyle)
    % Using NaN breaks to plot each mesh edge in one line call.
    for e=1:size(conn,1)
        ids=conn(e,[1 2 3 4 1]);
        plot(X(ids,1),X(ids,2),'Color',color, ...
            'LineStyle',linestyle,'LineWidth',0.6);
    end
end
