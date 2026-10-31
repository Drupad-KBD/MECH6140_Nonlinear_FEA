%% ================================================================
% HW3 - CANTILEVER: ABAQUS MESH + MATLAB FEA
% MECH 6140
% ================================================================

clearvars;
clc;
close all;

%% 1. USER INPUTS
tic
inpFile = 'Job-1.inp';

% Material
Eyoung = 1;
nu = 0.3;
thickness = 1;

materialModel = 'StVK';
%materialModel = 'LogElastic';

% Applied force
Ptotal = -0.01;       % Negative = downward

% Load stepping
nSteps = 30;

% Solver
nonlinearMethod = 'full';
%nonlinearMethod = 'modified';
%nonlinearMethod = 'secant';
%nonlinearMethod = 'numerical';

%bcMethod = 'elimination';
bcMethod = 'penalty';

penaltyFactor = 1e6;

maxIter = 30;
tolResidual = 1e-9;
showIterations = true;

useLoadBisection = true;

% Plot settings
plotLinearContours = true;
plotNonlinearContours = true;
showNodes = false;
contourPalette = 'jet';
contourScale = 1;
showAxisLabels = false;

%% 2. IMPORT ABAQUS MESH

mesh = fea.mesh.AbaqusINP.read(inpFile);

X = mesh.X;
conn = mesh.connectivity;

nNodes = size(X,1);
nElements = size(conn,1);
nDOF = 2*nNodes;

fprintf('\n--- IMPORTED MESH ---\n');
fprintf('Nodes          : %d\n',nNodes);
fprintf('Elements       : %d\n',nElements);
fprintf('Global DOFs    : %d\n',nDOF);
fprintf('Element type   : %s\n',mesh.elementType);

if endsWith(mesh.elementType,'R')
    warning(['Abaqus uses reduced integration, ' ...
             'but MATLAB uses full 2x2 integration.']);
end

%% 3. DEFINE BOUNDARY CONDITIONS IN MATLAB

tol = 1e-8*max([ ...
    max(X(:,1))-min(X(:,1)), ...
    max(X(:,2))-min(X(:,2)), ...
    1]);

% Left-edge nodes
leftNodes = find( ...
    abs(X(:,1)-min(X(:,1))) < tol);

% Fix UX and UY on left edge
fixedDOF = sort([ ...
    2*leftNodes-1; ...
    2*leftNodes]);
 
%% 4. DEFINE LOADING IN MATLAB

%loadType = 'traction';
loadType = 'point';

% Identify right-edge nodes
rightNodes = find( ...
    abs(X(:,1)-max(X(:,1))) < tol);

% Upper-right corner node
[~,index] = max(X(rightNodes,2));
loadNode = rightNodes(index);
loadDOF = 2*loadNode;

Ftotal = zeros(nDOF,1);

switch lower(loadType)

    case 'point'

        % Original downward point load
        Ftotal(loadDOF) = Ptotal;

    case 'traction'

        % Same total resultant force, distributed over right edge
        edgeLength = max(X(rightNodes,2)) ...
                   - min(X(rightNodes,2));

        traction = [0, Ptotal/(edgeLength*thickness)];

        % Q4 right boundary is local edge 2: nodes 2-3
        % This also handles differently oriented Q4 connectivity.

        edgeNodes = [1 2; 2 3; 3 4; 4 1];

        for e = 1:nElements

            ids = conn(e,:);

            for edge = 1:4

                localNodes = edgeNodes(edge,:);
                globalNodes = ids(localNodes);

                if all(ismember(globalNodes,rightNodes))

                    element = fea.elements.Quad4( ...
                        X(ids,:), ...
                        fea.materials.LinearElastic(Eyoung,nu), ...
                        thickness);

                    fe = fea.loads.SurfaceTraction.EdgeForce_cal( ...
                        element,edge,traction);

                    dofs = reshape( ...
                        [2*ids-1; 2*ids],[],1);

                    Ftotal(dofs) = Ftotal(dofs)+fe;
                end
            end
        end

    otherwise
        error('Unknown load type.');

end

% Load history
loadFactors = linspace(0,1,nSteps+1);
loadHistory = Ftotal*loadFactors;

fprintf('\n--- LOADING ---\n');
fprintf('Load type      : %s\n',loadType);
fprintf('Total applied FY: %.6g\n',sum(Ftotal(2:2:end)));


%% 5. MATERIAL AND FORMULATION

switch lower(materialModel)

    case 'stvk'
        material = ...
            fea.materials.LinearElastic(Eyoung,nu);

    case 'logelastic'
        material = ...
            fea.materials.LogElastic(Eyoung,nu);

    otherwise
        error('Unknown material model.');
end

% Total Lagrangian element calculations are in fea.elements.ElementCalculations

%% 6. GLOBAL LINEAR STIFFNESS

fprintf('\nAssembling global stiffness...\n');

q0 = zeros(nDOF,1);

[~,KL0,KNL0] = ...
    fea.assembly.Assembler.Assemble_TL_cal( ...
    X,conn,material,thickness,q0);

Klinear = KL0 + KNL0;

%% 7. LINEAR SOLUTION

linearOpts = struct( ...
    'bcMethod',bcMethod, ...
    'penaltyFactor',penaltyFactor);

LIN = fea.solvers.LinearStatic.solve( ...
    Klinear,loadHistory,fixedDOF,linearOpts);

qLinear = LIN.q;

%% 8. NONLINEAR SOLUTION

response = @(q) globalResponse( ...
    X,conn,material,thickness,q);

nonlinearOpts = struct( ...
    'method',nonlinearMethod, ...
    'maxIter',maxIter, ...
    'tolResidual',tolResidual, ...
    'verbose',showIterations, ...
    'printEvery',1, ...
    'cutback',useLoadBisection, ...
    'maxCutbacks',8, ...
    'bcMethod',bcMethod, ...
    'penaltyFactor',penaltyFactor);

NL = fea.solvers.NewtonRaphson.solve( ...
    response,loadHistory,fixedDOF,nonlinearOpts);

qNonlinear = NL.q;
toc
%% 9. PRINT RESULTS

tipUxDOF = 2*loadNode-1;
tipUyDOF = 2*loadNode;

fprintf('\n--- FINAL TIP DISPLACEMENTS ---\n');

fprintf('                UX              UY\n');

fprintf('Linear     % .6e    % .6e\n', ...
    qLinear(tipUxDOF),qLinear(tipUyDOF));

fprintf('Nonlinear  % .6e    % .6e\n', ...
    qNonlinear(tipUxDOF),qNonlinear(tipUyDOF));

%% 10. FORCE-DISPLACEMENT

figure('Name','Tip Force-Displacement','Color','w');

plot( ...
    LIN.qHistory(tipUyDOF,:), ...
    Ptotal*loadFactors, ...
    'r-','LineWidth',1.3);
hold on;

plot( ...
    NL.qHistory(tipUyDOF,:), ...
    Ptotal*loadFactors, ...
    'b-','LineWidth',1.3);

xlabel('Tip displacement U_Y');
ylabel('Applied force F_Y');

legend('Linear','Nonlinear TL','Location','best');

title('Cantilever Force-Displacement', ...
    'FontWeight','normal');

grid on;
box on;

%% 11. DEFORMED MESH COMPARISON

uL = reshape(qLinear,2,[])';
uNL = reshape(qNonlinear,2,[])';

Xlinear = X + contourScale*uL;
Xnonlinear = X + contourScale*uNL;

figure('Name','Deformed Mesh Comparison','Color','w');
hold on;

drawMesh(X,conn,[0.4 0.4 0.4],':');
drawMesh(Xlinear,conn,[0.85 0.1 0.1],'--');
drawMesh(Xnonlinear,conn,[0.1 0.3 0.85],'-');

axis equal;
grid on;
box on;

xlabel('X');
ylabel('Y');
title('Cantilever Deformation','FontWeight','normal');

%% 12. DISPLACEMENT FIELDS

fieldsLinear = [ ...
    uL(:,1), ...
    uL(:,2), ...
    hypot(uL(:,1),uL(:,2))];

fieldsNonlinear = [ ...
    uNL(:,1), ...
    uNL(:,2), ...
    hypot(uNL(:,1),uNL(:,2))];

plotOptions = struct();

plotOptions.location = 'nodal';
plotOptions.colormap = contourPalette;
plotOptions.scale = contourScale;
plotOptions.showNodes = showNodes;
plotOptions.showAxisLabels = showAxisLabels;
plotOptions.fieldNames = {'u_x','u_y','|u|'};

% Linear contours
if plotLinearContours

    plotOptions.q = qLinear;
    plotOptions.figureName = 'Linear Displacement Fields';

    fea.results.FieldPlotter.plot( ...
        X,conn,fieldsLinear,plotOptions);
end

% Nonlinear contours
if plotNonlinearContours

    plotOptions.q = qNonlinear;
    plotOptions.figureName = 'Nonlinear Displacement Fields';

    fea.results.FieldPlotter.plot( ...
        X,conn,fieldsNonlinear,plotOptions);
end

%% 13. NEWTON CONVERGENCE

figure('Name','Newton Convergence','Color','w');
hold on;

for step = 1:nSteps

    rr = NL.residualHistory{step};

    semilogy( ...
        0:numel(rr)-1, ...
        max(rr,eps), ...
        'o-', ...
        'LineWidth',1.1, ...
        'DisplayName',sprintf('Load Step %d',step));
end

xlabel('Newton Iteration');
ylabel('Relative Residual');

title('Newton Convergence','FontWeight','normal');

legend('Location','bestoutside');

grid on;
box on;

if any(NL.cutbacks > 0)
    warning(['Cutback histories contain multiple Newton attempts ' ...
             'and should be plotted separately.']);
end


%% 14. SAVE RESULTS

results.mesh.X = X;
results.mesh.connectivity = conn;

results.material.E = Eyoung;
results.material.nu = nu;
results.material.thickness = thickness;
results.material.model = materialModel;

results.loading.type = loadType;
results.loading.Ftotal = Ftotal;
results.loading.loadFactors = loadFactors;
results.loading.loadHistory = loadHistory;

results.boundary.fixedDOF = fixedDOF;
results.boundary.method = bcMethod;

results.linear = LIN;
results.nonlinear = NL;

% Displacement history [node x direction x load step]
results.displacements.linear = reshape( ...
    LIN.qHistory,2,nNodes,[]);

results.displacements.nonlinear = reshape( ...
    NL.qHistory,2,nNodes,[]);

% Tip displacement history
results.tip.UxLinear = LIN.qHistory(tipUxDOF,:);
results.tip.UyLinear = LIN.qHistory(tipUyDOF,:);

results.tip.UxNonlinear = NL.qHistory(tipUxDOF,:);
results.tip.UyNonlinear = NL.qHistory(tipUyDOF,:);

% Save MATLAB file
save('FEA_Results.mat','results');

fprintf('\nResults saved to FEA_Results.mat\n');
%% ================================================================
% LOCAL FUNCTIONS
% ================================================================

function [fint,Kt] = globalResponse( ...
    X,conn,material,thickness,q)

    [fint,~,~,Kt] = ...
        fea.assembly.Assembler.Assemble_TL_cal( ...
        X,conn,material,thickness,q);
end

function drawMesh(X,conn,color,lineStyle)

    for e = 1:size(conn,1)

        ids = conn(e,[1 2 3 4 1]);

        plot(X(ids,1),X(ids,2), ...
            'Color',color, ...
            'LineStyle',lineStyle, ...
            'LineWidth',0.6);
    end
end























































% %% ================================================================
% % HW3 - VERIFY ABAQUS MESH IMPORT
% % ================================================================
% 
% clearvars;
% clc;
% close all;
% 
% %% 1. READ ABAQUS MESH
% 
% inpFile = 'Job-1.inp';
% 
% mesh = fea.mesh.AbaqusINP.read(inpFile);
% 
% X = mesh.X;
% conn = mesh.connectivity;
% 
% %% 2. PRINT MESH INFORMATION
% 
% fprintf('\n========== ABAQUS MESH IMPORT ==========\n');
% fprintf('Nodes        : %d\n',mesh.nNodes);
% fprintf('Elements     : %d\n',mesh.nElements);
% fprintf('Element type : %s\n',mesh.elementType);
% fprintf('Total DOFs   : %d\n',2*mesh.nNodes);
% 
% fprintf('\nX range: %.4f to %.4f\n',min(X(:,1)),max(X(:,1)));
% fprintf('Y range: %.4f to %.4f\n',min(X(:,2)),max(X(:,2)));
% 
% %% 3. PRINT NODE SETS
% 
% fprintf('\n========== NODE SETS ==========\n');
% 
% for k = 1:numel(mesh.nsets)
% 
%     fprintf('%-20s | %-12s | %d nodes\n', ...
%         mesh.nsets(k).name, ...
%         mesh.nsets(k).scope, ...
%         numel(mesh.nsets(k).nodes));
% end
% 
% %% 4. PRINT ELEMENT SETS
% 
% fprintf('\n========== ELEMENT SETS ==========\n');
% 
% for k = 1:numel(mesh.elsets)
% 
%     fprintf('%-20s | %-12s | %d elements\n', ...
%         mesh.elsets(k).name, ...
%         mesh.elsets(k).scope, ...
%         numel(mesh.elsets(k).elements));
% end
% 
% %% 5. PLOT IMPORTED MESH
% 
% figure('Name','Abaqus Imported Mesh','Color','w');
% 
% patch('Faces',conn, ...
%       'Vertices',X, ...
%       'FaceColor','none', ...
%       'EdgeColor','k', ...
%       'LineWidth',0.7);
% 
% axis equal;
% axis tight;
% grid on;
% box on;
% 
% xlabel('X');
% ylabel('Y');
% title('Imported Abaqus Q4 Mesh','FontWeight','normal');
% 
% %% 6. VERIFY CONNECTIVITY
% 
% assert(all(conn(:)>=1) && ...
%        all(conn(:)<=mesh.nNodes), ...
%        'Invalid mesh connectivity.');
% 
% fprintf('\nPASS: Mesh imported successfully.\n');