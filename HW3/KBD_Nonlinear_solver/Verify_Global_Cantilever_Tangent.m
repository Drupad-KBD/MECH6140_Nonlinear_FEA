%% ================================================================
% VERIFY: FULL CANTILEVER GLOBAL TANGENT (Q4 OR Q8)
% Uses a saved NONZERO displacement state from FEA_Results.mat.
% Compares the assembled physical tangent with a centered finite
% difference of the assembled GLOBAL internal-force vector.
% No boundary-condition penalty is applied to either matrix.
% ================================================================
clearvars; clc; close all;

%% SETTINGS
resultsFile = 'FEA_Results.mat';
stateIndex = [];            % [] = final converged saved state
hRelative = 1e-6;           % h = hRelative * characteristic mesh length
saveFigures = true;
figureFolder = 'FEA_Figures';

%% LOAD SAVED MODEL AND GLOBAL DISPLACEMENTS
s = load(resultsFile,'results');
r = s.results;
X = r.mesh.X;
conn = r.mesh.connectivity;
qHistory = r.nonlinear.qHistory;

nNodes = size(X,1);
nDOF = 2*nNodes;
assert(size(qHistory,1)==nDOF,'Saved displacement DOF count does not match mesh.');

if isempty(stateIndex)
    stateIndex = size(qHistory,2);
end
assert(isscalar(stateIndex) && stateIndex>=2 && ...
       stateIndex<=size(qHistory,2) && stateIndex==floor(stateIndex), ...
       'Choose a valid nonzero load-state index.');

q = qHistory(:,stateIndex);
assert(norm(q)>0,'The selected displacement state is zero.');

switch lower(r.material.model)
    case 'stvk'
        material = fea.materials.LinearElastic(r.material.E,r.material.nu);
    case 'logelastic'
        material = fea.materials.LogElastic(r.material.E,r.material.nu);
    otherwise
        error('Unrecognized saved material: %s',r.material.model);
end
thickness = r.material.thickness;

% Global perturbation step, in the same length units as mesh coordinates.
Lchar = max(max(X,[],1)-min(X,[],1));
h = hRelative*Lchar;
assert(h>0 && isfinite(h),'Invalid perturbation size.');

fprintf('\nFULL GLOBAL TANGENT CHECK\n');
fprintf('Nodes: %d | Elements: %d | Global DOFs: %d\n', ...
    nNodes,size(conn,1),nDOF);
fprintf('Saved state index: %d | ||q||_2 = %.8e\n',stateIndex,norm(q));
fprintf('Central difference perturbation h = %.8e\n',h);

%% 1. ANALYTICAL ASSEMBLED GLOBAL TANGENT
% Identical geometry, material, and q for both comparisons.
[fint,KL,KNL,Kanalytical] = ...
    fea.assembly.Assembler.Assemble_TL_cal(X,conn,material,thickness,q);

%% 2. NUMERICAL ASSEMBLED GLOBAL TANGENT
% Column j = [fint(q+h*ej)-fint(q-h*ej)]/(2*h)
% Perturb EVERY global displacement DOF, including constrained DOFs.
% Constrained DOFs are still physical kinematic variables for this check.
Knumerical = zeros(nDOF,nDOF);
for j = 1:nDOF
    qPlus = q;
    qMinus = q;
    qPlus(j) = qPlus(j)+h;
    qMinus(j) = qMinus(j)-h;

    fPlus = fea.assembly.Assembler.Assemble_TL_cal( ...
        X,conn,material,thickness,qPlus);
    fMinus = fea.assembly.Assembler.Assemble_TL_cal( ...
        X,conn,material,thickness,qMinus);

    Knumerical(:,j) = (fPlus-fMinus)/(2*h);

    if mod(j,25)==0 || j==nDOF
        fprintf('Checked %d / %d global displacement DOFs\n',j,nDOF);
    end
end

%% 3. ERROR NORMS AND HEATMAP
D = full(Kanalytical)-Knumerical;
Kref = full(Kanalytical);

relativeFrobenius = norm(D,'fro')/max(norm(Kref,'fro'),eps);
absoluteFrobenius = norm(D,'fro');
maximumAbsoluteError = max(abs(D(:)));
maximumStiffness = max(abs(Kref(:)));
maximumNormalizedEntryError = maximumAbsoluteError/max(maximumStiffness,eps);

% Compare physical tangent on unconstrained DOFs as an additional check.
fixedDOF = unique(r.boundary.fixedDOF(:));
freeDOF = setdiff((1:nDOF)',fixedDOF);
Dfree = D(freeDOF,freeDOF);
Kfree = Kref(freeDOF,freeDOF);
relativeFreeError = norm(Dfree,'fro')/max(norm(Kfree,'fro'),eps);

% A single global scale avoids blow-ups from structural zero entries.
relativeMap = abs(D)/max(maximumStiffness,eps);

fprintf('\n--- GLOBAL TANGENT VERIFICATION ---\n');
fprintf('||Kanalytical - Knumerical||_F      = %.10e\n',absoluteFrobenius);
fprintf('Relative GLOBAL Frobenius error    = %.10e\n',relativeFrobenius);
fprintf('Relative FREE-DOF Frobenius error  = %.10e\n',relativeFreeError);
fprintf('Maximum absolute entry error      = %.10e\n',maximumAbsoluteError);
fprintf('Maximum normalized entry error    = %.10e\n',maximumNormalizedEntryError);

%% 4. FULL GLOBAL MATRIX HEATMAP
figure('Name','Full Global Tangent Relative Difference', ...
       'Color','w','Units','pixels','Position',[100 100 900 800]);
imagesc(relativeMap);
axis image;
colormap(parula);
cb = colorbar;
cb.Label.String = '|K_{analytical}-K_{numerical}| / max|K_{analytical}|';
xlabel('Global displacement DOF j');
ylabel('Global internal-force DOF i');
title(sprintf('Full Cantilever Tangent Error, ||Delta K||_F / ||K||_F = %.3e', ...
      relativeFrobenius),'FontWeight','normal');
set(gca,'YDir','normal');

%% OPTIONAL: DISPLAY MATRIX VALUES SIDE BY SIDE
% Note: the heatmap is of the entire nDOF-by-nDOF system.
% Both matrices are stored for further inspection.
verification = struct();
verification.Kanalytical = Kanalytical;
verification.Knumerical = Knumerical;
verification.difference = D;
verification.relativeMap = relativeMap;
verification.relativeFrobenius = relativeFrobenius;
verification.relativeFreeError = relativeFreeError;
verification.maximumAbsoluteError = maximumAbsoluteError;
verification.maximumNormalizedEntryError = maximumNormalizedEntryError;
verification.stateIndex = stateIndex;
verification.perturbation = h;
verification.nDOF = nDOF;

if saveFigures
    if ~exist(figureFolder,'dir'), mkdir(figureFolder); end
    exportgraphics(gcf,fullfile(figureFolder,'Full_Global_Tangent_Heatmap.png'), ...
        'Resolution',300);
end
save('FEA_Global_Tangent_Verification.mat','verification','-v7.3');
fprintf('Saved FEA_Global_Tangent_Verification.mat\n');
