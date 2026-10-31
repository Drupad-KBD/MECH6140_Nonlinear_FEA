%% Plot verification and experiment results from saved FEA data
% Q3.1: small-deformation limit (load X, tip displacement Y)
% Q4.1: load-displacement
% Q4.2: displacement percent difference
% Q4.3: three nonlinear deformed states
% No solver or +fea package is needed to run this script.
clearvars; clc; close all;

%% SETTINGS
resultsFile = 'FEA_Results.mat';
plotFolder = 'FEA_Figures';
saveFigures = true;
lineWidth = 2;
shapeScale = 1;                 % 1 = actual deformation
smallSteps = 3;                 % first three NONZERO saved load states
shapeSteps = [];                % [] = first, middle, last nonzero state

%% LOAD DATA
s = load(resultsFile,'results');
results = s.results;
X = results.mesh.X;
conn = results.mesh.connectivity;
qL = results.linear.qHistory;
qNL = results.nonlinear.qHistory;

if isfield(results.loading,'totalFy')
    P = results.loading.totalFy(:);
else
    P = sum(results.loading.loadHistory(2:2:end,:),1).';
end

if isfield(results.tip,'UyLinear') && isfield(results.tip,'UyNonlinear')
    uyL = results.tip.UyLinear(:);
    uyNL = results.tip.UyNonlinear(:);
else
    % Select the top-right node from the saved geometry.
    [~,tipNode] = max(X(:,1) + 1e-8*X(:,2));
    uyL = qL(2*tipNode,:).';
    uyNL = qNL(2*tipNode,:).';
end

% Use positive magnitudes for the downward-loading plots.
% Original signed P, uyL, and uyNL remain unchanged in the saved results.
loadMag = abs(P);
dispLinear = abs(uyL);
dispNonlinear = abs(uyNL);

nStates = numel(P);
assert(nStates==size(qL,2) && nStates==size(qNL,2), ...
    'Load and displacement history lengths do not agree.');
assert(nStates>=4,'At least three nonzero load states are needed for the shape plot.');
if saveFigures && ~exist(plotFolder,'dir'), mkdir(plotFolder); end

%% Q3.1 - SMALL DEFORMATION LIMIT
% These are the first available nonzero load steps, not necessarily
% sufficiently small. Inspect printed difference and rerun smaller loads
% if they do not closely approach the linear solution.
kSmall = 2:min(nStates,smallSteps+1);
figure('Name','Q3.1 Small Deformation','Color','w');
plot(loadMag(1:kSmall(end)),dispLinear(1:kSmall(end)),'r-','LineWidth',lineWidth); hold on;
plot(loadMag(1:kSmall(end)),dispNonlinear(1:kSmall(end)),'b-','LineWidth',lineWidth);
xlabel('Applied load magnitude |F_y|'); ylabel('Tip displacement magnitude |u_y|');
xlim([0,1*max(loadMag(kSmall))]);
ylim([0,1*max([dispLinear(kSmall);dispNonlinear(kSmall)])]);
title('Q3.1  Small-deformation limit','FontWeight','normal');
legend('Linear','Nonlinear TL','Location','best'); grid off; box on;
if saveFigures, exportgraphics(gcf,fullfile(plotFolder,'Q3_1_SmallDeformation.png'),'Resolution',300); end

smallDiff = 100*abs(uyNL(kSmall)-uyL(kSmall))./max(abs(uyL(kSmall)),eps);
fprintf('\nQ3.1: first %d nonzero states, max difference = %.6g %%\n', ...
    numel(kSmall),max(smallDiff));
if max(smallDiff)>1
    fprintf('NOTE: These states may not be in the small-deformation limit.\n');
    fprintf('      Run a separate analysis with smaller total applied force.\n');
end

%% Q4.1 - FULL LOAD-DISPLACEMENT
figure('Name','Q4.1 Load Displacement','Color','w');
plot(loadMag,dispLinear,'r-','LineWidth',lineWidth); hold on;
plot(loadMag,dispNonlinear,'b-','LineWidth',lineWidth);
xlabel('Applied load magnitude |F_y|'); ylabel('Tip displacement magnitude |u_y|');
xlim([0,1.05*max(loadMag)]);
ylim([0,1.05*max([dispLinear;dispNonlinear])]);
title('Q4.1  Load-displacement response','FontWeight','normal');
legend('Linear','Nonlinear TL','Location','best'); grid off; box on;
if saveFigures, exportgraphics(gcf,fullfile(plotFolder,'Q4_1_LoadDisplacement.png'),'Resolution',300); end

%% Q4.2 - PERCENT DIFFERENCE (LINEAR REFERENCE)
% 100 |uNL-uL|/|uL|; exclude zero load since denominator is zero.
percentDiff = NaN(nStates,1);
valid = abs(uyL)>max(1e-14,1e-12*max(abs(uyL)));
percentDiff(valid) = 100*abs(uyNL(valid)-uyL(valid))./abs(uyL(valid));
assert(any(valid),'No nonzero displacement values for percentage comparison.');
figure('Name','Q4.2 Percent Difference','Color','w');
plot(loadMag(valid),percentDiff(valid),'b-','LineWidth',lineWidth);
xlim([0,1.05*max(loadMag)]);
ylim([0,1.05*max(percentDiff(valid))]);
xlabel('Applied load magnitude |F_y|');
ylabel('Displacement difference (%)');
title('Q4.2  Nonlinear vs linear displacement','FontWeight','normal');
grid off; box on;
if saveFigures, exportgraphics(gcf,fullfile(plotFolder,'Q4_2_PercentDifference.png'),'Resolution',300); end


%% ================================================================
% 4.3 - DEFORMED SHAPES: LINEAR AND NONLINEAR
% ================================================================

load('FEA_Results.mat','results');

X = results.mesh.X;
conn = results.mesh.connectivity;

nNodes = size(X,1);
nStates = size(results.nonlinear.qHistory,2);

% Three load states
shapeSteps = [2, round((nStates+1)/2), nStates];

% Deformation scale
shapeScale = 1;

% Colors
originalColor = [0.55 0.55 0.55];

linearColors = [ ...
    1.00 0.65 0.65;
    0.85 0.30 0.30;
    0.55 0.00 0.00];

nonlinearColors = [ ...
    0.65 0.80 1.00;
    0.25 0.50 0.85;
    0.00 0.15 0.55];

originalWidth = 2.5;
shapeWidth = 2;

% Element boundary ordering
if size(conn,2) == 4
    edgeOrder = [1 2 3 4 1];
elseif size(conn,2) == 8
    edgeOrder = [1 5 2 6 3 7 4 8 1];
else
    error('Only Quad4 and Quad8 supported.');
end

%% Determine common axis limits for Linear and Nonlinear

allX = X(:,1);
allY = X(:,2);

for model = 1:2

    if model == 1
        qHistory = results.linear.qHistory;
    else
        qHistory = results.nonlinear.qHistory;
    end

    for k = 1:numel(shapeSteps)

        q = qHistory(:,shapeSteps(k));
        U = reshape(q,2,nNodes).';

        Xdef = X + shapeScale*U;

        allX = [allX; Xdef(:,1)];
        allY = [allY; Xdef(:,2)];
    end
end

xPadding = 0.05*max(max(allX)-min(allX),eps);
yPadding = 0.05*max(max(allY)-min(allY),eps);

commonXLim = [min(allX)-xPadding, max(allX)+xPadding];
commonYLim = [min(allY)-yPadding, max(allY)+yPadding];

%% Four figures: shared limits and individually fitted limits

for axisMode = 1:2

    for model = 1:2

        if model == 1
            qHistory = results.linear.qHistory;
            shapeColors = linearColors;
            modelName = 'Linear';
        else
            qHistory = results.nonlinear.qHistory;
            shapeColors = nonlinearColors;
            modelName = 'Nonlinear';
        end

        figure('Name',sprintf('%s Shapes - Mode %d', ...
            modelName,axisMode),'Color','w');
        hold on;

        % Original geometry
        for e = 1:size(conn,1)

            ids = conn(e,edgeOrder);

            plot(X(ids,1),X(ids,2), ...
                'Color',originalColor, ...
                'LineWidth',originalWidth, ...
                'HandleVisibility','off');
        end

        % Three deformed shapes
        for k = 1:numel(shapeSteps)

            q = qHistory(:,shapeSteps(k));
            U = reshape(q,2,nNodes).';

            Xdef = X + shapeScale*U;

            for e = 1:size(conn,1)

                ids = conn(e,edgeOrder);

                plot(Xdef(ids,1),Xdef(ids,2), ...
                    'Color',shapeColors(k,:), ...
                    'LineWidth',shapeWidth, ...
                    'HandleVisibility','off');
            end
        end

        % Legend
        plot(nan,nan,'-', ...
            'Color',originalColor, ...
            'LineWidth',originalWidth, ...
            'DisplayName','Undeformed');

        for k = 1:numel(shapeSteps)

            plot(nan,nan,'-', ...
                'Color',shapeColors(k,:), ...
                'LineWidth',shapeWidth, ...
                'DisplayName',sprintf('Step %d',shapeSteps(k)-1));
        end

        axis equal;
        grid off;
        box on;

        xlabel('');
        ylabel('');

        if axisMode == 1

            % Same axis limits for both models
            xlim(commonXLim);
            ylim(commonYLim);

            % No axis numbers
            xticks([]);
            yticks([]);

            title([modelName ' - Common Axes'], ...
                'FontWeight','normal');

        else

            % Individually fitted limits, with numbers
            axis tight;

            title([modelName ' - Auto Axes'], ...
                'FontWeight','normal');
        end

        legend('Location','bestoutside');
    end
end

%% ================================================================
% 4.3 - DEFORMED SHAPES: LINEAR AND NONLINEAR
% ================================================================

load('FEA_Results.mat','results');

X = results.mesh.X;
conn = results.mesh.connectivity;

nNodes = size(X,1);
nStates = size(results.nonlinear.qHistory,2);

% Three load states
shapeSteps = [2, round((nStates+1)/2), nStates];

% Deformation scale
shapeScale = 1;

% Colors
originalColor = [0.55 0.55 0.55];

linearColors = [ ...
    1.00 0.65 0.65;
    0.85 0.30 0.30;
    0.55 0.00 0.00];

nonlinearColors = [ ...
    0.65 0.80 1.00;
    0.25 0.50 0.85;
    0.00 0.15 0.55];

originalWidth = 2.5;
shapeWidth = 2;

% Element boundary ordering
if size(conn,2) == 4
    edgeOrder = [1 2 3 4 1];
elseif size(conn,2) == 8
    edgeOrder = [1 5 2 6 3 7 4 8 1];
else
    error('Only Quad4 and Quad8 supported.');
end

%% Determine common axis limits

allX = X(:,1);
allY = X(:,2);

for model = 1:2

    if model == 1
        qHistory = results.linear.qHistory;
    else
        qHistory = results.nonlinear.qHistory;
    end

    for k = 1:numel(shapeSteps)

        U = reshape(qHistory(:,shapeSteps(k)),2,nNodes).';
        Xdef = X + shapeScale*U;

        allX = [allX; Xdef(:,1)];
        allY = [allY; Xdef(:,2)];
    end
end

xPadding = 0.05*range(allX);
yPadding = 0.05*range(allY);

commonXLim = [min(allX)-xPadding, max(allX)+xPadding];
commonYLim = [min(allY)-yPadding, max(allY)+yPadding];

%% Plot linear and nonlinear separately

for model = 1:2

    if model == 1
        qHistory = results.linear.qHistory;
        shapeColors = linearColors;
        modelName = 'Linear';
    else
        qHistory = results.nonlinear.qHistory;
        shapeColors = nonlinearColors;
        modelName = 'Nonlinear';
    end

    figure('Name',[modelName ' Deformed Shapes'], ...
        'Color','w', ...
        'Units','pixels', ...
        'Position',[100 100 900 900]);

    hold on;

    % Original geometry
    for e = 1:size(conn,1)

        ids = conn(e,edgeOrder);

        plot(X(ids,1),X(ids,2), ...
            'Color',originalColor, ...
            'LineWidth',originalWidth, ...
            'HandleVisibility','off');
    end

    % Three deformed states
    for k = 1:numel(shapeSteps)

        step = shapeSteps(k);

        U = reshape(qHistory(:,step),2,nNodes).';
        Xdef = X + shapeScale*U;

        for e = 1:size(conn,1)

            ids = conn(e,edgeOrder);

            plot(Xdef(ids,1),Xdef(ids,2), ...
                'Color',shapeColors(k,:), ...
                'LineWidth',shapeWidth, ...
                'HandleVisibility','off');
        end
    end

    % Legend
    plot(nan,nan,'-', ...
        'Color',originalColor, ...
        'LineWidth',originalWidth, ...
        'DisplayName','Undeformed');

    for k = 1:numel(shapeSteps)

        step = shapeSteps(k);

        Fy = sum(results.loading.loadHistory(2:2:end,step));

        plot(nan,nan,'-', ...
            'Color',shapeColors(k,:), ...
            'LineWidth',shapeWidth, ...
            'DisplayName',sprintf( ...
            'Step %d (F_y = %.4g N)',step-1,Fy));
    end

    % Same geometry limits for both models
    xlim(commonXLim);
    ylim(commonYLim);

    axis equal;
    axis off;
    grid off;

    title([modelName ' Deformed Shapes'], ...
        'FontWeight','normal');

    legend('Location','southoutside');

end