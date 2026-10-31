
classdef FieldPlotter
    %FIELDPLOTTER General Q4 contour visualization.
    %
    % Plots one solution at a time.
    % Supports displacement, strain, stress, and other scalar fields.
    %
    % INPUTS:
    % X      : nNodes x 2 reference coordinates
    % conn   : nElements x 4 Q4 connectivity
    % fields : nValues x nFields
    % opts   : plotting options
    %
    % opts.location = 'nodal' or 'element'
    %
    % Nodal fields are interpolated.
    % Element fields are piecewise constant.

    methods (Static)

        function fig = plot(X, conn, fields, opts)

            if nargin < 4 || isempty(opts)
                opts = struct();
            end

            %% DEFAULT SETTINGS

            if ~isfield(opts,'location')
                opts.location = 'nodal';
            end

            if ~isfield(opts,'colormap')
                opts.colormap = 'jet';
            end

            if ~isfield(opts,'scale')
                opts.scale = 1;
            end

            if ~isfield(opts,'showNodes')
                opts.showNodes = false;
            end

            if ~isfield(opts,'figureName')
                opts.figureName = 'FEA Field Contours';
            end
            if ~isfield(opts,'showAxisLabels')
                opts.showAxisLabels = true;
            end

            %% VALIDATE INPUTS

            nNodes = size(X,1);
            nElements = size(conn,1);
            nFields = size(fields,2);

            assert(size(X,2)==2, ...
                'Coordinates must have two columns.');

            assert(size(conn,2)==4, ...
                'Q4 connectivity must have four columns.');

            assert(all(conn(:)>=1 & conn(:)<=nNodes & ...
                conn(:)==fix(conn(:))), ...
                'Invalid connectivity.');

            assert(all(isfinite(X(:))) && ...
                all(isfinite(fields(:))), ...
                'Coordinates and field values must be finite.');

            assert(nFields>=1, ...
                'At least one scalar field is required.');

            location = lower(char(opts.location));

            switch location

                case 'nodal'
                    assert(size(fields,1)==nNodes, ...
                        'Nodal fields require one row per node.');

                case 'element'
                    assert(size(fields,1)==nElements, ...
                        'Element fields require one row per element.');

                otherwise
                    error('Location must be nodal or element.');
            end

            %% FIELD NAMES

            if ~isfield(opts,'fieldNames')
                opts.fieldNames = arrayfun( ...
                    @(k)sprintf('Field %d',k), ...
                    1:nFields, ...
                    'UniformOutput',false);
            end

            assert(numel(opts.fieldNames)==nFields, ...
                'Incorrect number of field names.');

            %% DISPLACEMENT FOR PLOTTING GEOMETRY

            if ~isfield(opts,'q')
                opts.q = zeros(2*nNodes,1);
            end

            assert(numel(opts.q)==2*nNodes, ...
                'Incorrect displacement vector length.');

            u = reshape(opts.q(:),2,nNodes).';

            Xplot = X + opts.scale*u;

            %% OPTIONAL MANUAL COLOR LIMITS

            useManualLimits = ...
                isfield(opts,'colorLimits') && ...
                ~isempty(opts.colorLimits);

            if useManualLimits
                assert(isequal(size(opts.colorLimits),[nFields 2]), ...
                    'colorLimits must be nFields x 2.');
            end

            %% CREATE FIGURE

            fig = figure( ...
                'Name',opts.figureName, ...
                'Color','w');

            layout = tiledlayout(fig,1,nFields, ...
                'TileSpacing','compact', ...
                'Padding','compact');

            %% PLOT EACH FIELD

            for k = 1:nFields

                ax = nexttile(layout,k);
                values = fields(:,k);

                %% Plot nodal or element field

                if strcmp(location,'nodal')

                    patch(ax, ...
                        'Faces',conn, ...
                        'Vertices',Xplot, ...
                        'FaceVertexCData',values, ...
                        'FaceColor','interp', ...
                        'EdgeColor',[0.15 0.15 0.15], ...
                        'LineWidth',1.1);

                else

                    hold(ax,'on');

                    for e = 1:nElements

                        patch(ax, ...
                            'Faces',1:4, ...
                            'Vertices',Xplot(conn(e,:),:), ...
                            'FaceVertexCData', ...
                                repmat(values(e),4,1), ...
                            'FaceColor','flat', ...
                            'EdgeColor',[0.15 0.15 0.15], ...
                            'LineWidth',1.1);
                    end
                end

                %% COLORMAP AND COLOR LIMITS

                colormap(ax,opts.colormap);

                if useManualLimits
                    limits = opts.colorLimits(k,:);
                else
                    limits = [min(values),max(values)];
                end

                % Avoid identical color limits for constant fields
                if limits(2)<=limits(1)

                    padding = max(1e-12, ...
                        0.01*max(1,abs(limits(1))));

                    limits = limits(1)+[-padding padding];
                end

                clim(ax,limits);

                %% COLORBAR WITH EXACT MINIMUM AND MAXIMUM LABELS

                cb = colorbar(ax);

                ticks = linspace(limits(1),limits(2),6);

                cb.Ticks = ticks;
                cb.TickLabels = compose('%.3g',ticks);

                %% OPTIONAL NODE NUMBERS

                if opts.showNodes

                    hold(ax,'on');

                    for node = 1:nNodes

                        text(ax, ...
                            Xplot(node,1), ...
                            Xplot(node,2), ...
                            sprintf(' %d',node), ...
                            'FontSize',9, ...
                            'Color','k', ...
                            'Clipping','on');
                    end
                end

                %% AXES AND LABELS

                axis(ax,'equal');
                axis(ax,'tight');

                if opts.showAxisLabels
                    xlabel(ax,'X');
                    ylabel(ax,'Y');
                    ax.XTickMode = 'auto';
                    ax.YTickMode = 'auto';
                else
                    xlabel(ax,'');
                    ylabel(ax,'');
                    ax.XTick = [];
                    ax.YTick = [];
                end
                title(ax,opts.fieldNames{k}, ...
                    'FontWeight','normal');

                box(ax,'on');
            end
        end
    end
end
