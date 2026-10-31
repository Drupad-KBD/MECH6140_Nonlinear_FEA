classdef AbaqusINP
    %ABAQUSINP Read 2D Q4 geometry and named node/element sets ONLY.
    % Supports a single meshed instance (as in Job-1.inp).
    % Does NOT import loads, boundary conditions, materials or steps.
    %
    % mesh = fea.mesh.AbaqusINP.read('Job-1.inp');
    % mesh.X             : nNodes x 2 coordinates
    % mesh.connectivity  : nElements x 4 MATLAB node-row indices
    % mesh.nodeLabels    : original Abaqus node labels
    % mesh.elementLabels : original Abaqus element labels
    % mesh.nsets         : struct array with name, scope, labels, nodes
    % mesh.elsets        : struct array with name, scope, labels, elements
    %
    % Sets in different scopes remain separate (e.g. Set-1).

    methods (Static)
        function mesh = read(filename)
            fid = fopen(filename,'r');
            assert(fid >= 0,'Cannot open input file: %s',filename);
            cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>

            nodes = zeros(0,3);
            elems = zeros(0,5);
            nsets = struct('name',{},'scope',{},'labels',{},'nodes',{});
            elsets = struct('name',{},'scope',{},'labels',{},'elements',{});
            section = '';
            scope = 'model';
            setName = '';
            setScope = '';
            generated = false;
            elementType = '';
            instanceName = '';

            while ~feof(fid)
                line = fgetl(fid);
                if ~ischar(line), break; end
                line = strtrim(line);
                if isempty(line) || startsWith(line,'**'), continue; end

                if startsWith(line,'*')
                    section = '';
                    low = lower(line);
                    if startsWith(low,'*instance,')
                        token = regexp(line,'(?i)name\s*=\s*([^,\s]+)','tokens','once');
                        assert(~isempty(token),'Instance without name.');
                        instanceName = token{1};
                        scope = ['instance:' instanceName];
                    elseif startsWith(low,'*end instance')
                        scope = 'assembly';
                    elseif startsWith(low,'*assembly,')
                        scope = 'assembly';
                    elseif startsWith(low,'*end assembly')
                        scope = 'model';
                    elseif startsWith(low,'*node')
                        section = 'node';
                    elseif startsWith(low,'*element,')
                        section = 'element';
                        token = regexp(line,'(?i)type\s*=\s*([^,\s]+)','tokens','once');
                        if ~isempty(token), elementType = upper(token{1}); end
                    elseif startsWith(low,'*nset,') || startsWith(low,'*elset,')
                        if startsWith(low,'*nset,')
                            section = 'nset'; key = 'nset';
                        else
                            section = 'elset'; key = 'elset';
                        end
                        token = regexp(line,['(?i)' key '\s*=\s*([^,\s]+)'],'tokens','once');
                        assert(~isempty(token),'Unnamed %s encountered.',key);
                        setName = token{1};
                        generated = contains(low,'generate');
                        inst = regexp(line,'(?i)instance\s*=\s*([^,\s]+)','tokens','once');
                        if isempty(inst)
                            setScope = scope;
                        else
                            setScope = ['assembly:' inst{1}];
                        end
                    end
                    continue;
                end

                parts = strtrim(strsplit(line,','));
                parts = parts(~cellfun('isempty',parts));
                v = str2double(parts);
                switch section
                    case 'node'
                        assert(numel(v)>=3 && all(isfinite(v(1:3))), ...
                            'Invalid node row: %s',line);
                        nodes(end+1,:) = v(1:3); %#ok<AGROW>
                    case 'element'
                        assert(numel(v)==5 && all(isfinite(v)), ...
                            'Only Q4 connectivity is supported: %s',line);
                        elems(end+1,:) = v; %#ok<AGROW>
                    case {'nset','elset'}
                        assert(all(isfinite(v)),'Invalid set row: %s',line);
                        if generated
                            assert(numel(v)==3 && v(3)>0 && v(2)>=v(1), ...
                                'Invalid generated set row: %s',line);
                            labels = (v(1):v(3):v(2)).';
                        else
                            labels = v(:);
                        end
                        if strcmp(section,'nset')
                            idx = find(strcmpi({nsets.name},setName) & ...
                                strcmpi({nsets.scope},setScope),1);
                            if isempty(idx)
                                idx = numel(nsets)+1;
                                nsets(idx).name = setName;
                                nsets(idx).scope = setScope;
                                nsets(idx).labels = [];
                                nsets(idx).nodes = [];
                            end
                            nsets(idx).labels = unique([nsets(idx).labels(:);labels]);
                        else
                            idx = find(strcmpi({elsets.name},setName) & ...
                                strcmpi({elsets.scope},setScope),1);
                            if isempty(idx)
                                idx = numel(elsets)+1;
                                elsets(idx).name = setName;
                                elsets(idx).scope = setScope;
                                elsets(idx).labels = [];
                                elsets(idx).elements = [];
                            end
                            elsets(idx).labels = unique([elsets(idx).labels(:);labels]);
                        end
                end
            end

            assert(~isempty(nodes) && ~isempty(elems),'Missing mesh nodes/elements.');
            assert(ismember(elementType,{'CPE4','CPE4R','CPS4','CPS4R'}), ...
                'Unsupported element type: %s',elementType);
            nodeLabels = nodes(:,1);
            elementLabels = elems(:,1);
            assert(numel(unique(nodeLabels))==numel(nodeLabels),'Duplicate node labels.');
            assert(numel(unique(elementLabels))==numel(elementLabels),'Duplicate element labels.');
            [ok,conn] = ismember(elems(:,2:5),nodeLabels);
            assert(all(ok(:)),'Element references unknown node labels.');
            X = nodes(:,2:3);

            % Ensure positive signed area / counterclockwise local Q4 order.
            for e=1:size(conn,1)
                xy = X(conn(e,:),:);
                twiceArea = sum(xy(:,1).*xy([2 3 4 1],2) - ...
                                xy([2 3 4 1],1).*xy(:,2));
                assert(abs(twiceArea)>1e-13,'Zero-area element %d.',elementLabels(e));
                if twiceArea<0
                    conn(e,:) = conn(e,[1 4 3 2]);
                end
            end

            % Convert original Abaqus labels in sets to MATLAB row indices.
            for k=1:numel(nsets)
                [found,rows] = ismember(nsets(k).labels,nodeLabels);
                assert(all(found),'Undefined node in set %s (%s).', ...
                    nsets(k).name,nsets(k).scope);
                nsets(k).nodes = rows(:);
            end
            for k=1:numel(elsets)
                [found,rows] = ismember(elsets(k).labels,elementLabels);
                assert(all(found),'Undefined element in set %s (%s).', ...
                    elsets(k).name,elsets(k).scope);
                elsets(k).elements = rows(:);
            end

            mesh.X = X;
            mesh.connectivity = conn;
            mesh.nodeLabels = nodeLabels;
            mesh.elementLabels = elementLabels;
            mesh.nsets = nsets;
            mesh.elsets = elsets;
            mesh.elementType = elementType;
            mesh.nNodes = size(X,1);
            mesh.nElements = size(conn,1);
            fprintf('Mesh OK: %d nodes, %d Q4 elements (%s), %d node sets, %d element sets.\n', ...
                mesh.nNodes,mesh.nElements,elementType,numel(nsets),numel(elsets));
        end
    end
end
