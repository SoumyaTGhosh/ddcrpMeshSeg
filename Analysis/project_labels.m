function labels = project_labels(mniw,param,data,meshpath)
    % Using left_right body symmetry project the half body segmentation 
    % results to the full body.
    
    % load symmetric mapping
    load(fullfile(meshpath,'tri2opposite.mat'));
    
    
    left_bodies = 1:2:2*data.num_bodies;
    right_bodies = 2:2:2*data.num_bodies;
    data_left.X = data.X(left_bodies);
    data_left.Y = data.Y(left_bodies);
    data_right.X = data.X(right_bodies);
    data_right.Y = data.Y(right_bodies);
    data_left.num_bodies = length(left_bodies);
    data_right.num_bodies = length(right_bodies);
    
    % concatenate into a whole body dataset
    clear data;
    for i = 1:data_left.num_bodies
        data.X{i} = [data_left.X{i};data_right.X{i}];
        data.Y{i} = [data_left.Y{i} data_right.Y{i}];
    end
    data.num_bodies = data_left.num_bodies;
    
    for i = 1:param.T
        table_members = param.t==i;
        %independent likelihoods
        ilik(i) = compute_table_lik(mniw,data_left,table_members) + compute_table_lik(mniw,data_right,table_members);
      %  plot3(data_left.X{1}(find(table_members),1),data_left.X{1}(find(table_members),3),data_left.X{1}(find(table_members),2),'r*');
        %joint likelihoods
        % project labels onto left half
        whole_table_members = (zeros(1,size(data.X{1},1)));
        whole_table_members(find(table_members)) = 1;
        %whole_table_members(tri2opposite.tri2opposite(find(table_members))) = 1;
        whole_table_members(find(table_members)+param.num_data) = 1;
        %labels(whole_table_members) = i;
       % figure;plot3(data.X{1}(find(whole_table_members),1),data.X{1}(find(whole_table_members),3),data.X{1}(find(whole_table_members),2),'r*');
        jlik(i) = compute_table_lik(mniw,data,logical(whole_table_members));

    end
    parts_to_split = find(ilik>jlik);
    %%% randpermute the order of labels to make colormap more meaningful.
    sa_m = 1;
    randt = zeros(1,length(param.t(sa_m,:)));
    mc1 = randperm(length(unique(param.t(sa_m,:))));
    for ii = 1:length(unique(param.t(sa_m,:)))
        randt(param.t(sa_m,:)==ii) = mc1(ii);
    end;
    
    right_labels = randt;
    left_labels = randt;
    for p = parts_to_split
        %right_labels(param.t>p) = right_labels(param.t>p) + 2;
        left_labels(param.t==p) = left_labels(param.t==p) + 2;
    end
    
    left_tris = tri2opposite.left_tris;
    labels = zeros(1,2*length(left_tris));
    labels(left_tris) = left_labels;
    labels(tri2opposite.tri2opposite(left_tris)) = right_labels;
    
    fnames =  dir(fullfile(meshpath,'meshes','*.ply'));
    
    [vrts,faces] = readply(fullfile(meshpath,'meshes',fnames(1).name));
    vrts = vrts';
    
    
    figure(2);t =trisurf(faces,vrts(:,1),vrts(:,3),vrts(:,2),labels,'EdgeColor','none');%title(sprintf('%d',sa_m));
    set(t,'FaceLighting','phong',...
        'AmbientStrength',0.5);
    
    light('Position',[5 0 5],'Style','infinite', 'Color', [.25 .25 .25]);
    light('Position',[0 5 5],'Style','infinite', 'Color', [.25 .25 .25]);
    light('Position',[-5 -5 5],'Style','infinite', 'Color', [.25 .25 .25]);
    
    axis equal;
    
    axis off;
    set(gcf, 'Color', [1 1 1]);
    view(-130,25);
   
    title('Segmentation after post processing');

