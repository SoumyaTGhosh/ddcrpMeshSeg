function [error,tri_err] = measure_nonrigidity(test_set,labels)
    test_pose = test_set.Y;
    reference_pose = test_set.X;
    error = Inf*ones(1,length(test_pose));
    per_trierr = zeros(length(labels),length(test_pose));
    for i = 1:length(test_pose)
        [error(i),per_trierr(:,i)] = evaluate_segmentation_perpose(test_pose{i}(:,:,1),reference_pose{i},labels);
    end
%     for i = 1:size(test_pose,3)
%         [error(i),per_trierr(:,i)] = evaluate_segmentation_perpose(test_pose{1}(:,:,i),reference_pose{1},labels);
%     end 
    tri_err = sum(per_trierr,2); 
end

function [score,pertri_err] = evaluate_segmentation_perpose(test_pose,reference_pose,labels)
% test_pose, reference_pose are N*3 matrices, where N is the number of
% triangles constituiting the mesh. labels is N*1.
  
    num_parts = length(unique(labels));
    part_score = zeros(1,num_parts);
    pertri_err  = NaN(1,length(labels));
    design = [1 1 1 1 0 0 0 0 0 0 0 0;...
              0 0 0 0 1 1 1 1 0 0 0 0;...
              0 0 0 0 0 0 0 0 1 1 1 1];
    ctr = 1;    
    for k = unique(labels)
        % select the triangles belonging to part k
        ind = labels==k;
        zeropadx = zeros(3*length(find(ind>0)),12);  
        y = test_pose(ind,:);
        x = reference_pose(ind,:);
        % stack y row wise
        y = y';
        y = y(:);
        % construct zero padded x matrix
        for i = 1: size(x,1)
            x_row = [x(i,:) 1]; 
            zeropadx(3*(i-1)+1:3*(i-1)+3,:) = repmat(x_row,[3,3]).*design;            
        end
        % least squares minimization
        a = zeropadx\y;
        % compute residuals
        residual = zeropadx*a-y;
        part_score(ctr) = norm(residual,2);
        % compute per face residuals 
        tri_ctr = 1;
        for tri = find(ind)
            pertri_err(tri) = norm(residual((tri_ctr-1)*3+(1:3)),2);
            tri_ctr = tri_ctr +1;
        end
        
        ctr = ctr + 1;
    end
    % sum over all parts.
    score = sum(part_score);
end