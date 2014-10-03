function [table_lik] = compute_table_lik(mniw,data,table_members)
%   computes table marginal likelihood.

    table_lik_bodies = zeros(1,length(data.num_bodies));
    for i = 1:data.num_bodies % each body refers to a different person
        % different body shape and hence different reference body.
%         X = data.X{i}(table_members,:)';
%         Y = data.Y{i}(:,table_members,:);

        table_lik_bodies(i) = mniw.computeMarginalLik(data.Y{i}(:,table_members,:),data.X{i}(table_members,:)');
%         [table_lik_bodies(i),~] = computeMultiPoseMNIWlik(Y,X,param);


        %part likelihood is computed by summing over all bodies.
        table_lik = sum(table_lik_bodies);

    end
end
