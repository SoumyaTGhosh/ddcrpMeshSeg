function [li,bookkeeper]=compute_likelihood_ratio(mniw,data,bookkeeper,t_prop,t_curr,t_prop_members,t_curr_members)
%t_curr is the table that super-pixel i currently belongs to, t_prop is the
%table being proposed for i.

if(t_prop==t_curr)
    %the new link doesn't merge tables
    %li (log-likelihood) will be zero. and we will sample only on the basis of the prior
    li = 0;
else
    % if we don't already have table likelihoods compute and return them along with their counts.
    if(~bookkeeper.table_lik(bookkeeper.valid_clusters(t_prop)))
        bookkeeper.table_lik(bookkeeper.valid_clusters(t_prop)) = compute_table_lik(mniw,data,t_prop_members);
        
    end
    if(~bookkeeper.table_lik(bookkeeper.valid_clusters(t_curr)))
        bookkeeper.table_lik(bookkeeper.valid_clusters(t_curr)) = compute_table_lik(mniw,data,t_curr_members);
    end
    min_t = min(bookkeeper.valid_clusters(t_prop),bookkeeper.valid_clusters(t_curr));
    max_t = max(bookkeeper.valid_clusters(t_prop),bookkeeper.valid_clusters(t_curr));
    if(bookkeeper.pairwise_table_lik(min_t,max_t))
        merged_lik = bookkeeper.pairwise_table_lik(min_t,max_t);
    else
        %merged_lik not on record. recompute.
        %note that for merged_lik computation the table_counts have already been computed.

        merged_lik_bodies = zeros(1,data.num_bodies);
        for i = 1:data.num_bodies
            %compute likelihoods across different body shapes, each
            %having its own reference pose.
      
            merged_lik_bodies(i) = mniw.computeMarginalLik([data.Y{i}(:,t_prop_members,:) data.Y{i}(:,t_curr_members,:)],...
                [data.X{i}(t_prop_members,:)',data.X{i}(t_curr_members,:)']);
        end
        %compute part likelihood by summing over all bodies.
        merged_lik = sum(merged_lik_bodies);
        bookkeeper.pairwise_table_lik(min_t,max_t) = merged_lik;
    end
    li= merged_lik...
        -bookkeeper.table_lik(bookkeeper.valid_clusters(t_prop)) -bookkeeper.table_lik(bookkeeper.valid_clusters(t_curr)); % p(t_old U t_new)/p(t_old)*p(t_new);
    
end


end
