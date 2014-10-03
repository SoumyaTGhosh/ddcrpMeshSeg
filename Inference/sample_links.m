function [sampler_state,bookkeeper] = sample_links(iter,i_s,sampler_state,bookkeeper,prior_param,mniw,data)
%%%  Sample p(c_i|c_ni,X,Y)

runCC = check_shallow_cycles(sampler_state.c(iter,:),i_s);

temp_c=sampler_state.c(iter,:);
%reset the current link.
temp_c(i_s) = i_s;
sampler_state.c(iter,i_s) = i_s;

if(runCC)

    % should only have to compute c_to_t for the current table.
    table = sampler_state.t(iter,i_s);
    table_members = find(sampler_state.t(iter,:)==table);
    inverted_index = sparse(1,table_members,1:length(table_members));
    
    %temp_c(table_members) are the customers tab_members are pointing to.
    %remapped_table_members contains these customers mapped to 1:table_members
    remapped_table_members = inverted_index(temp_c(table_members));

    %Perform a connected components operation and determine if restting the
    %current link splits the corresponding table.
     [tables,split] = ...
        fast_cc(length(remapped_table_members),remapped_table_members);
else
    %shallow cycle found
    split = 1;
end

if(split==2)
    split_customers = table_members(tables==2);
    sampler_state.t(iter,split_customers) = max(sampler_state.t(iter,:))+1;
    num_tables = sampler_state.T(iter)+1;
    sampler_state.T(iter) = num_tables;
    bookkeeper.valid_clusters(num_tables) = num_tables;
    %reset the likelihood of the split tables
    if (bookkeeper.table_lik(bookkeeper.valid_clusters(sampler_state.T(iter))))       
       bookkeeper.table_lik(bookkeeper.valid_clusters(sampler_state.T(iter))) = 0; 
    end
    bookkeeper.table_lik(table) = 0;
    
    %reset pairwise table_lik    
    bookkeeper.pairwise_table_lik(bookkeeper.valid_clusters(sampler_state.T(iter)),:) = bookkeeper.reset_vec;
    bookkeeper.pairwise_table_lik(bookkeeper.valid_clusters(table),:) = bookkeeper.reset_vec;
    bookkeeper.pairwise_table_lik(:,bookkeeper.valid_clusters(sampler_state.T(iter))) = bookkeeper.reset_vec';
    bookkeeper.pairwise_table_lik(:,bookkeeper.valid_clusters(table)) = bookkeeper.reset_vec';
    
else
%     param.t = param.t;
%     %table_lik remains unchanged as well
%     num_tables = param.T;
end
if(split>2)
    disp('Err - Cant split table more than two ways');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sample new link for the current mesh face 
[nbor,~,prior] = find(prior_param.d2_mat(:,i_s));
logposterior = log(prior);

t_curr = sampler_state.t(iter,i_s);
t_curr_members = sampler_state.t(iter,:)==t_curr;
for ii = 1:length(nbor)
    t_prop=sampler_state.t(iter,nbor(ii));
    t_prop_members = sampler_state.t(iter,:)==t_prop;
    [loglik_ratio,bookkeeper]=compute_likelihood_ratio(mniw,data,bookkeeper,...
        t_prop,t_curr,t_prop_members,t_curr_members);
    logposterior(ii) = loglik_ratio+logposterior(ii);
end

logposterior=logposterior-max(logposterior);
logposterior=exp(logposterior);
logposterior=logposterior./sum(logposterior);


c_ji_new = nbor(sample(logposterior,1));
%update customer links
sampler_state.c(iter,i_s) = c_ji_new;

%update tables
link_table = sampler_state.t(iter,c_ji_new);
curr_tab = sampler_state.t(iter,i_s); %note this won't be the same as tab, if a split occurs.
if(link_table==curr_tab)
    % do nothing, linking to a customer in the current table

else
    large_idx = max(link_table,curr_tab);
    small_idx = min(link_table,curr_tab);
    ids = (sampler_state.t(iter,:)==large_idx);
    sampler_state.t(iter,ids) = small_idx;
   
    ids = sampler_state.t(iter,:)>large_idx;
    sampler_state.t(iter,ids) = sampler_state.t(iter,ids) - 1;
    bookkeeper.valid_clusters(bookkeeper.valid_clusters>large_idx) = bookkeeper.valid_clusters(bookkeeper.valid_clusters>large_idx) - 1;
%   Modify table_lik to reflect the merge.
    % reset the likelihood of the merged table.
    loc = bookkeeper.valid_clusters(small_idx);
    bookkeeper.table_lik(loc) = 0;

    bookkeeper.valid_clusters(large_idx:sampler_state.T(iter)-1) = bookkeeper.valid_clusters(large_idx+1:sampler_state.T(iter)); 
    bookkeeper.pairwise_table_lik(:,loc) = bookkeeper.reset_vec';
    bookkeeper.pairwise_table_lik(loc,:) = bookkeeper.reset_vec;

    sampler_state.T(iter) = sampler_state.T(iter)-1;
       
end


if (sampler_state.T(iter)==0)
    disp('something is wrong');
end




end

function runCC = check_shallow_cycles(curr_c,i_s)
% check for shallow cycles. Presence of a cycle ensures that removing the
% current link won't cause a split. Saves us the expensive connected
% components operation.
%       
runCC = 1;
if(curr_c(i_s)==i_s)
    % if the existing link is a self-link, removing it can't split a
    % component.
    runCC = 0;
else
    %%% check shallow cycles
    % A->B; B->A, removing either link can't split the component.
    if(curr_c(curr_c(i_s))==i_s)
        runCC = 0;
        
        % A->B->C; C->A, removing any of these links can't split the component.
    elseif(curr_c(curr_c(curr_c(i_s)))==i_s)
        runCC = 0;
        
        % A->B->C->D; D->A, removing any of these links can't split the component.
    elseif(curr_c(curr_c(curr_c(curr_c(i_s))))==i_s)
        runCC = 0;
        
        % A->B->C->D->E; E->A, removing any of these links can't split the component.
    elseif(curr_c(curr_c(curr_c(curr_c(curr_c(i_s)))))==i_s)
        runCC = 0;
        
        % 6 cycle, removing any of these links can't split the component.
    elseif(curr_c(curr_c(curr_c(curr_c(curr_c(curr_c(i_s))))))==i_s)
        runCC = 0;
        
        % 7 cycle
    elseif(curr_c(curr_c(curr_c(curr_c(curr_c(curr_c(curr_c(i_s)))))))==i_s)
        runCC = 0;
        
        % 8 cycle
    elseif(curr_c(curr_c(curr_c(curr_c(curr_c(curr_c(curr_c(curr_c(i_s))))))))==i_s)
        runCC = 0;
        
        % 9 cycle
    elseif(curr_c(curr_c(curr_c(curr_c(curr_c(curr_c(curr_c(curr_c(curr_c(i_s)))))))))==i_s)
        runCC = 0;
        
        % 10 cycle
    elseif(curr_c(curr_c(curr_c(curr_c(curr_c(curr_c(curr_c(curr_c(curr_c(i_s)))))))))==i_s)
        runCC = 0;
    end
    
    
end
end