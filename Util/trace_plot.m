function [ln_map_e,ll,lp]=trace_plot(iter,mniw,data,sampler_state,D)

ln_map_e_store = zeros(1,sampler_state.T(iter));
p = zeros(sampler_state.num_data,1);
ctr = 1;
for t=unique(sampler_state.t(iter,:))
    tab_members = sampler_state.t(iter,:)==t;
    ln_map_e_store(ctr) = compute_table_lik(mniw,data,tab_members);
    ctr = ctr + 1;
end


nc = sum(D,2);
for i = 1:sampler_state.num_data
    p(i) = D(i,sampler_state.c(iter,i));
end
prior = log(p./nc);
% for i = 1:length(param.c)
%     c = param.c(i);
%     nc = nnz(D(i,:));
%     prior(i) = log(D(i,c)/nc);
% end
lp = sum(prior);

ll = sum(ln_map_e_store);

ln_map_e = ll + lp ;%- 1.5*sampler_state.num_data*log(pi);

end