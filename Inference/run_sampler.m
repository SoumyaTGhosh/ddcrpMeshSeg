function sampler_state=run_sampler(prior_param,mniw,data,sampler_state,bookkeeper,num_iter)


savedir = bookkeeper.savedir;



for iter=2:num_iter
    
    sampler_state.c(iter,:)=sampler_state.c(iter-1,:);
    sampler_state.t(iter,:)=sampler_state.t(iter-1,:);
    sampler_state.T(iter)=sampler_state.T(iter-1);
    
    tic
    for i_s=1:sampler_state.num_data
        [sampler_state,bookkeeper] = sample_links(iter,i_s,sampler_state,bookkeeper,...
            prior_param,mniw,data);
    end
    t =toc;
    
    fprintf('. ');
    %fprintf('Iteration - %d took %0.03f seconds to complete\n',iter,t);
    if(mod(iter,20)==0)
        save('-v7.3',fullfile(savedir,sprintf('Intermediate_MCMC_sample%d',iter)),'sampler_state');
        fprintf('Iteration - %d took %0.03f seconds to complete\n',iter,t);
    end
    
    
end


