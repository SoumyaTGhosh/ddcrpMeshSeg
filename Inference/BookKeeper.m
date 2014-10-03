classdef BookKeeper
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Access = public)
        valid_clusters
        table_lik
        pairwise_table_lik
        reset_vec
        savedir
    end
    
    methods
        function obj = BookKeeper(sampler_state,savedir)
            obj.valid_clusters = 1:sampler_state.num_data;
            obj.table_lik = sparse(1,sampler_state.num_data);
            obj.pairwise_table_lik = sparse(sampler_state.num_data,sampler_state.num_data);
            obj.reset_vec = sparse(1,size(obj.pairwise_table_lik,1));
            
            %
            obj.savedir = savedir;
        end
    end
    
end

