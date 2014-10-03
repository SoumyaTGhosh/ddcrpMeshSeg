classdef SamplerState
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        c
        t
        T
        num_data
        
    end
    
    methods
        function obj = SamplerState(data,num_iter)
            % number of mesh faces (data instances)
            obj.num_data=size(data.X{1},1);
            obj.c = NaN(num_iter,obj.num_data);
            obj.t = NaN(num_iter,obj.num_data);
            obj.T = NaN(1,num_iter);
      
        end
        
    end
    methods
        function obj = initialize(obj)
            % initialize each meshface as its own part.
            obj.c(1,:)=1:obj.num_data;
            obj.t(1,:)=1:obj.num_data;
            % number of active tables
            obj.T=obj.num_data;
        end
    end
    
end

