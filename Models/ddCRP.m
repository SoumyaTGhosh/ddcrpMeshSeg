classdef ddCRP
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        %self link probability
        alpha
        %pairwise distance matrix
        d2_mat
    end
    
    methods
        function obj = ddCRP(alpha,d2_mat)
            obj.alpha = alpha;
            obj.d2_mat = d2_mat;
        end
    end
    
end

