classdef MNIWlik
    % Matrix Normal Inverse Wishart Likelihood class
    %   Detailed explanation goes here
    
    properties
        %%---Inverse Wishart parameters---%%
        %Scale
        scale;
        % Expected variance
        S0;
        % Dof
        n0 = 5;
        
        %%---Matrix Normal paramters---%%
        % Mean transformation
        M;
        
        %PRECISON matrix. Higher = smaller var.
        K; 
        
        %%---Precomputations useful for Marginal likelihood computation---%%
        %Precomputed Ratio of Gammas (\Gamma_3((N_k+n0)/2)/\Gamma_3(n0/2))
        rofg;
        %
        I = eye(4);
        MK;
        MKMt;
        sumKS0;
        
        
    end
    
    methods
        function mniw = MNIWlik(scale,S0,n0,M,K,rofg)
            
            mniw.S0 = S0;
            mniw.n0 = n0;
            mniw.M = M;
            mniw.K = K;
            mniw.rofg = rofg;
            mniw.scale = scale;
            % precompute values needed for Marginal likelihood computations
            L1 = chol(mniw.S0,'lower');
            logdet_S0 = 2*sum(log(diag(L1)));
            L1 = chol(mniw.K,'lower');
            logdet_K = 2*sum(log(diag(L1)));
            mniw.MKMt = M*mniw.K*M';
            mniw.MK = M*mniw.K;
            mniw.sumKS0 = (3/2)*logdet_K+ (n0/2)*logdet_S0;
            
        end
    end
    methods 
        function llik = computeMarginalLik(mniw,Y_p,X)
            % I/P - hyper-params and data
            % O/P - log marginal likelihood
            
            
            %X \in R^4*N, Y_p \in R^3*N*num_poses
           
            num_poses = size(Y_p,3);
           
            N = size(Y_p,2);
            
            
            % Sxx = X*X'+ I;
            Sxx = X*X'+ mniw.K;
            L = chol(Sxx,'lower');
            invL = L\mniw.I;
            T = invL*X;
            b = invL*mniw.MK';
            B = mniw.MKMt - b'*b;
            
            
            %compute pose independent terms of mlik;
            
            C1 = mniw.sumKS0 ...
                -(3/2)* 2*sum(log(diag(L)));
            %ratio %of 3 dimensional gamma functions
            C2 = mniw.rofg(N);
            C = C1+C2;
            %compute pose dependent terms of the mlik.
            poselik = zeros(1,num_poses);
         
            for p = 1:num_poses
                t = Y_p(:,:,p)*T';
                btcrossterm = t*b;
                Syx = Y_p(:,:,p)*Y_p(:,:,p)' + B -t*t' -btcrossterm' - btcrossterm; % This is what needs to be computed for each pose Syx = Y*Y' + MKMt -t*t' -b'*b - (t*b)'-(t*b);
%                 Syx = Y_p(:,:,p)*Y_p(:,:,p)' + B -t*t'; % This is what needs to be computed for each pose Syx = Y*Y' + MKMt -t*t' -b'*b;
                
                L1 = chol(mniw.S0+Syx,'lower');
                logdet_S0Syx = 2*sum(log(diag(L1)));
                
                poselik(p) = C-((N+mniw.n0)/2)*logdet_S0Syx;
            end
   
            llik = sum(poselik,2); 
            llik = llik - 1.5*N*num_poses*log(pi);
        end
    end
    
end

