function demo_segmentation
% Demo for mesh segmentation
% Reads from disk (./Data) meshes preprocessed into point clouds.
% Writes to disk (./Results) the resulting segmentation.

%reset random number seed
rng shuffle;

%%%%%%%%% SET PARAMETERS %%%%%%%%%%%%%%%%%%%%%%
% 1) number of sampler iterations. Higher the better.
% Recommended setting > 150. Reasonable segmentations can often
% be achieved within the first 30 iterations. While these premature samples
% do not correspond to samples from the posterior distribution, they might 
% be sufficient to produce reasonable results. The segmentation quality
% gradually increases with increasing number of iterations.
NUMITER = 35; 
% 2) Select dataset to analyze
% We include two sets of meshes from the TOSCA dataset (Centaur and Horse).
% Each dataset consists of a single synthetic mesh in several different poses.
% We will segment 
DATASET = 'Centaur'; % could also select 'Horse'
% 3) Noise parameter, this controls the amount of deviation from the
% predicted affine transformation.
% Noise --- HIGHER values will result in LARGER segments. 
NOISESCALE = 9e-04;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% Setup paths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pathtolightspeed = '#### TODO: POINT TO LIGHTSPEED TOOLBOX #####';
Pathtomatlab_bgl = '#### TODO: POINT TO MATLAB_BGL TOOLBOX #####';
addpath(genpath('.'));
addpath(genpath(Pathtolightspeed));
addpath(genpath(Pathtomatlab_bgl));

savedir = './Results'; 
if(~exist(savedir,'dir'));mkdir(savedir);end

datapath = './Data';

%%%%%%%% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load Mesh face neighbors, symmetric mapping along the midsagittal plane , 
%and mesh faces as point clouds.  
[A,left_tris,pcf] = load_data(datapath,DATASET);

A = A(left_tris,left_tris);
% create a ddCRP object.
ddcrp = ddCRP(double(full(A(1,1))),double(A));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for i = 1:length(pcf.X)
    % reference coords.
    X{i} = [pcf.X{i} ones(size(pcf.X{i},1),1)]; 
end
Y = pcf.Y; % these are the coordinates of triangles in different body poses.

% visualize the reference mesh.
figure;plot3(pcf.X{1}(:,1),pcf.X{1}(:,2),pcf.X{1}(:,3),'r*');title('Reference mesh visualized as a point cloud');


view(94,10);axis equal;
pause(5.0);
close all


[ data ] = package_data(X,Y,1);

% load precomputed ratio of gammas to speed up likelihood computation.
load(fullfile(datapath,'Ratio_of_Gammas.mat'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%Set hyper-parameters governing the affine likelihoods
%(See http://cs.brown.edu/~sghosh/papers/GhoshSudderthLoperBlack12NIPS.pdf for details)%%%%%%%% 

% Expected variance
S0 = NOISESCALE*eye(3);
n0 = 5;

% mean transformation
M = [1 0 0 0;0 1 0 0;0 0 1 0];

%PRECISON matrix. Higher = smaller var.
K = NOISESCALE*[1 0 0 0; 0 1 0 0;0 0 1 0;0 0 0 0.1]; 

% Matrix Normal Inverse Wishart likelihood object
mniw = MNIWlik(NOISESCALE,S0,n0,M,K,Ratio_of_Gammas.rofg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%% Initialize and Run Sampler%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Gibbs Sampler will be run for %d iterations.\n',NUMITER);

% Create sampler state object
sampler_state = SamplerState(data,NUMITER);
% Initialize Sampler
sampler_state=sampler_state.initialize();
% Create bookkeeping object
bookkeeper = BookKeeper(sampler_state,savedir);

% Run Sampler
converged_sampler_state=run_sampler(ddcrp,mniw,data,sampler_state...
                                    ,bookkeeper,NUMITER);
% save results
save('-v7.3',fullfile(savedir,'MCMC_Chain.mat'),'converged_sampler_state');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%% Diagnostics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% joint lik trace
ll_c = NaN(1,NUMITER);
ll = NaN(1,NUMITER);
prior = NaN(1,NUMITER);
for i = 1:NUMITER
    [ll_c(i) ll(i) prior(i)] = trace_plot(i,mniw,data,converged_sampler_state,ddcrp.d2_mat);
end;
figure(3);plot(ll_c,'r*-');title('joint log-lik')
[~,MAP_SAMPLE] = max(ll_c);

% visualize result 
visualize_result(2,converged_sampler_state,MAP_SAMPLE,left_tris,DATASET);
save(fullfile(savedir,'State.mat'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

