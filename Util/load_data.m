function [A,left_tris,pcf] = load_data(datapath,dataset_id)

%load pairwise distance matrix%%%%%%%%%%%%%%
load(fullfile(datapath,sprintf('adjacency_matrix_Tosca%s.mat',dataset_id)),'A');

%load symmetry mapping%%%%%%%%%%%%%
load(fullfile(datapath,sprintf('tri2opposite_Tosca%s.mat',dataset_id)));
left_tris = tri2opposite.left_tris;

%load mesh-face point clouds. Each mesh face is represented as a single 3D
%point -- the mean of its vertices.
load(fullfile(datapath,sprintf('PointCloudFeatures_Tosca%s.mat',dataset_id)));