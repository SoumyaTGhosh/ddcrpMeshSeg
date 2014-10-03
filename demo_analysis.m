function demo_analysis(DATASET)
    % Measures and visualizes non rigid regions given a segmentaiton.
    % We will load a precomupted_MCMC_sample and extract the segmentation 
    % from it. We will then measure non-rigidity exhibited by this
    % segmentation.
    
    % Reads from disk (./Data and ./Results) meshes preprocessed into 
    % point clouds as well as precomputed MCMC samples.
    
    % I/P: DATASET - 'Centaur' or 'Horse'
    % USAGE : demo_analysis('Centaur')
    
    close all
    
    datapath = './Data';
    savedir = './Results';
    [~,left_tris,pcf] = load_data(datapath,DATASET);
    load(fullfile(datapath,sprintf('tri2opposite_Tosca%s.mat',DATASET)));

    
    % load precomputed labels
    load(fullfile(savedir,sprintf('precomputed_converged_MCMC_sample_%s.mat',DATASET)));
    labels = converged_sampler_state.t;
   
    % Measure non rigidity
    [~,average_tri_err] = measure_nonrigidity(pcf,labels);
    
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % visualize non rigidity %%%%%%%%%%%%%%%
    labels= zeros(1,size(left_tris,2)*2);
    labels(left_tris) = average_tri_err;
    % optionally enforce symmetry by projecting the right mesh segentation on the left.
    right_tris = setdiff(1:size(left_tris,1)*2,left_tris);
    labels(right_tris) = average_tri_err;
    
    if(strcmp(DATASET,'Centaur'))
        fnames =  dir(fullfile(datapath,'meshes','cen*.mat'));
    else
        fnames =  dir(fullfile(datapath,'meshes','hor*.mat'));
    end
    data = load(fullfile(datapath,'meshes',fnames(1).name));
    surface = data.surface;
    vrts(:,1) = surface.X;
    vrts(:,2) = surface.Y;
    vrts(:,3) = surface.Z;
    faces = surface.TRIV;
   
   
   figure(1);t =trisurf(faces,vrts(:,1),vrts(:,3),vrts(:,2),labels,'EdgeColor','none');title('Non rigidity visualization. Red/Yellow = high nonrigidity');
   set(t,'FaceLighting','phong','AmbientStrength',0.5);
   
   light('Position',[5 0 5],'Style','infinite', 'Color', [.25 .25 .25]);
   light('Position',[0 5 5],'Style','infinite', 'Color', [.25 .25 .25]);
   light('Position',[-5 -5 5],'Style','infinite', 'Color', [.25 .25 .25]);
   
   axis equal;
   
   axis off;
   set(gcf, 'Color', [1 1 1]);
   view(-180,-80);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  non_rigidity = labels;
  save(fullfile(savedir,'PerPart_Nonrigidity.mat'),'non_rigidity');