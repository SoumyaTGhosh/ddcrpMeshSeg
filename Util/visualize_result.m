function [labels,t] = visualize_result(f,c_MCMC_param,sa_m,left_tris,dataset)
meshpath = './Data';

param = c_MCMC_param;
%%% randpermute the order of labels to make colormap more meaningful.
randt = zeros(1,length(param.t(sa_m,:)));
mc1 = randperm(length(unique(param.t(sa_m,:))));
for ii = 1:length(unique(param.t(sa_m,:)))
    randt(param.t(sa_m,:)==ii) = mc1(ii);
end;

labels= zeros(1,size(left_tris,2)*2);
labels(left_tris) = randt;
% optionally enforce symmetry by projecting the right mesh segentation on the left.
%  right_tris = setdiff(1:size(left_tris,1)*2,left_tris);
%  labels(right_tris) = randt;

if(strcmp(dataset,'Centaur'))
    fnames =  dir(fullfile(meshpath,'meshes','cen*.mat'));
else
    fnames =  dir(fullfile(meshpath,'meshes','hor*.mat'));
end
data = load(fullfile(meshpath,'meshes',fnames(1).name));
surface = data.surface;
vrts(:,1) = surface.X;
vrts(:,2) = surface.Y;
vrts(:,3) = surface.Z;
faces = surface.TRIV;

figure(f);t =trisurf(faces,vrts(:,1),vrts(:,3),vrts(:,2),labels,'EdgeColor','none');%title(sprintf('%d',sa_m));
set(t,'FaceLighting','phong','AmbientStrength',0.5);

light('Position',[5 0 5],'Style','infinite', 'Color', [.25 .25 .25]);
light('Position',[0 5 5],'Style','infinite', 'Color', [.25 .25 .25]);
light('Position',[-5 -5 5],'Style','infinite', 'Color', [.25 .25 .25]);

axis equal;

axis off;
set(gcf, 'Color', [1 1 1]);
view(-180,-80);

end