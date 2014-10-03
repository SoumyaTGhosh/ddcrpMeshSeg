function [ data ] = package_data(X,Y,numSubj)
% package data for later processing 
% numSubj - number of subjects to be processed (For the datasets included
% with this package, this is 1. If you have scans of multiple subjects set
% this parameter appropriately. It could take any value between 1 and
% length(data.X)).
   
% massage data
data.X = {X{1:numSubj}};


for ii = 1:numSubj
    % transpose Y.
    tY = [];
    for jj = 1:size(Y{ii},3)
        tY(:,:,jj) = Y{ii}(:,:,jj)';
    end
    data.Y{ii} = tY;
    
end


    % number of distinct mesh subjects per partition
    % each subject has her own reference mesh
    data.num_bodies = length(data.X);



end

