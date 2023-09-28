

% define subject directory;
Subdir = '/media/charleslynch/storage_1/SIMD/data/SIMD04';
str = strsplit(Subdir,'/');
Subject = str{end};

% define the midthickness surfaces;
MidthickSurfs{1} = [Subdir '/anat/T1w/fsaverage_LR32k/' Subject '.L.midthickness.32k_fs_LR.surf.gii'];
MidthickSurfs{2} = [Subdir '/anat/T1w/fsaverage_LR32k/' Subject '.R.midthickness.32k_fs_LR.surf.gii'];

% load the resting-state data;
C = ft_read_cifti_mod([Subdir '/func/rest/ConcatenatedCiftis/Rest_OCME+MEICA+MGTR_Concatenated+SubcortRegression+SpatialSmoothing2.55.dtseries.nii']);
load([Subdir '/func/rest/ConcatenatedCiftis/FD.mat']);
C.data = single(C.data(:,FD < 0.3)); % scrub & convert to single type;

% read in the infomap communities / subnetworks
N = ft_read_cifti_mod([Subdir '/pfm2/Bipartite_PhysicalCommunities+SpatialFiltering.dtseries.nii']);
N.data = N.data(:,6); % extract the column of interest

% load priors;
load('priors.mat');

% run the network identification algorithm;
identify_networks(C,N,MidthickSurfs,[],Priors,...
'Bipartite_PhysicalCommunities+AlgorithmicLabeling',[pwd '/SIMD04']);


