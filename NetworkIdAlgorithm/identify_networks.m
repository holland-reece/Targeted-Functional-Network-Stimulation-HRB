function identify_networks(C,Ic,MidthickSurfs,Col,Priors,OutFile,OutDir)
% cjl; cjl2007@med.cornell.edu;
rng(44); % for reproducibility.

% make output directory ;
if ~exist(OutDir,'dir')
    mkdir(OutDir);
end

% read in the 
% resting-state data
if ischar(C)
    C = ft_read_cifti_mod(C);
end

% count the number of cortical vertices; note: this should be 59412;
nCorticalVertices = nnz(C.brainstructure==1) + nnz(C.brainstructure==2);

% fc matrix;
m = paircorr_mod(C.data(1:nCorticalVertices,:)');
m(eye(size(m,1))==1) = 0; % remove the diagonal;
m(isnan(m)) = 0; % remove nans

% read in 
% infomap
% communities 
if ischar(Ic)
    Ic = ft_read_cifti_mod(Ic); % Ic == infomap communities
end

% extract
% graph density
% of interest;
if ~isempty(Col)
    Ic.data = Ic.data(:,Col);
end

O = Ic; % preallocate output
O.data = zeros(size(Ic.data));

% unique communities;
uCi = unique(nonzeros(Ic.data));

% preallocate functional connectivity of each community
uCi_FC = zeros(nCorticalVertices,length(uCi)); %

% sweep each of the
% unique communities
for i = 1:length(uCi)
    
    % calculate this communit functional connectivity profile
    uCi_FC(:,i) = nanmean(m(:,Ic.data(1:nCorticalVertices)==uCi(i)),2);
    
end

% calculate spatial
% similarity of FC with templates
uCi_rho = corr(uCi_FC,Priors.FC);

% preallocate a measure of spatial
% overlap (each community & probability maps)
uCi_prob = zeros(length(uCi),size(Priors.Spatial,2)); %

% sweep each of the
% unique communities
for i = 1:length(uCi)
    
    % loop through the networks;
    for ii = 1:size(Priors.Spatial,2)
        
        % average probability of community "i" belonging to functional network "ii"
        uCi_prob(i,ii) = mean(Priors.Spatial(Ic.data(1:nCorticalVertices)==uCi(i),ii));
        
    end
    
end

% preallocate;
rd = zeros(1,length(uCi)); % rd = relative difference

% preallocate some variables;
Community = zeros(length(uCi),1);
Networks = cell(length(uCi),1);
R = zeros(length(uCi),1);
G = zeros(length(uCi),1);
B = zeros(length(uCi),1);
Confidence = zeros(length(uCi),1);
Alternatives = cell(length(uCi),1);

% sweep through
% the communities;
for i = 1:length(uCi)
    
    % sort from best match to worst match; 
    [x,y] = sort(uCi_rho(i,:) .* uCi_prob(i,:),'Descend'); % x = similarity values, y = index of networks
    
    % calculate relative difference (rd) between similarity with best network x(1) vs runner-up x(2)
    rd(i) = abs((x(1)-x(2))/x(2)); % rd == measure of  "confidence"; higher values == better match
    Ci_best_match = y(1);
   
    % log all
    % the info;
    Community(i) = i;
    R(i) = Priors.NetworkColors(Ci_best_match,1);
    G(i) = Priors.NetworkColors(Ci_best_match,2);
    B(i) = Priors.NetworkColors(Ci_best_match,3);
    Networks{i} = Priors.NetworkLabels{Ci_best_match};
    Confidence(i) = rd(i); % confidence
    Alternatives{i} = Priors.NetworkLabels{y(2)};

    % log the network identity;
    O.data(Ic.data(:,1)==uCi(i),1) = Ci_best_match;
    
end

% write out .xls sheet summarizing the algorithms decisions;
T = table(Community,Networks,R,G,B,Confidence,Alternatives); %
writetable(T,[OutDir '/NetworkLabels.xls']);

% write out the temporary cifti file;
ft_write_cifti_mod([OutDir '/Tmp.dtseries.nii'],O);

% write out the first network;
system(['echo ' char(Priors.NetworkLabels{1}) ' > ' OutDir '/LabelListFile.txt ']);
system(['echo 1 ' num2str(round(Priors.NetworkColors(1,1)*255)) ' ' num2str(round(Priors.NetworkColors(1,2)*255)) ' ' num2str(round(Priors.NetworkColors(1,3)*255)) ' 255 >> ' OutDir '/LabelListFile.txt ']);

% sweep through the networks;
for i = 2:length(Priors.NetworkLabels)
    
    system(['echo ' char(Priors.NetworkLabels{i}) ' >> ' OutDir '/LabelListFile.txt ']);
    system(['echo ' num2str(i) ' ' num2str(round(Priors.NetworkColors(i,1)*255)) ' ' num2str(round(Priors.NetworkColors(i,2)*255)) ' ' num2str(round(Priors.NetworkColors(i,3)*255)) ' 255 >> ' OutDir '/LabelListFile.txt ']);
    
end

% make dense label file + network borders;
system(['wb_command -cifti-label-import ' OutDir '/Tmp.dtseries.nii ' OutDir '/LabelListFile.txt ' OutDir '/' OutFile '.dlabel.nii -discard-others']);
system(['wb_command -cifti-label-to-border ' OutDir '/' OutFile '.dlabel.nii -border ' MidthickSurfs{1} ' ' OutDir '/' OutFile '.L.border']); % LH
system(['wb_command -cifti-label-to-border ' OutDir '/' OutFile '.dlabel.nii -border ' MidthickSurfs{2} ' ' OutDir '/' OutFile '.R.border']); % RH

% remove some intermediate files;
system(['rm ' OutDir '/Tmp.dtseries.nii']);
system(['rm ' OutDir '/LabelListFile.txt']);

O = Ic; % preallocate output;
O.data = zeros(size(O.data,1),size(uCi_FC,2)); % blank slate
O.data(1:nCorticalVertices,:) = uCi_FC; % log community FC profiles
ft_write_cifti_mod([OutDir '/CommunityFC.dtseries.nii'],O);

end

% version #1 ; not used
% % sweep through
% % the communities;
% for i = 1:length(uCi)
%     
%     % sort from best match to worst match  (based on similarity of FC profile)
%     [x,y] = sort(uCi_rho(i,:),'Descend'); % x = similarity values, y = index of networks
%     
%     % calculate relative difference (rd) between similarity with best network x(1) vs runner-up x(2)
%     rd(i) = abs((x(1)-x(2))/x(2)); % rd == measure of  "confidence"; higher values == better match
%     
%     % note: if the initial classification is ambiguous (either because there is a
%     % close "runner-up" identity or if the match is poor in an abolsute sense...
%     if rd(i) < rd_threshold || x(1) < absolute_threshold %
%         
%         % which of the top two possibilities
%         % have the best spatial overlap
%         % with the  network probability map?
%         [~,z] = max(uCi_prob(i,y(1:2)));
%         Ci_best_match = y(z); % opt for runner-up, if needed.
%         
%         % log the new runner-up; gives user option to easily manually
%         % review after and override algorithm's decision if needed;
%         Alternatives{i} = Priors.NetworkLabels{y(1:2~=z)};
%         
%     else
%         Ci_best_match = y(1);
%     end
%     
%     % log all the info;
%     Community(i) = i;
%     R(i) = Priors.NetworkColors(Ci_best_match,1);
%     G(i) = Priors.NetworkColors(Ci_best_match,2);
%     B(i) = Priors.NetworkColors(Ci_best_match,3);
%     Networks{i} = Priors.NetworkLabels{Ci_best_match};
%     Confidence(i) = rd(i); % confidence
%     
%     % log the network identity;
%     O.data(Ic.data(:,1)==uCi(i),1) = Ci_best_match;
%     
% end

