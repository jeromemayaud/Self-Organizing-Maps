function [N,nx,ny] = getNumberOfPatterns(data)

% NOTE: Based off of scipts provided by Valentina Radic.  
%
% Syntax: [N,nx,ny] = getNumberOfPatterns(data)
%
% this function determines the optimal number of patterns (ie: number of
% nodes to use in SOM).  User can determine Nx and Ny from N.  This function
% uses heirarchical clustering method.
%
% N = double; total number of patterns
% nx = double; number of rows to input to SOM 
% ny = double; number of columns to input to SOM
% data = matrix of doubles; rows: stations; columns: time series
%
% Further Reading: V. Radic, K. Unglert, and M. Jellinek, "Machine Learning Method for Pattern Recognition in
% Volcano Seismic Spectra"

sD = som_data_struct(data);
sM = som_make(sD); %this is the SOM
ny_som = sM.topol.msize(1);
nx_som = sM.topol.msize(2);
en=ny_som*nx_som;

Bmus = som_bmus(sM,sD); %best matching units of each node

myhits = som_hits(sM,sD); %number of nodes in each BMU
[inbNaN dummy]=find(myhits == 0);


% plot colored SOM (with actual hexagonal nodes)
%  and distance among neigbouring nodes
U = som_umat(sM);
Um = U(1:2:size(U,1),1:2:size(U,2));
C = som_colorcode(sM,'rgb2');

C0=C;
C0(inbNaN,1)=1;
C0(inbNaN,2)=1;
C0(inbNaN,3)=1;

% % plot colored SOM (with actual hexagonal nodes)
% % and SOM patterns
% figure('Renderer','Painters')
% som_cplane(sM,C0);
% hold on
% som_plotplane(sM,sM.codebook)
% title('SOM with pattern: Large Size from SOM')
% set(gca,'xticklabel',[])
% hold off

%do PCA on the SOM nodes
[eigvec, pcs, eigvals] = pca(sM.codebook);
pc_topo = sqrt(pcs(:,1).^2 + pcs(:,2).^2);

pc_topo_mat = reshape(pc_topo,ny_som,nx_som);

% figure %visualize PC1^2 + PC2^2 field
% s = surf(pc_topo_mat);
% s.EdgeColor = 'none';
% title('PC1^2 + PC2^2')
% 
% figure('Renderer','Painters')
% som_cplane(sM,pc_topo/max(pc_topo));
% title('PC1^2 + PC2^2')
% set(gca,'xticklabel',[])
% colorbar

localmax_mat = imregionalmax(pc_topo_mat); %1s where local maxima exist, 0s elsewhere
localmax_lin = reshape(localmax_mat,length(pc_topo(:,1)),length(pc_topo(1,:))); %convert to a line for plotting

globalmin = min(pc_topo_mat(:)); %global minimum value in pc topography space
globalmin_mat = zeros(length(pc_topo_mat(:,1)),length(pc_topo_mat(1,:)));
globalmin_mat(find(pc_topo_mat==globalmin)) = 1;
globalmin_lin = reshape(globalmin_mat,length(pc_topo(:,1)),length(pc_topo(1,:))); %convert to a line for plotting

maxmin_mat = globalmin_mat + localmax_mat;
maxmin_lin = globalmin_lin + localmax_lin;

% figure('Renderer','Painters') %visualize locations of maxima/minima
% som_cplane(sM,maxmin_lin/2);
% title('Locations of Local Maxima and Global Minima')
% set(gca,'xticklabel',[])
% colorbar

maxmin_patterns = sM.codebook(find(maxmin_lin==1),:); %patterns at the global min/local max locations

%now, find which of the max/min patterns most closely resemble each of the other nodes

for kk = 1:length(sM.codebook(:,1)) %for each node
    
    for jj = 1:length(maxmin_patterns(:,1)) %look at each of the max/min patterns    
        RMSE(jj,1) = sqrt(sum((sM.codebook(kk,:) - maxmin_patterns(jj,:)).^2)/length(sM.codebook(kk,:))); %RMSE between node kk and the jjth max/min mode       
    end
    
    bestind(kk,1) = find(RMSE == min(RMSE)); %
    
end

% figure('Renderer','Painters') %visualize the patterns and their cluster
% som_cplane(sM,bestind);
% hold on
% som_plotplane(sM,sM.codebook)
% title('Optimized SOM with patterns')
% set(gca,'xticklabel',[])
% colorbar

%%

N = length(maxmin_patterns(:,1)); %'ideal' number patterns: now -- what should Nx and Ny be? Concept: keep aspect ratio most similar
if N==factor(N) %if N is prime
    N = N+1; %add one more pattern?
end
factors = factor(N);

best_aspect_ratio = ny_som/nx_som; %try to match this aspect ratio with grid size

index = 1;
for kk = 1:length(factors)
    
    combinations = combnk(factors,kk); %ie: fCkk (from factors, choose kk)
    
    for jj = 1:length(combinations(:,1)) %for each combination
              
        Nx_possibilities(index) = prod(combinations(jj,:));
        Ny_possibilities(index) = N/Nx_possibilities(index);
        index = index+1;
        
    end
    
end

aspect_ratio_possibilities = Ny_possibilities./Nx_possibilities;
RMSE_aspect_ratio = (best_aspect_ratio - aspect_ratio_possibilities).^2 / length(aspect_ratio_possibilities);
[err ind] = min(RMSE_aspect_ratio);
nx = Nx_possibilities(ind);
ny = Ny_possibilities(ind);

