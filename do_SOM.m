function [sM, sMap, bmus] = do_SOM(data, nx_som, ny_som)

% NOTE: Based off of scipts provided by Valentina Radic.
%
% Syntax: [sM, sMap, bmus] = do_SOM(data, nx_som, ny_som)
% 
% sM: This is the self-organized map.  Note: sM.codebook contains the cluster patterns.
% sMap: This is the original map (before training) -- used in some optional figures below.  
% bmus: This is a vector containing the number of the best matching unit for each row in 'data'.  The pattern of bmu 'x' is sM.codebook(x,:).
% data: matrix of doubles; rows: stations; columns: 'time series'.
% nx_som: number of columns in SOM.
% ny_som: number of rows in SOM.

en = nx_som*ny_som;

msize=[ny_som nx_som];
% performing linear initialization of nodes (i.e. the nodes are distibuted in the space of
% first two eignevectors from PCA
display('initialization')
data(isnan(data))=0; %SOM can't deal with NaN's -- set them to zero
sMap=som_lininit(data,'msize',msize,'hexa','sheet');

% training SOM
display('training')
[sM,sT] = som_batchtrain(sMap,data,'ep','hexa','sheet','radius',[2 1],'trainlen',200);

% calulating error
[q,t]=som_quality(sM,data);

% calulating hits (frequencies) of occurences of each pattern, for each seasn
hi=som_hits(sM,data);
hi=100*hi/sum(hi);

bmus=som_bmus(sM,data);

%% Plots you might use:

% % plot inital pattterns
% figure;
% for i=1:en
%     subplot(ny_som,nx_som,i);
%     plot(1:365,sMap.codebook(index(i),:));
%     xlabel('time');
%     title(['initial node ' num2str(index(i))])
% end
%
% %% plot colored SOM (with actual hexagonal nodes)
% %  and distance among neighboring nodes
% U = som_umat(sM);
% Um = U(1:2:size(U,1),1:2:size(U,2));
% C = som_colorcode(sM,'rgb2');
% 
% figure('Renderer','Painters')
% som_cplane(sM,C,1-Um(:)/max(Um(:)));
% title('Color coding + distance matrix')
% 
% % plot colored SOM (with actual hexagonal nodes)
% % and visualize the SOM patterns inside each node
% figure('Renderer','Painters')
% som_cplane(sM,C);
% hold on
% som_plotplane(sM,sM.codebook)
% title('SOM with patterns')
% set(gca,'xticklabel',[])
% hold off