%% Clear everything
clear
clc

%% Input data

filename = '/Users/jeromemayaud/Documents/University/BritishColumbia/Modelling/MachineLearning/bin_midpoints_surrey.csv';
bins_all = cell2mat(table2cell(readtable(filename)));
bins_one = bins_all(1,:);


filename = '/Users/jeromemayaud/Documents/University/BritishColumbia/Modelling/MachineLearning/income_surrey_2016_nansremoved.csv';
data = readtable(filename);
data = data(:,2:end);
data = cell2mat(table2cell(data));

%remove rows in data that are all zero
tot = sum(data'); 
rowsWithData = find(tot~=0);
data = data(rowsWithData,:);


%convert to input for kde -- currently, 'data' is the output of a
%histogram; we want to reproduce the input to this histogram, and then use that to
%make a kde; do this by making a vector for each row in data that has the
%values of the bin centre the number of times as is its frequency in data
for kk = 1:length(data(:,1)) %for each grid cell
    
    clear input2kde
    start = 1; %note: start:fin will be the indices in input2kde for each bin's occurrances 
    input2kde = [];
    
    for jj = 1:length(data(1,:)) %for each 'frequency' of the bin centre (ie: each histogram value for this grid cell)
        
        if data(kk,jj)>0 %if there is a non-zero frequency for this bin
            fin = start+data(kk,jj)-1; 
            input2kde(start:fin) = bins_all(kk,jj)*ones(1,fin-start+1);
            start = fin+1;
        end
        
    end
    
    [kde_vec,kde_income(kk,:)] = ksdensity(input2kde);
    data_kde(kk,:) = kde_vec;  
    
end

%% interpolate data
income = linspace(min(min(kde_income)),max(max(kde_income)),100); 
data_interp = zeros(length(data(:,1)),length(income));
for jj = 1:length(data(:,1)) %for each spatial cell with data
    data_interp(jj,:) = interp1(kde_income(jj,:),data_kde(jj,:),income,'linear');
end

data_interp(isnan(data_interp)) = 0; %the interpolation sets values outside the original domain to NaN; set those values to zeros
data = data_interp;

%% normalize data
for kk = 1:length(data(:,1)) %standardize each row
    data_norm(kk,:) = (data(kk,:) - mean(data(kk,:))) / std(data(kk,:));
end

%% If you're comparing a second date (e.g. 2022) with one whose PCA clusters you've already calculated (e.g. 2016):
%% Need to have the 'PCs, eigvecs and eigvals' variables from the already-calculated year (e.g. 2016)
%% ONLY RUN THIS SECTION IF YOU ALREADY HAVE SOM CLUSTERS FROM PREVIOUS YEAR - OTHERWISE SKIP

k = 8; %This is number of clusters you want
meas = [PCs(:,1),PCs(:,2),PCs(:,3)];
Z = linkage(meas,'ward','euclidean');
c = cluster(Z,'maxclust',k);

% Calculate Cluster Centres
for ii = 1:k
    [ind dummy] = find(c == ii);
    PCs_avg(ii,1) = mean(PCs(ind,1));
    PCs_avg(ii,2) = mean(PCs(ind,2));
    PCs_avg(ii,3) = mean(PCs(ind,3));
end

% Use cluster centres to find most representative distributions
for jj = 1:k
    repPat(jj,:) = PCs_avg(jj,1)*eigvecs(:,1)' + PCs_avg(jj,2)*eigvecs(:,2)' + PCs_avg(jj,3)*eigvecs(:,3)';
end

data_mat = data; 

PCA_table = zeros((length(data_mat(:,1))), 1);

for yy = 1:length(data_mat(:,1))
    comparison_table = zeros((length(repPat(:,1))), 1);
    current_pdf = data_norm(yy,:);
    for zz = 1:length(repPat(:,1))
        RMSE = sqrt(mean((current_pdf - repPat(zz,:)).^2));
        comparison_table(zz,1) = RMSE;
    end
    %Find smallest RMSE of all the clusters, and assign the corresponding BMU to the cell
    [value, index] = min(comparison_table(:)); 
    [row, col] = ind2sub(size(comparison_table), index);
    PCA_table(yy,1) = row; %The row number is effectively the PCA
end


%% now, do principal component analysis to see what modes capture the variance in this dataset
 
[PCs, eigvecs, eigvals] = pca(data_norm'); 
% columns of eigvecs are the linearly independant modes of income
% columns of PCs are the weightings used to reconstruct the data
% eigvals is related to the variance explained by each mode

%% visualize fraction of variance explained by each mode

totalVar = sum(eigvals); %total variance in the data
fracVar = eigvals/totalVar*100; %percentage of total variance accounted for by each mode

figure
scatter(1:length(fracVar),fracVar)
xlabel('Mode')
ylabel('% Variance Explained')
title('Fraction of Variance Explained by Each PCA Mode')


%% visualize modes and PCs

numModes = 3; %number of modes to visualize
figure
for kk = 1:numModes
    
    subplot(numModes,2,2*kk-1)
    plot(income,eigvecs(:,kk))
    xlabel('Income')
    titlestr = sprintf('Eigenvector #%d',2*kk-1);
    title(titlestr)
    
    subplot(numModes,2,2*kk)
    plot(PCs(:,kk))
    xlabel('Grid Cell')
    titlestr = sprintf('PCs');
    title(titlestr)
    
end

%% reconstruct data to see how it matches up
for kk = 1:length(data(:,1))
    data_rec(kk,:) = eigvecs(:,1)*PCs(kk,1) + eigvecs(:,2)*PCs(kk,2) + eigvecs(:,3)*PCs(kk,3); %each row = reconstructed data_norm using first few modes
    rho(kk) = corr(data_rec(kk,:)',data_norm(kk,:)');
end

% figure %visualize correlation bw real data and reconstructed data
% plot(1:length(rho),rho)
% xlabel('Grid Cell')
% ylabel('Correlation')
% title('Correlation b/w Real Data and PCA-Reconstructed Data')

figure %visualize real data and reconstructed data
for kk = 1:12
    subplot(3,4,kk)
    hold on
    plot(income, data_norm(kk,:),'k')
    plot(income, data_rec(kk,:),'b')
    legend('Original Data','Reconstructed Data')
end

%% cluster in PC-space to find representative patterns

% make dendrogram
meas = [PCs(:,1),PCs(:,2),PCs(:,3)];
Z = linkage(meas,'ward','euclidean');

figure;
dendrogram(Z)
title('Dendrogram');

%from dendrogram, choose four clusters for each temp and precip

%% cluster PCs for temp and precip

k = 8;
c = cluster(Z,'maxclust',k);

% plot clustered data, each cluster in different color
figure;
scatter3(PCs(:,1),PCs(:,2),PCs(:,3),30,c,'filled')
xlabel('PC1','fontsize',14);
ylabel('PC2','fontsize',14);
zlabel('PC3','fontsize',14);
title('Clusters','fontsize',18);

% Calculate Cluster Centres
for ii = 1:k
    [ind dummy] = find(c == ii);
    PCs_avg(ii,1) = mean(PCs(ind,1));
    PCs_avg(ii,2) = mean(PCs(ind,2));
    PCs_avg(ii,3) = mean(PCs(ind,3));
end

%% use cluster centres to reconstruct 'most representative' distributions

for jj = 1:k
    repPat(jj,:) = PCs_avg(jj,1)*eigvecs(:,1)' + PCs_avg(jj,2)*eigvecs(:,2)' + PCs_avg(jj,3)*eigvecs(:,3)';
end

figure %visualize representative patterns
for kk = 1:k
    subplot(4,2,kk),plot(income,repPat(kk,:))
    xlabel('Income')
    titlestr = sprintf('Mean Pattern of Cluster #%d',kk);
    title(titlestr)
end

