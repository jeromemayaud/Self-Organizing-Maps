clear
clc

filename = '/Users/jeromemayaud/Documents/University/BritishColumbia/Modelling/MachineLearning/SOM_income_2016.csv';
bins_all = cell2mat(table2cell(readtable(filename)));
bins_one = bins_all(1,:);
income = bins_one;

filename = '/Users/jeromemayaud/Documents/University/BritishColumbia/Modelling/MachineLearning/SOM_cluster_distributions_2016.csv';
data = readtable(filename);
data = data(:,1:end);
data = cell2mat(table2cell(data));

%remove rows in data that are all zero
tot = sum(data'); 
rowsWithData = find(tot~=0);
data = data(rowsWithData,:);

%% now, do principal component analysis to see what modes capture the variance in this dataset
 
[PCs, eigvecs, eigvals] = pca(data'); 

%% visualize fraction of variance explained by each mode

totalVar = sum(eigvals); %total variance in the data
fracVar = eigvals/totalVar*100; %percentage of total variance accounted for by each mode

figure
scatter(1:length(fracVar),fracVar)
xlabel('Mode')
ylabel('% Variance Explained')
title('Fraction of Variance Explained by Each PCA Mode')

%NOTE: It is unsurprising that the first mode captures so much variance.
%To explain ~80% of the variance is massive, and I think this comes from
%the fact that all data should look like a distribution similar to this.
%To get modes that explain more minute features, we could look at income
%anomalies.

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