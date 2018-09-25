%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sam Anderson         %%
%% April 30, 2018       %%
%% Data --> KDE --> SOM %%
%% Clear everything
clear
clc

%% Input data
filename = '/Users/jeromemayaud/Documents/University/BritishColumbia/Modelling/MachineLearning/bin_midpoints_surrey.csv';
bins_all = cell2mat(table2cell(readtable(filename)));
bins_one = bins_all(1,:);

filename = '/Users/jeromemayaud/Documents/University/BritishColumbia/Modelling/MachineLearning/income_surrey_2016_nansremoved.csv';
data_mat = readtable(filename);
data_mat = data_mat(:,2:end);
data_mat = cell2mat(table2cell(data_mat));

%remove rows in data that are all zero
tot = sum(data_mat'); 
rowsWithData = find(tot~=0);
data_mat = data_mat(rowsWithData,:);

%convert to input for kde -- currently, 'data' is the output of a
%histogram; we want to reproduce the input to this histogram, and then use that to
%make a kde; do this by making a vector for each row in data that has the
%values of the bin centre the number of times as is its frequency in data

for kk = 1:length(data_mat(:,1)) %for each grid cell
    
    clear input2kde
    start = 1; %note: start:fin will be the indices in input2kde for each bin's occurrances 
    input2kde = [];
    
    for jj = 1:length(data_mat(1,:)) %for each 'frequency' of the bin centre (ie: each histogram value for this grid cell)
        
        if data_mat(kk,jj)>0 %if there is a non-zero frequency for this bin
            fin = start+data_mat(kk,jj)-1; 
            input2kde(start:fin) = bins_all(kk,jj)*ones(1,fin-start+1);
            start = fin+1;
        end
        
    end
    
    [kde_vec,kde_income(kk,:)] = ksdensity(input2kde);
    data_kde(kk,:) = kde_vec;  
    
end

numSpatialCells = size(data_mat,1); %number of grid cells in space, where each grid cell has its own data
numrows = 20; %number of spatial rows you want to display in plots
numcols = floor(numSpatialCells/numrows); %number of spatial columns

%% Linearly interpolate all data onto the same x-axis
income = linspace(min(min(kde_income)),max(max(kde_income)),100); 
data_mat_interp = zeros(length(data_mat(:,1)),length(income));
for jj = 1:length(data_mat(:,1)) %for each spatial cell with data
    data_mat_interp(jj,:) = interp1(kde_income(jj,:),data_kde(jj,:),income,'linear');
end

data_mat_interp(isnan(data_mat_interp)) = 0; %the interpolation sets values outside the original domain to NaN; set those values to zeros
data_mat = data_mat_interp;

%% Normalize the data so we're comparing the shapes of the distribution, and not the values themselves

for kk = 1:length(data_mat(:,1)) %standardize each row
    data_norm(kk,:) = (data_mat(kk,:) - mean(data_mat(kk,:))) / std(data_mat(kk,:));
end

%% If you're comparing a second date (e.g. 2022) with one whose SOM clusters you've already calculated (e.g. 2016):
%% Need to have the 'sM' variable from the already-calculated year (e.g. 2016)
%% ONLY RUN THIS SECTION IF YOU ALREADY HAVE SOM CLUSTERS FROM PREVIOUS YEAR - OTHERWISE SKIP

cluster_pdfs = sM.codebook; %This is the PDFs of the clusters from the baseline year (e.g. 2016)
BMU_table = zeros((length(data_mat(:,1))), 1);

for yy = 1:length(data_mat(:,1))
    comparison_table = zeros((length(cluster_pdfs(:,1))), 1);
    current_pdf = data_norm(yy,:);
    for zz = 1:length(cluster_pdfs(:,1))
        RMSE = sqrt(mean((current_pdf - cluster_pdfs(zz,:)).^2));
        comparison_table(zz,1) = RMSE;
    end
    %Find smallest RMSE of all the clusters, and assign the corresponding BMU to the cell
    [value, index] = min(comparison_table(:)); 
    [row, col] = ind2sub(size(comparison_table), index);
    BMU_table(yy,1) = row; %The row number is effectively the BMU
end

%% Make a self-organizing map of the data

[N,nx,ny] = getNumberOfPatterns(data_norm); %use heirarchical clustering method to get number of patterns in the data
[sM, sMap, bmus] = do_SOM(data_norm,nx,ny); %make the SOM

figure %visualize SOM patterns
for kk = 1:N
    subplot(ny,nx,kk)
    plot(income,sM.codebook(kk,:))
    bmu_freq(kk) = length(find(bmus==kk))/length(bmus);
    title(sprintf('Node #%d \n Frequency = %0.1f%%',kk,bmu_freq(kk)*100))
end

figure %visualize best matching units in space
imagesc(reshape(bmus,numrows,numcols))
colorbar


%% ALTERNATIVELY: do all the above but use anomalies, rather than the raw data
%% THE INTERPOLATION OF THIS SECTION ISN'T CORRECT - NEED TO DO KDE PROPERLY (see above)

%next, linearly interpolate all data onto the same x-axis
x_all = linspace(min(min(x)),max(max(x)),1000); 
data_mat_interp = zeros(length(data_mat(:,1)),length(x_all));
for jj = 1:length(data_mat(:,1)) %for each spatial cell with data
    data_mat_interp(jj,:) = interp1(x(jj,:),data_mat(jj,:),x_all,'linear');
end
data_mat_interp(isnan(data_mat_interp)) = 0; %the interpolation sets values outside the original domain to NaN; set those values to zeros
data_mat = data_mat_interp;

%now, compute average signal and remove it from each row
avg_dist = mean(data_mat);
data_anomaly = zeros(size(data_mat));
for kk = 1:length(data_mat(:,1))
    data_anomaly(kk,:) = data_mat(kk,:) - avg_dist;
end

%% now, normalize the data so we're comparing the shapes of the distribution, and not the values themselves

for kk = 1:length(data_mat(:,1))   
    data_anomaly_norm(kk,:) = ( data_anomaly(kk,:) - mean(data_anomaly(kk,:)) ) / std(data_anomaly(kk,:));   
end

% figure %visualize normalized data
% for kk = 1:length(data_mat_interp(:,1))
%     yyaxis left, plot(x_all,data_anomaly_norm(kk,:))
%     yyaxis right, plot(x_all,data_anomaly(kk,:))
%     legend('Normalized','Interpolated');
%     pause(0.1)
% end

%% next, make a self-organizing map of the data

[N,nx,ny] = getNumberOfPatterns(data_anomaly_norm); %use heirarchical clustering method to get number of patterns in the data
[sM, sMap, bmus] = do_SOM(data_anomaly_norm,nx,ny); 

figure %visualize SOM patterns
for kk = 1:N
    subplot(ny,nx,kk)
    plot(x_all,sM.codebook(kk,:))
    bmu_freq(kk) = length(find(bmus==kk))/length(bmus);
    title(sprintf('Node #%d \n Frequency = %0.1f%%',kk,bmu_freq(kk)*100))
end

figure %visualize best matching units in space
imagesc(reshape(bmus,numrows,numcols))
colorbar



