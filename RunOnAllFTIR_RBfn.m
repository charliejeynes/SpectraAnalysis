
clc
clear
close all

%%

FTIRmouseData_RBfn = struct('filename',{}, 'kmeancluster',{}, 'Ctransposed',{}, 'dataOfInterest', {}); 

% filename = '/Users/jcgj201/Documents/MATLAB/FTIR data/17-10-2016/25-02-16-WT/25-02-16-wt.dmt'; 

% filename = dir('/Users/jcgj201/Documents/MATLAB/FTIR data/17-10-2016'); 

filename = dir('/Users/jcgj201/Documents/MATLAB/FTIR data/allresults/*.dmt');

% filename = dir('/Users/jcgj201/Documents/MATLAB/FTIR data/17-10-2016/02-03-16-TG/*.dmt');

%%
for i = 1:length(filename)
    
    img = [filename(i).folder , '/',  filename(i).name ]; 

    [kmeancluster, wavenumbers3T, Ctransposed] = kMeansROI_RB(img); 
    
    FTIRmouseData_RBfn(i).kmeancluster = kmeancluster; 
    FTIRmouseData_RBfn(i).filename = filename(i).name; 
    FTIRmouseData_RBfn(i).Ctransposed = Ctransposed; 
    
%     h = warndlg('press to continue', 'do you want the next image?'); 
%     drawnow; 
%     waitfor(h); 
    
    
    
    input = inputdlg({'dataOfInterest', 'UserID', 'age', 'sex: M/F'}, 'Participent Details', [1 40; 1 40; 1 40; 1 40], ...
    {'', '', '', ''}); 


    session  = input{1}; 
    FTIRmouseData_RBfn(i).dataOfInterest = input{1} ; 
%     = struc(kmeancluster,{}, wavenumbers3T,{}, Ctransposed,{},  filename,{});
    
    
end

%% save the structure file


save('17-10-2016.mat', 'FTIRmouseData_RBfn')
%%

% f  = figure(gcf); 
% savefig(f, '5') ;


%%
% clustercolumn = []; 
for k = 1:length(FTIRmouseData_RBfn)
   
    
    clustercolumn(k) = str2num(FTIRmouseData_RBfn(k).dataOfInterest); 
       
end

%%

for k = 1:length(FTIRmouseData_RBfn)

    for j = 1:length(clustercolumn)
        f = clustercolumn(j); 
        spectraForPCA(:, j) = FTIRmouseData_RBfn(k).Ctransposed(:, f); 
    end
       
end
%% plot picked out spectra from the kmeans Clustering

figure()
plot(wavenumbers3T,spectraForPCA)
% xlabel('1st Principal Component')
% ylabel('2nd Principal Component')
%% rubber band in a for loop on all the spectrum THIS NEEDS TO BE SORTED

% wavenumbersT = wavenumbers'; length(data3_pixels_wn_amide)
% tic
% spectraPCAamide = spectraForPCA(1:417, :); 
% amideLength = size(spectraPCAamide); 
% 
% for k = 1:amideLength(2)
%     spectraPCAamideT = spectraPCAamide'; 
%     p = spectraPCAamideT; 
%      
%     n = wavenumbers3T(1:417); 
%     n = n'; 
%  
%     [z,a,it,ord,s,fct] = backcor(n,p, 2, 1, 'atq');
%     z = z'; 
%     spectraPCAamideT_RB(k, :) = z;
% 
% %     subtracted(k, :) = data3_pixels_wn_amide(k, :) - data3_pixels_wn_amide_RB(k, :); 
%     spectraPCAamideT(k, :) = spectraPCAamideT(k, :) - spectraPCAamideT_RB(k, :);
% %     subtracted = y' - z; 
% end
toc
%% Plots the region
% figure, 
%region = wavenumbers3(1:417); 
figure, plot(n, spectraPCAamideT(:, :))
set(gca,'xdir','reverse')
figure, plot(n', spectraPCAamide(:, :))
set(gca,'xdir','reverse')

%% run the PCA analysis on Ctransposed

% allKmeanSpectrum = FTIRmouseData(2).Ctransposed; 
% 
% [coeff, score, latent]  = pca(allKmeanSpectrum); 
% 
% figure()
% plot(score(:,1),score(:,2),'+')
% xlabel('1st Principal Component')
% ylabel('2nd Principal Component')



