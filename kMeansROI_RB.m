function [kmeancluster, wavenumbers3T, Ctransposed, data3, filename] = kMeansROI_RB(filename)

image3 =  filename; 
% image3 = '/Users/jcgj201/Documents/MATLAB/FTIR data/17-10-2016/25-02-16-WT/25-02-16-wt.dmt'; 
[wavenumbers3, data3, width3, height3, filename3, acqdate3] = readvarianmosaic_v4_1(image3);

[pathstr, name, ext] = fileparts(filename3);

%% reshape the data so that the Co2 can be cut out

sizeData3 = size(data3); % puts the elements of the array into a []

allHeight = sizeData3(1);

allWidth = sizeData3(2);

allPixels = prod(sizeData3(1:2)); 

spectrum = sizeData3(3); 

data3toBeCut = reshape(data3, prod(sizeData3(1:2)), sizeData3(3)); % MUST BE PIXELS first for some reason YES CHECKED

% figure,  plot(wavenumbers3, data3toBeCut(60000,:))

%%
for i=1:length(data3toBeCut)
    
%     data3toBeCut(679:728, i) = data3toBeCut(678, i); 
    data3toBeCut( i, 679:728) = data3toBeCut( i, 678);
    %data3toBeCut( i, 679:728) = 0.1;
end

%%
% figure,  plot(wavenumbers3, data3toBeCut(60000,:))
%%

data3 = reshape(data3toBeCut, allHeight, allWidth, 1506);

%% select a region
intoQuarterH = allHeight / 4;
intoQuarterW= allWidth / 4;
y = allHeight / 4 : intoQuarterH + 100; % height
x = intoQuarterW : intoQuarterW + 200 ; %width

data3ROI = squeeze(data3(y,x,:)); 
%%

pixels = length(y) * length(x); 
% reshapedData = reshape(data3ROI, 1506, pixels); %wrong way does something
% funny with the spectrum
reshapedData = reshape(data3ROI, pixels, 1506);  % the right way round




%% reshape the data so that the Co2 can be cut out SPLIT

% sizeData3 = size(data3); % puts the elements of the array into a []
% 
% %prod(sizeData3(1:2)); 
% 
% data3_pixels_wn = reshape(data3, prod(sizeData3(1:2)), sizeData3(3)); % MUST BE PIXELS first for some reason YES CHECKED
% 
% figure,  plot(wavenumbers3, data3_pixels_wn(34000,:))
% 
% test = data3_pixels_wn(34000,:); 

%% split the data into 2 regions 

data3_pixels_wn_amide = reshapedData(:, 1:417);
% middlebit = ; 
data3_pixels_wn_CH = reshapedData(:, 920:1236); 
% endbit = ; 

%% check the cut out

figure,  plot(wavenumbers3(1:417), data3_pixels_wn_amide(1, :)); 

%% rubber band normalise on 1 spectrum

% % wavenumbersT = wavenumbers';
% 
% p = data3_pixels_wn_amide(1, :); 
% n = wavenumbers3(1:417); 
% 
% [z,a,it,ord,s,fct] = backcor(n,p,2, 0.01, 'atq'); 
% %[z,a,it,ord,s,fct] = backcor(n,y); 
% subtracted = p' - z; 
% 
% %% check the RB on single spectrum
% 
% figure,  plot(wavenumbers3(1:417), data3_pixels_wn_amide(1, :)); 
% 
% subtractedT = subtracted'; 
% 
% figure,  plot(wavenumbers3(1:417), p(1, :)); 

%% rubber band in a for loop on all the spectrum

% wavenumbersT = wavenumbers'; length(data3_pixels_wn_amide)
tic
disp ('Rubber band baseline')
amideLength = size(data3_pixels_wn_amide); 

for k = 1:amideLength(1)
    p = data3_pixels_wn_amide(k, :); 
    n = wavenumbers3(1:417); 
%  
    [z,a,it,ord,s,fct] = backcor(n,p, 2, 0.01, 'atq');
    z = z'; 
    data3_pixels_wn_amide_RB(k, :) = z;

%     subtracted(k, :) = data3_pixels_wn_amide(k, :) - data3_pixels_wn_amide_RB(k, :); 
    data3_pixels_wn_amide(k, :) = data3_pixels_wn_amide(k, :) - data3_pixels_wn_amide_RB(k, :);
%     subtracted = y' - z; 
end
toc

%% Plots the region
% figure, 
% %region = wavenumbers3(1:417); 
% plot(n, subtracted(1:200, :))
% set(gca,'xdir','reverse')
% figure, plot(n, data3_pixels_wn_amide(1:200, :))
% set(gca,'xdir','reverse')



%% try the normalisation routine
% tic
% 
% subLength = size(subtracted); 
% 
% for i=1:subLength(1)
% 
%     subtracted(i, :) = subtracted(i, :) / subtracted(i, 341) ; 
%     
% end 
% toc



%% image the data3 with sum
% 
% window_title = filename3;
% figure('Name',window_title,'NumberTitle','off');
% imagesc(sum(data3ROI,3));axis image;axis off;
% [pathstr, name, ext] = fileparts(filename3); 
% titlefilename = fullfile(pathstr,name);
% title(titlefilename, 'interpreter', 'none');

%% try the normalisation routine
tic
for i=1:length(reshapedData) 
    
%     disp(i)
    data3_pixels_wn_amide(i, :) = data3_pixels_wn_amide(i, :) / data3_pixels_wn_amide(i, 341) ; 
    
end 
toc

%% check the spectrum after normalisation to amide and Rubber band
figure, plot(wavenumbers3(1:417), data3_pixels_wn_amide(1:10, :))

%% check the look of the amideI normalised data - - works well

amideInorm_data3ROI = reshape(data3_pixels_wn_amide, length(y), length(x), 417); 
% 
figure, imagesc(sum(amideInorm_data3ROI,3));axis image;axis off;
colorbar

%% transpose the data for correct line up for k means
Y = data3_pixels_wn_amide; 
wavenumbers3T = wavenumbers3';
% figure, plot(wavenumbers3T, Y(:, 20000)) 

%%
tic
[dataROI_idx_6,C] = kmeans(Y,6); 
toc

%% for data 3 ROI

kmeancluster = reshape(dataROI_idx_6, length(y), length(x));

%%
figure, imagesc(kmeancluster);
title(filename)
colorbar

%% 
Ctransposed = C'; 

figure, plot(wavenumbers3T(1:417, :),  Ctransposed);
legend('show'); 

% f = figure('Units','normalized', ...
%         'Position',[0.25 0.25 0.5 0.5], 'Name', filename, 'NumberTitle', 'off');
%         ax = axes('Position',[0 0 0.8 1]);

       
        
f = figure;       
subplot(2, 2, 1); 
imagesc(sum(data3,3));axis image;axis on;
colorbar 

subplot(2, 2, 2); 
% title(filename);
imagesc(sum(data3ROI,3));axis image;axis on;
colorbar 

subplot(2, 2, 3); 
imagesc(sum(kmeancluster,3));axis image;axis on;
colorbar 

subplot(2, 2, 4); 
plot(1:417, Ctransposed);
set(gca,'xdir','reverse')
legend('show'); 


set(gcf, 'units','normalized','outerposition',[0 0 1 1]); % amke it full screen

end

%%
% mkdir('/Users/jcgj201/Documents/MATLAB/FTIR data/allresults')
% 
% name = [name, '.tif'];
% 
% saveas(f, name)

%end

