%%
clear 
clc
%%


% [wavenumbers, data, width, height, filename, acqdate] = readvarianmosaic_v4_1(filename, keepme); 
% images =  dir('/Users/jcgj201/Documents/MATLAB/FTIR data/17-10-2016/02-03-16-TG/*.dmt'); 
% for i=1:
image1 =  '/Users/jcgj201/Documents/MATLAB/FTIR data/allresults/02-03-16-TG.dmt'; 
image2 =  '/Users/jcgj201/Documents/MATLAB/FTIR data/allresults/25-02-16-WT.dmt'; 

[wavenumbers, data, width, height, filename, acqdate] = readvarianmosaic_v4_1(image1);
[wavenumbers2, data2, width2, height2, filename2, acqdate2] = readvarianmosaic_v4_1(image2);

%% data1

figure('Name',filename2,'NumberTitle','off');
imagesc(sum(data,3));
title('02-03-16-TG')

%% data2

figure('Name',filename2,'NumberTitle','off');
imagesc(sum(data2,3));

%% select a region data2
y = 110:120; % height
x =  300:350; %width

% dataROI = squeeze(data(y,x,1200:1300)); 
dataROI2 = squeeze(data2(y,x,:));

%% data 2

figure('Name',filename2,'NumberTitle','off');
imagesc(sum(dataROI2,3));
colorbar

%% data 2
ROIdata2reshaped = reshape(dataROI2, 561, 1506);
mean2 = mean(ROIdata2reshaped); 
figure, plot(wavenumbers2, mean2)

%% data 1
% 
% figure, 
% plot(wavenumbers, squeeze(dataROI(1, 1, :))); 

%% select a region data 1
y = 80:90; % height
x =  200:230; %width

% dataROI = squeeze(data(y,x,1200:1300)); 
dataROI = squeeze(data(y,x,:));

%% data1

figure('Name',filename,'NumberTitle','off');
imagesc(sum(dataROI,3));
colorbar

%% data 1
ROIdatareshaped = reshape(dataROI, 341, 1506);
mean = mean(ROIdatareshaped); 
figure, plot(wavenumbers, mean)

%% normalise with 2nd derivative 

ndDerivMean = gradient(mean); 
figure, plot(wavenumbers, ndDerivMean)


%% plot both

figure, plot(wavenumbers, mean, wavenumbers, mean2)
legend('WT', 'TG')

%% display just the 2 parts of the spectrum

amide = 1:417; 
lipid = 920:1120; 
combined = [amide, lipid]; 

%%

figure; plot(wavenumbers([1:417,920:1120]), mean([1:417,920:1120]), wavenumbers([1:417,920:1120]), mean2([1:417,920:1120]))

%% normalise with 1st derivative 

oneDerivMean = gradient(mean); 
oneDerivMean2 = gradient(mean2);

figure,
plot(wavenumbers(amide), oneDerivMean(amide), wavenumbers(amide), oneDerivMean2(amide))
hold on
plot(wavenumbers(lipid), oneDerivMean(lipid), wavenumbers(lipid), oneDerivMean2(lipid))
hold off

%% normalise with 2nd derivative 

twoDerivMean = gradient(gradient(mean)); 
twoDerivMean2 = gradient(gradient(mean2));

figure,
plot(wavenumbers(amide), twoDerivMean(amide), wavenumbers(amide), twoDerivMean2(amide))
hold on
plot(wavenumbers(lipid), twoDerivMean(lipid), wavenumbers(lipid), twoDerivMean2(lipid))
hold off



%% split the data into 2 regions data1

data3_pixels_wn_amide = mean(:, 1:417);

data3_pixels_wn_CH = mean(:, 920:1120); 

%% rubber band baseline on the amide peak in a for loop on all the spectrum in the ROI data1

% wavenumbersT = wavenumbers'; length(data3_pixels_wn_amide)
tic
disp ('Rubber band baseline Amide')
amideLength = size(data3_pixels_wn_amide); 

for k = 1:amideLength(1)
    p = data3_pixels_wn_amide(k, :); 
    n = wavenumbers(1:417); 
%  
    [z,a,it,ord,s,fct] = backcor(n,p, 2, 0.01, 'atq');
    z = z'; 
    data3_pixels_wn_amide_RB(k, :) = z;

%     subtracted(k, :) = data3_pixels_wn_amide(k, :) - data3_pixels_wn_amide_RB(k, :); 
    data3_pixels_wn_amide(k, :) = data3_pixels_wn_amide(k, :) - data3_pixels_wn_amide_RB(k, :);
%     subtracted = y' - z; 
end
toc

%% split the data2 into 2 regions data 2

data3_pixels_wn_amide2 = mean2(:, 1:417);

data3_pixels_wn_CH2 = mean2(:, 920:1120); 

%% rubber band baseline on the amide peak in a for loop on all the spectrum in the ROI data 2

% wavenumbersT = wavenumbers'; length(data3_pixels_wn_amide)
tic
disp ('Rubber band baseline Amide')
amideLength2 = size(data3_pixels_wn_amide2); 

for k = 1:amideLength2(1)
    p = data3_pixels_wn_amide2(k, :); 
    n = wavenumbers(1:417); 
%  
    [z,a,it,ord,s,fct] = backcor(n,p, 2, 0.01, 'atq');
    z = z'; 
    data3_pixels_wn_amide_RB2(k, :) = z;

%     subtracted(k, :) = data3_pixels_wn_amide(k, :) - data3_pixels_wn_amide_RB(k, :); 
    data3_pixels_wn_amide2(k, :) = data3_pixels_wn_amide2(k, :) - data3_pixels_wn_amide_RB2(k, :);
%     subtracted = y' - z; 
end
toc


%% plot RB baseline amide

difference = data3_pixels_wn_amide - data3_pixels_wn_amide2; 
 
 
figure(), plot(wavenumbers(1:417), data3_pixels_wn_amide, wavenumbers(1:417), data3_pixels_wn_amide2, wavenumbers(1:417), difference);
legend('02-03-16-TG', '25-02-16-WT', 'difference')
title('Background subtracted comparison between a rTG4510 and WT mouse')
ylabel('wavenumber')
xlabel('absorbance')
set(gca,'xdir','reverse')
set(findall(gcf,'-property','FontSize'),'FontSize',6)

%% normalise 

Y1 = msnorm(wavenumbers(1:417)', data3_pixels_wn_amide', 'Max', 1);
Y2 = msnorm(wavenumbers(1:417)', data3_pixels_wn_amide2', 'Max', 1);

figure, plot(wavenumbers(1:417), Y1', wavenumbers(1:417), Y2')

%% random number 

randNumbers = rand(1, 417);  

%% run PCA on 2 spectrum

% pca_data = [data3_pixels_wn_amide', data3_pixels_wn_amide2', randNumbers'];

pca_data = [data3_pixels_wn_amide', data3_pixels_wn_amide2']; 

% pca_data = [data3_pixels_wn_amide;  data3_pixels_wn_amide2];

%% plot PCA  RB normal data


% pca(v(1:3,7:10)
[wcoeff,score,latent,tsquared,explained] = pca(pca_data); 

figure()
plot(score(:,1),score(:,2),'+')
xlabel('1st Principal Component')
ylabel('2nd Principal Component')

%% PCA 1st derivative

oneDeriv = [oneDerivMean2', oneDerivMean2']; 
[wcoeff,score,latent,tsquared,explained] = pca(oneDeriv);

figure()
plot(score(:,1),score(:,2),'+')
xlabel('1st Principal Component')
ylabel('2nd Principal Component')

%% PCA 2nd derivative

twoDeriv = [twoDerivMean2', twoDerivMean2']; 
[wcoeff,score,latent,tsquared,explained] = pca(twoDeriv);

figure()
plot(score(10:20-100:30,1),score(:,2),'+')
xlabel('1st Principal Component')
ylabel('2nd Principal Component')
%%



%% get most relevant numbers

importantWaveNumbers = 335:345; 
WN417 = wavenumbers(1:417); 
WN417(:, importantWaveNumbers) 

%%
figure, 
plot(wavenumbers(1:417), score(:,3)); 




















%% reshape the data so that the Co2 can be cut out

sizeData = size(data); % puts the elements of the array into a []

allHeight = sizeData(1);

allWidth = sizeData3(2);

allPixels = prod(sizeData3(1:2)); 

spectrum = sizeData3(3); 

data3toBeCut = reshape(data3, prod(sizeData3(1:2)), sizeData3(3)); % MUST BE PIXELS first for some reason YES CHECKED

figure,  plot(wavenumbers3, data3toBeCut(60000,:))

%%
for i=1:length(data3toBeCut)
    
%     data3toBeCut(679:728, i) = data3toBeCut(678, i); 
    data3toBeCut( i, 679:728) = data3toBeCut( i, 678);
    %data3toBeCut( i, 679:728) = 0.1;
end

%%
figure,  plot(wavenumbers3, data3toBeCut(60000,:))
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


%% image the  data3 with sum

window_title = filename3;
figure('Name',window_title,'NumberTitle','off');
imagesc(sum(data3ROI,3));axis image;axis off;
[pathstr, name, ext] = fileparts(filename3); 
titlefilename = fullfile(pathstr,name);
title(titlefilename, 'interpreter', 'none');
