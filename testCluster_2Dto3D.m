%% Clear Workspace
clc;
clear all;
close all;

%% Variables
R_SE1 = 20;             % Top-Hat structuring element radius (counting)
SE1 = strel('disk', R_SE1);     % For Top-Hat (Counting)

filteredI = [];
for i = 1:size(catI,3)
    I = catI(:,:,i);
    I = mat2gray(I);
    tophatI = imtophat(I, SE1);
    eqI = adapthisteq(I,'numTiles',[8 8],'nBins',128);
    filteredI = cat(3,eqI,filteredI);
end

%% Simple Linear Iteration Clustering (SLIC)
Filt = [5 5 5];
FilteredI3d = medfilt3(filteredI, Filt);
FilteredI3d = imboxfilt3(FilteredI3d, Filt);
[L,N] = superpixels3(FilteredI3d,128,'Method','slic');

imSize = size(catI);
imPlusBoundaries = zeros(imSize(1),imSize(2),3,imSize(3),'uint8');
for slice = 1:imSize(3)
  BW = boundarymask(L(:,:,slice));
  % Create an RGB representation of this plane with boundary shown
  % in cyan.
  imPlusBoundaries(:,:,:,slice) = imoverlay(catI(:,:,slice),BW,'cyan');
end

implay(imPlusBoundaries,5)