function data = getData(folder, theFiles, condition)
%  ------------------------------------------------------------------------
%  Author: Siddhartha Dhiman
%  E-Mail: sdhiman@buffalo.edu
%  Created: 02/16/2018 using Matlab R2017b
%  ------------------------------------------------------------------------
%  Function to generate cell counts, deviation and process length
%  into the matrix structure 'data'. This matrix contains the
%  following information:
%       FileName:           Name of files processed
%       ROI:                ROI index number
%       TotalCells:         Total events/cells detected
%       LongProc:           Count of cells with long processes
%       ShortProc:          Count of cells with short processes
%       NoProc:             Count of cells with no processes
%       Deviation:          Row vector with devaiation of long processes
%       MajorAxisLength:    Length of all events
%   -----------------------------------------------------------------------
%   Input Arguments:
%       folder:     Path to folder containing images
%       theFiles:   Structural matrix with directory information. Can be
%                   constructed using 'dir(filePattern)', where
%                   'filePattern' is a file filtering criteria.
%       condition:  1 for Control
%                   2 for Schizophrenia
%   -----------------------------------------------------------------------
%   Output:
%       data: Structural matrix containing information useful for
%                generating biological networks.
%   -----------------------------------------------------------------------
%   Feel free to drop me an email on questions and concerns
%   -----------------------------------------------------------------------

%% Clearing Workspace
clc;
clear all;
close all;
warning off;

%% Define Variables
R_SE1 = 50;     %Tophat structuring element radius

%% Structuring Elements
SE1 = strel('disk', R_SE1);     % For Tophat

%% Open 3D Matrix
Image = load('MATLAB_Converted/Multiphoton/MS2.mat');
Image = Image.catI;
[X Y Z] = size(Image);
resX = 0.82; % X-Resolution
resY = 0.82; % Y-Resolution
resZ = 2.0;  % Z-resolution

%% Apply Filters
filteredI = zeros(size(Image));
I = zeros(size(Image,1),size(Image,2));
for i = 1:size(Image,3)
    I = Image(:,:,i);
    tophatI = imtophat(I,SE1);
    eqI = adapthisteq(tophatI);
    filteredI = cat(3,tophatI,filteredI);
end

Filt = [3 3 3];
FilteredI3d = medfilt3(filteredI, Filt);
FilteredI3d = imboxfilt3(FilteredI3d, Filt);

nSuperPx = ceil(sqrt(M*N*P))*3;
[L,N] = superpixels3(filteredI,nSuperPx,'Method','slic');

imSize = size(FilteredI3d);
imPlusBoundaries = zeros(imSize(1),imSize(2),3,imSize(3),'uint8');
for plane = 1:imSize(3)
    BW = boundarymask(L(:,:,plane));
    % Create an RGB representation of this plane with boundary shown
    % in cyan.
    imPlusBoundaries(:,:,:,plane) = imoverlay(FilteredI3d(:,:, plane),BW);
    pixelBW(:,:,plane) = BW;
end

implay(imPlusBoundaries,5)

pixelIdxList = label2idx(L);
stats = regionprops3(L,FilteredI3d,'MeanIntensity','SurfaceArea',...
    'Volume','EigenValues','VoxelIdxList');

shape = stats.SurfaceArea./stats.Volume;

eigen = stats.EigenValues;
for i = 1:length(eigen)
    e = eigen{i};
    for j = 1:length(e);
        eigenValues(i,j) = e(j);
    end
end

for superpixel = 1:N
        memberPixelIdx = pixelIdxList{superpixel};
        
% 
% for i = 1:size(eigenValues,1)
%     for j = 1:size(eigenValues,2)
%         orientation(i,j) = eigenValues(i,j)/numel(stats.VoxelIdxList{i});
%     end
% end

intensity = double(stats.MeanIntensity);

figure, scatter(shape,intensity);
figure, scatter3(orientation(:,1),orientation(:,2),orientation(:,3))

featureMatrix = horzcat(shape,intensity,eigenValues);


%% Define a Reference 3D Coordinate System for Image
coordRef = imref3d([X Y Z],resX,resY,resZ);

%% Plotting 3D Figure
[x,y,z] = meshgrid(1:X,1:Y,1:Z);

%% Slice
[a b c] = find()
V = patch(isosurface(pixelBW,0.5),'FaceColor','r','EdgeColor','none')
slice(V,X/2,Y/2,Z/2);
grid on, shading interp, colormap gray;
