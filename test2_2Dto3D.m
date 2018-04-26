%% Clear Workspace
clc;
clear all;
close all;

%% Variables
R_SE1 = 10;             % Top-Hat structuring element radius (counting)

%% Structuring Elements
SE1 = strel('disk', R_SE1);     % For Top-Hat (Counting)
SE2 = strel('disk', R_SE2);     % For morphological operations
SE3 = strel(ones(sq,sq));       % For marker cleaning

filteredI = [];
for a = 1:size(catI,3)
    I = catI(:,:,a);
%% Contrast Enhancment of gray image using CLAHE
I = mat2gray(I);
tophatI = imtophat(I, SE1);
eqI = adapthisteq(I,'numTiles',[8 8],'nBins',128);

%% Background Exclusion
% Apply Average Filter
h = fspecial('average', [7 7]);
JF = imfilter(eqI, h);

% Take the difference between the gray image and Average Filter
Z = imsubtract(JF, eqI);

%% Threshold using the IsoData Method
level=graythresh(Z); % this is our threshold level

%% Convert to Binary
BW = imbinarize(Z, level);

%% Remove small pixels
BW2 = bwareaopen(BW, 100);

%% Open
BW3 = imclose(BW2,strel('disk',1));

filteredI = cat(3,BW3,filteredI);
end

filteredI = smooth3(filteredI);
isosurface(filteredI,1/2);

[faces,vertices]= isosurface(filteredI,1/2);

figure;
patch('Faces',faces,'Vertices',vertices)
axis vis3d
grid on
colormap jet