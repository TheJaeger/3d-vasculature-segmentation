%% Clear Workspace
clc;
clear all;
close all;

%% Variables
R_SE1 = 50;             % Top-Hat structuring element radius (counting)
R_SE2 = 2;              % Morphological opening structuring element radius
sq = 5;                 % Size of 'sq x sq' square matrix for marker cleaning
open_AC = 50;            % Area for noise removal using 'bwareaopen'
cannyThresh = 0.30;     % Threshold of Canny filter
cannySigma = 1;         % Sigma of Canny filter

%% Structuring Elements
SE1 = strel('disk', R_SE1);     % For Top-Hat (Counting)
SE2 = strel('disk', R_SE2);     % For morphological operations
SE3 = strel(ones(sq,sq));       % For marker cleaning

%% Looping per Layer

sizeI = size(catI);
I = zeros(sizeI(1),sizeI(2));
skel3D = [];
filteredI = [];
for a = 1:size(catI,3)
    I = catI(:,:,a);
    %% Top-Hat Operation
tophatI = imtophat(I, SE1);
eqI = adapthisteq(tophatI);

%% Morphological Processing
Io = imopen(eqI, SE2);
erodeI = imerode(tophatI, SE2);
Iobr = imreconstruct(erodeI, tophatI);
Ioc = imclose(Io, SE2);
Iobrd = imdilate(Iobr, SE2);
Iobrcbr = imreconstruct(imcomplement(Iobrd), imcomplement(Iobr));
Iobrcbr = imcomplement(Iobrcbr);

%% Threshold Image
thresh = 0.20;
bw1 = imbinarize(Iobrcbr,thresh);
bw2 = bwmorph(bw1,'fill');
bw3 = bwmorph(bw2,'bridge');
bw4 = imclose(bw3,strel('disk',5));
bw5 = bwareaopen(bw4,open_AC);

% %% Vessel Segmentation
% vesselBoundary = bwboundaries(bw5);
% vesselMask=false(size(bw1));
% for i = 1:length(vesselBoundary)
%   for j = 1:length(vesselBoundary{i})
%       ind = vesselBoundary{i}(j,:);
%       vesselMask(ind(1),ind(2))=1;
%   end
% end
% vesselMask = bwmorph(vesselMask,'bridge');
% vesselMask = bwmorph(vesselMask,'clean');
% vesselMask = bwmorph(vesselMask,'fill');
% vesselMask = imfill(vesselMask,'holes');
% vesselMask = imclose(vesselMask, strel('disk',2));

%% Skeletonize VesselMask
vesselMask = bwmorph(bw5,'thin',Inf);
vesselmask = bwmorph(bw5,'spur',Inf);

filteredI = cat(3,bw5,filteredI);
skel3D = cat(3,vesselMask,skel3D);
end

filteredI = smooth3(filteredI);
interpI = interp3(filteredI,'spline');

[X, Y, Z] = meshgrid((0:size(filteredI,1)-1)*0.82,...
    (0:size(filteredI,2)-1)*0.82,...
    (0:size(filteredI,3)-1)*2.0 );
[faces,vertices]= isosurface(X,Y,Z,filteredI);

figure;
hold on
patch(isosurface(X,Y,Z,filteredI));
patch(isosurface(X,Y,Z,filteredI));
isonormals(X,Y,Z,filteredI,p)
isonormals(X,Y,Z,skel3D,q);
hold off
p.FaceColor = 'red';
p.EdgeColor = 'none';
q.FaceColor = 'black';
q.Edgecolor = 'none';
daspect([1 1 1])
view(3); 
axis tight
camlight 
lighting gouraud

figure;hold on
colormap([1 0 0;0 0 1]) %red and blue
surf(X,Y,Z,filteredI); %first color (red)
surf(X,Y,Z,skel3D); %first color (red)