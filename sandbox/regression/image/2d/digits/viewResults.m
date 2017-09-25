%close all;
addpath '../utilities/matlab/'

I = imread('data/I2.png');
CP =  load('output/Regression_ControlPoints.txt');
MOM = load('output/Regression_InitialMomenta.txt');

figure;
imagesc(I); colormap gray; colorbar;
hold on
quiver(CP(:,1),CP(:,2),MOM(:,1),MOM(:,2),0,'LineWidth',3);

