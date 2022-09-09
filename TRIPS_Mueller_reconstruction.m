clc
clear
close all
%% include the toolbox path 

addpath(genpath('./sspsoct_toolbox'));

%% load data
% load calibration data
load('CalStru_cor.mat');

% a single frame to generate the image
% struct name: frame
load('S123.mat');


%% reconstruct tomography, retardance and optic axis map

PSsize = 500; % the PS image size

PS = PSTriple_Manager_v2(CalStru);

PS.load_Stokes(frame.S1,frame.S2,frame.S3);
PS.checkSymmetrizationCorrection();     
PS.applySymmetrizationCorrection();
PS.applyPMDCorrection();
            
PS.flip_Mueller();

% A qualitative retardation map that applied polar decomposition to get
% retardation matrix mR
% ret = vector length of rotation matrix (mR(z)*inv(mR(z-1)))
% please look into the function to see details
ret = PS.get_local_retardance();

%% flatten the image along the RPE of the retina 

fm = FlattenManagerForHuman();
I = frame.I;

% TomographyManager_TM: a class with some functions to resize the images
tom = TomographyManager_TM(CalStru);            
Ips = tom.SNRmap(tom.avg_imresize(I,PSsize)); 

% find how to flat
fm.surface_flatten_setup(flipud(frame.dopu),Ips);

% apply the flattening
Ipsf = fm.flatten_image(Ips);
PS = fm.flatten_Mueller(PS);

%% remove system and cornea birefringence

% remove square root of the entire system mueller matrix alone b-scan
PS.remove_linear_birefringence_from_surface();

% local optic axis reconstruction
[retFinal,phi] = PS.local_optic_axis();
phi = fm.reverse_flatten_image(phi);




    
%% Image presentation
iw = ImageWriter('demo','_demo');
imex = iw.write_intensity_image(Ips,[0 30]);
figure(10)
imshow(imex(:,20:end-20,:)) % noisy edges is due to a small misaglignemnt of galvo and modulator

imex = iw.write_color_ret_image(Ips,ret,[0 20],[0.02 1.8],0);
figure(11)
imshow(imex(:,20:end-20,:))

%imex = iw.write_ret_oa_image(((ret).*(Ips>1.5)).^(0.5),phi,[0.15 0.75],hsv(256));
imex = iw.write_ret_oa_image(((ret).*(Ips>2)),phi,[0.05 0.6],hsv(256));
figure(12)
imshow(imex(:,20:end-20,:))