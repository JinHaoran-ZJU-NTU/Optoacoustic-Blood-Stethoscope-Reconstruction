%% Data load

close all
load Data.mat

%% Parameter

fs = 1/kgrid.dt;        % sampling rate 
px = dx;                % scan interval (In this 3D imaging demo, dx = dy ) 

focus_length = 0;       % unfocused transducer

density = 2;            % NUFFT interpolation density

fmin = 1e6;             % transducer bandwidth
fmax = 8e6;


% Layered Material Parameter


c = [1750,1450];        % Sound Speed; ( c = [c1,c2,c3,c4...,cn]);
layer = 2.1e-3;         % Layer thickness; (layer = [thick1,thick2,...thick(n-1)]; The final layer thickness is set as infinite.


disp = 0;               % Time offset; suggest to set it as 0


%% RawData
rfdata = permute(sensor_data(:,:,1:310),[3,1,2]);   


%% Reconstruction



tic;
migRF2 = PS_3D_NUFFT(rfdata,fs,px,disp,layer,c,fmin,fmax,focus_length,density);
toc;


%% Display
display_z = 1:105;      % set a reasonable display range
figure,imagesc(squeeze(max(abs(migRF2(display_z,:,:)))));       





