%% Data load

close all
load Data.mat


%% Parameter
rfdata = permute(sensor_data(:,:,1:310),[3,1,2]);   %
fs = 1/kgrid.dt;        % sampling rate
pitch = dx;             % scanning interval

focus_length = 0;       % unfocused transducer
density = 2;            % NUFFT interpolation density
fmin = 1e6;             % transducer bandwidth
fmax = 8e6;

c = [1750,1450];        % Sound Speed;
layer = 2.1e-3;         % First layer thickness (Nz/2*dz); determined by the simulation parameter: medium.sound_speed
disp = 0;
% Layered
tic;
migRF2 = Omega_K_3D_NUFFT(rfdata,fs,pitch,disp,layer,c,fmin,fmax,focus_length,density);
toc;

figure,imagesc(squeeze(max(abs(migRF2(1:105,:,:)))));





