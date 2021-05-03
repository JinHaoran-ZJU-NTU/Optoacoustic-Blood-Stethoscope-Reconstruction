%% Data load

close all
load Data.mat

%% Parameter

% Scan and Sampling Parameters

fs = 1/kgrid.dt;        % sampling rate 
px = dx;                % scan interval (In this 3D imaging demo, dx = dy ) 

focus_length = 0;       % unfocused transducer

fmin = 1e6;             % transducer bandwidth
fmax = 8e6;


% Layered Material Parameters
c = [1750,1450];        % Sound Speed; ( c = [c1,c2,c3,c4...,cn]);
layer = 2.1e-3;         % Layer thickness; (layer = [thick1,thick2,...thick(n-1)]; The final layer thickness is set as infinite.


disp = 0;               % Time offset; If the sampling does not start from t =0, then set the disp = time offset (s).


% Reconstruction Parameters
density = 2;            % NUFFT interpolation density (available density = 1, 2, 4, ...)

%% RawData
rfdata = permute(sensor_data(:,:,1:310),[3,1,2]);   % 3D axis ------- (t, x, y)


%% Reconstruction

tic;
migRF2 = PS_3D_NUFFT(rfdata,fs,px,disp,layer,c,fmin,fmax,focus_length,density);
toc;


%% Display
display_t = 100:200;    % set a reasonable display range 
display_z = 1:105;      % set a reasonable display range

figure(1),imagesc(squeeze(max(abs(rfdata(display_t,:,:))))); title('raw data');
figure(2),imagesc(squeeze(max(abs(migRF2(display_z,:,:))))); title('PS-NUFFT');






