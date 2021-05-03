function  migRF = PS_3D_NUFFT(rfdata,fs,pitch,disp,layer,c,fmin,fmax,focus_length,density)

%%
[L,M,N]= size(rfdata);      %L is time axis, M,N are x-y axis

%% Experiment Condition
% General
fs  = fs;             % sample rate: param.fs (MHz) --> fs (Hz)
dt  = 1/fs;
dx = pitch;         % axis X pitch (mm) --> dx (m)
dy = pitch;         % axis Y pitch  (mm)--> dy (m)  //dx = dy
% Material
ERMv = c;
th = layer;             % layertk (mm) --> th (m)

% Transducer 
% fc      = param.transducer.fc*1e6;      % center frequency (MHz) --> fc (Hz)
fmin  = fmin;  % bandwidth lower limit --> fmin (Hz)
fmax  = fmax; % bandwidth upper limit --> fmax (Hz)
% Record Time Delay
disp = disp;            % param.delayt (us) --> delayt (s)

%% FFT
newM = 2^nextpow2(M);
newN = 2^nextpow2(N);
newL = 2^nextpow2(L);

P0wk = fftshift( fftn(rfdata,[newL,newM,newN]) );
kx = 2*pi*(-newM/2 : newM/2-1)/dx/newM;             
ky = 2*pi*(-newN/2 : newN/2-1)/dy/newN;             
w  = 2*pi*(-newL/2 : newL/2-1)/dt/newL;  

df = fs/newL;
Nw1=floor( fmin/df )+(newL/2+1);
Nw2=ceil( fmax/df )+(newL/2+1);
P0wk = P0wk(Nw1:Nw2,:,:);
w = w(Nw1:Nw2);
%% Wave field extrapolation
focus_time = focus_length/ERMv(1);  
Trunc_P0wk = P0wk;                                          
[gkx,gw,gky] = meshgrid(kx,w,ky);                                 

gkz2 = gw.^2/ERMv(1)^2 - gkx.^2-gky.^2; mask = gkz2 >= 0;
% All elements of kz for which (w.^2/ERMv(2)^2 - kx.^2) <0 should therefore be set to zero
gkz = sign(gw).*sqrt(gkz2.*mask);                            
Trunc_P0wk = Trunc_P0wk.*(mask~=0);                         
Trunc_P0wk = Trunc_P0wk.*exp(1i.*gkz*(th-focus_length));
                                                            

E_P0wk = Trunc_P0wk.*exp(-1i*gw*(disp-focus_time));                  

%% NUFFT migration


f = gw;
fkz = ERMv(2)*sign(f).*sqrt(gkx.^2+gky.^2 + f.^2/ERMv(2)^2); 
                                       
%fftRF = interp3(gkx,f,gky,E_P0wk,gkx,fkz,gky,'linear',0);
fftRF = InterpNUFFT(f,kx,ky,E_P0wk,fkz,gkx,gky,density);


% Jacobian (optional)
%kz = (-newN/2:newN/2-1)'/ERMv/fs/newN;
gkz = f/ERMv(2);     
fftRF = fftRF.*gkz./(fkz+eps);
%MASK = (gkx.^2 + gky.^2 + (gkz./(fkz+eps)).^2 )< ((f/ERMv(2)).^2);
%fftRF = bsxfun(@times,fftRF,gkz)./(fkz+eps); 

%% IFFT 
fftdata = zeros(newL,newM,newN);
fftdata(Nw1:Nw2,:,:) = (fftRF);
fftRF = fftdata;

% IFFT & Migrated RF
migRF = ifftn(ifftshift(fftRF),'nonsymmetric');
migRF = migRF(1:L,1:M,1:N);
end