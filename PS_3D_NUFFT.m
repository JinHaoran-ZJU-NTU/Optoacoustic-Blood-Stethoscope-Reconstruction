function  migRF = PS_3D_NUFFT(rfdata,fs,pitch,disp,layer,c,fmin,fmax,density)

%%
[L,M,N]= size(rfdata);      %L is time axis, M,N are x-y axis

%% Experiment Condition
% General
fs  = fs;                   % sample rate: param.fs (MHz) --> fs (Hz)
dt  = 1/fs;
dx = pitch;                 % axis X pitch (mm) --> dx (m)
dy = pitch;                 % axis Y pitch  (mm)--> dy (m)  //dx = dy

% Material
ERMv = [c(1),c];
th = [0,layer];                 % layertk (mm) --> th (m)

% Transducer 
% fc      = param.transducer.fc*1e6;      % center frequency (MHz) --> fc (Hz)
fmin  = fmin;               % bandwidth lower limit --> fmin (Hz)
fmax  = fmax;               % bandwidth upper limit --> fmax (Hz)
% Record Time Delay
disp = disp;                % param.delayt (us) --> delayt (s)

%% Reconstruction Parameter
dz = c*dt;  
cmax = max(c);
maxdz = max(dz);            % The reconstructed dz precison is determined by the maxdz

% Reconstructed Image z-axis
z_end = cumsum(th./c/dt); z_end =  round((L - round(z_end(end)))*dt*c(end)/maxdz);
z_indexs = cumsum(round(th/maxdz))+1; z_indexs = [z_indexs,z_indexs(end)+z_end];

maxlayer = length(th);


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
E_P0wk = P0wk(Nw1:Nw2,:,:);         % Truncate the band
w = w(Nw1:Nw2);             

[gkx,gw,gky] = meshgrid(kx,w,ky);   % build kx,ky,w matrice


for layers = 1:maxlayer
%% Wave field extrapolation                                                                                         
gkz2 = gw.^2/ERMv(layers)^2 - gkx.^2-gky.^2; mask = gkz2 >= 0;
% All elements of kz for which (w.^2/ERMv(2)^2 - kx.^2) <0 should therefore be set to zero
gkz = sign(gw).*sqrt(gkz2.*mask);                            
E_P0wk = E_P0wk.*(mask~=0);                         
E_P0wk = E_P0wk.*exp(1i.*gkz*(th(layers)));
                 

%% NUFFT migration
f = gw;
fkz = ERMv(layers+1)*sign(f).*sqrt(gkx.^2+gky.^2 + f.^2/ERMv(layers+1)^2); 

if(density ==1)
fftRF = interp3(gkx,f,gky,E_P0wk,gkx,fkz,gky,'linear',0);
else
fftRF = InterpNUFFT(f,kx,ky,E_P0wk,fkz,gkx,gky,density);
end
% Jacobian (optional)
%kz = (-newN/2:newN/2-1)'/ERMv/fs/newN;
gkz = f/ERMv(layers+1);     
fftRF = fftRF.*gkz./(fkz+eps);
%MASK = (gkx.^2 + gky.^2 + (gkz./(fkz+eps)).^2 )< ((f/ERMv(2)).^2);
%fftRF = bsxfun(@times,fftRF,gkz)./(fkz+eps); 

%% IFFT
newL = round(newL/cmax*c(layers));
fftdata = zeros(newL,newM,newN);
fftdata(Nw1:Nw2,:,:) = (fftRF);
fftRF = fftdata;

% IFFT & Migrated RF
Lz = 1:z_indexs((layers+1))-z_indexs(layers)+1;
Temp = ifftn(ifftshift(fftRF),'nonsymmetric');
migRF(z_indexs(layers):z_indexs(layers+1),:,:) = Temp(Lz,1:M,1:N);
end
end