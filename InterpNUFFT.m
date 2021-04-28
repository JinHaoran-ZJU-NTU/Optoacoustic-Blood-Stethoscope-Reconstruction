function yi = InterpNUFFT(gf,kx,ky,y,gfkz,gkx,gky,density)
% -- Sinc interpolation --
% Harlan, 1982, "Avoiding interpolation artifacts in Stolt migration"
%                cf. equation (4) page 109
% --
% A p-point sinc interpolator is used
% --

%-- SINC approximation around 0:
% sinc0(x) ~ sinc(x) for x in [-.5,.5]
% sinc0(1) = 1 and sinc0(0.5) = sinc(0.5)
f = gf(:,1,1);f = f(:);             % f axis
df = f(2)-f(1);                     % df;
newdf = df/density;                 % newdf;
newf = f(1):newdf:f(end)+df-newdf;  % newf axis
[L,M,N] = size(y);

% sinc interpolation (sinc convolution)
[Orgf,Newf] = meshgrid(f,newf);
detalgrid = (Newf - Orgf)/df;
Interp = exp(-1i*pi*detalgrid).*sinc(detalgrid);
SincFT = Interp*y(1:L,:);
SincFT = reshape(SincFT,L*density,M,N);

% interpolation(nearest or linear)
[NMmesh,Nfmesh,NNmesh] = meshgrid(kx,newf,ky);
yi = interp3(NMmesh,Nfmesh,NNmesh,SincFT,gkx,gfkz,gky,'near',0);
end