function yi = InterpSinc(gf,y,gfkz)
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
[L,M,N] = size(gfkz);
for k = 1:M
    for l = 1:N
        fkz = gfkz(:,k,l); fkz = fkz(:);    % fkz axis
        [fgrid,fkzgrid] = meshgrid(f,fkz);
        detalgrid = (fkzgrid - fgrid)/df;
        Interp = exp(-1i*pi*detalgrid).*sinc(detalgrid);
        yi(:,k,l) = Interp*y(:,k,l);
    end
end
end