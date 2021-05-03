%% Parameter Checks

% NUFFT Density
if(density ~= 2^nextpow2(density))
    error('density setup is incorrect');
end


% c and layer check
len_c = length(c);
len_l = length(layer);

if((len_c-1)~=len_l)
    error('c and layer setup is incorrect');
end

