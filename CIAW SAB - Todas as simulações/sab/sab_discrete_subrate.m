function [zi, zm] = sab_discrete_subrate( z, dt, newdt )
    zm = ceil(newdt/dt);    
    zi = mod(z, zm);
end