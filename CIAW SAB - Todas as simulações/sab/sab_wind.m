% computes random horizontal wind speed vector (vwii)
% weather: multiples of 10 knots uniform horizontal wind speed

function vwii = sab_wind(weather)
    dir = rand * 2*pi;
    wind = rand * weather * 10*(1852/3600);
    vwii = wind * [cos(dir), sin(dir), 0]';
end