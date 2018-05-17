% SAB - Parâmetros atmosféricos e gravitacionais locais 
%
% Modelos simples:
%   Atmosfera padrão ISA até 20km / 60kft.
%   Gravidade constante.
% 
% Entradas: 
%   altitude    = altitude ASL (Above Sea Level)
%   deltaISA    = variação adicional relativa à temperatura ISA (15 graus)
%
% Saída: estrutura air com:
%   air.temperature     = temperatura local
%   air.pressure        = pressão estática local
%   air.density         = densidade do ar local
%   air.sound_speed     = velocidade do som no local
%   air.gravity         = aceleração da gravidade constante de 1901   

function [air] = sab_air_simple(altitude, deltaISA) 
    gamma = 1.4;
    R = 8.314472 / 0.029;    
    Tmsl = 288.15;
    p0 = 101325;
    T0 = Tmsl + deltaISA;
    
    air.gravity = 9.80665;

    if altitude <= 11000
        air.temperature = T0 - 0.0065*altitude;
        air.pressure = p0 * (air.temperature / T0) ^ (air.gravity / (0.0065*R));        
    
    elseif altitude <= 20000
        air.temperature = T0 - 71.55;
        p11000 = p0 * ( (T/T0) ^ (air.gravity / (0.0065*R)) );
        air.pressure = p11000 * exp( -air.gravity * (altitude-11000) / (R*air.temperature) );
    else
        error("Modelo atmosférico não definido acima de 20 km");
    end 

    air.density = air.pressure / (R * air.temperature);
    air.sound_speed = sqrt( gamma * R * air.temperature );
end