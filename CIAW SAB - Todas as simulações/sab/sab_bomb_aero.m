% Cálculo dos coeficientes aerondâmicos da bomba na condição atual
%   dadas as "derivadas de estabilidade aerodinâmicas"

function [aero] = sab_bomb_aero(speed, Mach, dynamic_pressure, aoa, aos, waib)
    
    % posição, comprimento e superfície de referência aerodinâmica
    aero.Xref = 1.3;
    aero.Lref = 0.36;
    aero.Sref = 0.1;
    aero.ref = [aero.Xref, 0, 0]';

    % deflexão de cada empena
    delta = 3.0 * pi/180;

    % adimensionalização do 
    damping = aero.Lref / (2 * speed);
    p = waib(1);
    q = waib(2);
    r = waib(3);      

    if abs(aoa) > 30*(pi/180) || abs(aos) > 30*(pi/180)
        warning("Ângulo de ataque fora do modelo aerodinâmico.");
    end

    % força no eixo longitudinal (sentido forward)
    CxM = interp1([    0     0.6     0.8     0.9    1.0    1.1    1.2    1.5    2.0   ], ...
                  [-0.12   -0.14   -0.17   -0.29  -0.49  -0.55  -0.52  -0.37  -0.33   ], ...
                  Mach, 'pchip', 'extrap');
    Cx = CxM;

    % força no eixo lateral (sentido right)
    Cy0 = 0;
    Cyb = -5.6;
    Cyr = +42.4;
    Cy = Cy0   +   aos*Cyb   +   r*Cyr*damping;    

    % força no eixo vertical (sentido down)
    Cz0 = 0;
    Cza = -5.6;
    Czq = -42.4;
    Cz = Cz0   +   aoa*Cza   +   q*Czq*damping;

    % torque no eixo longitudinal (rolling moment)
    Cl0 = 0;
    Cla = 0;
    Clb = 0;
    Cld = 0.2;
    Clp = -1.4;
    Cl = Cl0  +  aoa*Cla  +  aos*Clb  +  p*Clp*damping  +  4*Cld*delta;

    % torque no eixo lateral (pitching moment)
    Cm0 = 0;
    Cma = -6.8;
    Cmq = -177;
    Cm = Cm0  +  aoa*Cma  +  q*Cmq*damping;

    % torque no eixo vertical (yawing moment)
    Cn0 = 0;
    Cnb = +6.8;
    Cnr = -177;
    Cn = Cn0  +  aos*Cnb  +  r*Cnr*damping;                                    

    % calculando forças e momentos (torques)
    aero.force  = dynamic_pressure * aero.Sref * [Cx Cy Cz]';
    aero.torque = dynamic_pressure * aero.Sref * aero.Lref * [Cl Cm Cn]'; 
    
    aero.Fx = aero.force(1);
    aero.Fy = aero.force(2);
    aero.Fz = aero.force(3);
    
    aero.Mx = aero.torque(1);
    aero.My = aero.torque(2);
    aero.Mz = aero.torque(3);
    
end