% Simula��o Bal�stica de Bomba em 3DOF
% Instrutor Jozias del Rios - 11.jun.2016
    
% Lan�amento nas condi��es iniciais:
%   Altitude de 10 kft
%   Velocidade inicial horizontal de 300 knots
%   Atitude de 2 graus cabrado (levantado em eleva��o)

%{

    [sim,t,s,aux] = balistico_3dof(0, -10000*0.3048, 300*0.5144,   0, 2*pi/180, 0);
    sab_plot_ballistic_3dof( sim,t,s,aux );

%}
 
function [sim,t,s,aux] = balistico_3dof(x0, z0, vx0, vz0, theta0, q0)

    tic;
    addpath('sab');
    
    %% prepara��o da estrutura de simula��o (sim)
    % par�metros do corpo da bomba
    sim.param.mass = 490;
    sim.param.Jyy = 155;
    sim.param.Xcm = 1.3;

    % estados iniciais de posi��o (alcance x, altitude z)
    sim.s.x = x0;
    sim.s.z = z0;
       
    % convertendo velocidades iniciais (vx0,vz0) para ref. corpo (u0,w0)
    sim.s.u = +vx0*cos(theta0) - vz0*sin(theta0);
    sim.s.w = +vx0*sin(theta0) + vz0*cos(theta0);
    
    % estados iniciais de orienta��o (theta) e velocidade de arfagem (q)
    sim.s.theta = theta0;
    sim.s.q = q0;
    
    % tempo m�ximo de simula��o: 5min. passo de simula��o: 10ms
    sim.tmax = 300;   
    sim.dt = 0.01;
    
    % fun��es do usu�rio:
    sim.dsdt = @dsdt;       % derivada dos estados
    sim.fstop = @fstop;     % condi��o de parada do simulador
    
    %% executar a simula��o
    sim.method = "RK4";
    [t,s,aux,~,~] = sab_sim(sim);
    
    %% resultados num�ricos da simula��o
    r2d = 180/pi;    
    fprintf('Tempo final: %5.2f segundos.\n', t(end));
    fprintf('Bomb-range: %4.2f metros\n', s(end).x );
    fprintf('�ngulo de impacto: %3.2f graus\n', 90 + r2d*s(end).theta );
    
    toc;
    
    %% fun��o de c�lculo derivada dos estados
    function [ds, aux] = dsdt( t, dt, s, dsp, auxp, param, out )       
        
        % velocidades no referencial inercial
        aux.vx = +s.u * cos(s.theta) + s.w * sin(s.theta);
        aux.vz = -s.u * sin(s.theta) + s.w * cos(s.theta);
               
        % par�metros atmosf�ricos locais (altitude atual, deltaISA=zero)
        aux.air = sab_air_simple( -s.z, 0);
        
        % velocidade aerodin�mica, press�o din�mica, Mach, AoA, ...
        aux.speed = sqrt(s.u^2 + s.w^2);
        aux.aoa = atan2(s.w, s.u);
        aux.dynpress = aux.air.density * aux.speed^2 / 2;
        aux.Mach = aux.speed / aux.air.sound_speed;
        
        %% c�lculo do somat�rio das for�as e torques externos
    
            % for�a peso
            Fweight.x = param.mass * aux.air.gravity * (-sin(s.theta));
            Fweight.z = param.mass * aux.air.gravity * (+cos(s.theta));

            %% esfor�os aerodin�micos
                % c�lculo dos coeficientes aerond�micos da bomba na condi��o atual
                %   dadas as "derivadas de estabilidade aerodin�micas"
                Xref = 1.3;
                Lref = 0.36;
                Sref = 0.1;

                if abs(aux.aoa) > 30*(pi/180)
                    warning("�ngulo de ataque fora do modelo aerodin�mico.");
                end

                damping = Lref / (2 * aux.speed);

                % for�a axial (arrasto) no eixo do corpo (x)
                Cx = interp1([    0     0.6     0.8     0.9    1.0    1.1    1.2    1.5    2.0   ], ...
                             [-0.12   -0.14   -0.17   -0.29  -0.49  -0.55  -0.52  -0.37  -0.33   ], ...
                             aux.Mach, 'pchip', 'extrap');

                % for�a vertical no eixo do corpo (z)
                Cz0 = 0;
                Cza = -5.6;
                Czq = -42;
                Cz = Cz0   +   aux.aoa*Cza   +   s.q * Czq * damping;

                % momento (torque) em pitch
                Cm0 = 0;
                Cma = -6.8;
                Cmq = -177;
                Cm = Cm0   +   aux.aoa*Cma   +   s.q * Cmq * damping;            
            
            % aplicando esfor�os (for�as e momentos) aerodin�micos
            Faer.x = aux.dynpress * Sref * Cx;
            Faer.z = aux.dynpress * Sref * Cz;
            Maer.y = aux.dynpress * Sref * Lref * Cm;
            
            % torque aerodin�mico adicional devido � ponto de ref. aer.
            Maer.y = Maer.y + (Xref - param.Xcm) * Faer.z;
            
        %% somat�rio das for�as e torques no corpo
        Fx = Fweight.x + Faer.x;
        Fz = Fweight.z + Faer.z;
        My = Maer.y;

        %% equa��es de din�mica
        ds.x = +s.u * cos(s.theta)  +  s.w * sin(s.theta);
        ds.z = -s.u * sin(s.theta)  +  s.w * cos(s.theta);
    
        ds.u = Fx/param.mass - s.q * s.w;
        ds.w = Fz/param.mass + s.q * s.u;
    
        ds.theta = s.q;
        ds.q = My/param.Jyy;
    end

    %% crit�rio de parada do simulador
    function stop = fstop(t, s, ds, aux, param)        
        
        % interrompa a simula��o se a bomba colidir com o solo (MSL h = 0)
        if s.z >= 0 
            stop = true;
        else
            stop = false;
        end
    end
end
