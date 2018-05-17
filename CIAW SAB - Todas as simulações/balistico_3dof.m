% Simulação Balística de Bomba em 3DOF
% Instrutor Jozias del Rios - 11.jun.2016
    
% Lançamento nas condições iniciais:
%   Altitude de 10 kft
%   Velocidade inicial horizontal de 300 knots
%   Atitude de 2 graus cabrado (levantado em elevação)

%{

    [sim,t,s,aux] = balistico_3dof(0, -10000*0.3048, 300*0.5144,   0, 2*pi/180, 0);
    sab_plot_ballistic_3dof( sim,t,s,aux );

%}
 
function [sim,t,s,aux] = balistico_3dof(x0, z0, vx0, vz0, theta0, q0)

    tic;
    addpath('sab');
    
    %% preparação da estrutura de simulação (sim)
    % parâmetros do corpo da bomba
    sim.param.mass = 490;
    sim.param.Jyy = 155;
    sim.param.Xcm = 1.3;

    % estados iniciais de posição (alcance x, altitude z)
    sim.s.x = x0;
    sim.s.z = z0;
       
    % convertendo velocidades iniciais (vx0,vz0) para ref. corpo (u0,w0)
    sim.s.u = +vx0*cos(theta0) - vz0*sin(theta0);
    sim.s.w = +vx0*sin(theta0) + vz0*cos(theta0);
    
    % estados iniciais de orientação (theta) e velocidade de arfagem (q)
    sim.s.theta = theta0;
    sim.s.q = q0;
    
    % tempo máximo de simulação: 5min. passo de simulação: 10ms
    sim.tmax = 300;   
    sim.dt = 0.01;
    
    % funções do usuário:
    sim.dsdt = @dsdt;       % derivada dos estados
    sim.fstop = @fstop;     % condição de parada do simulador
    
    %% executar a simulação
    sim.method = "RK4";
    [t,s,aux,~,~] = sab_sim(sim);
    
    %% resultados numéricos da simulação
    r2d = 180/pi;    
    fprintf('Tempo final: %5.2f segundos.\n', t(end));
    fprintf('Bomb-range: %4.2f metros\n', s(end).x );
    fprintf('Ângulo de impacto: %3.2f graus\n', 90 + r2d*s(end).theta );
    
    toc;
    
    %% função de cálculo derivada dos estados
    function [ds, aux] = dsdt( t, dt, s, dsp, auxp, param, out )       
        
        % velocidades no referencial inercial
        aux.vx = +s.u * cos(s.theta) + s.w * sin(s.theta);
        aux.vz = -s.u * sin(s.theta) + s.w * cos(s.theta);
               
        % parâmetros atmosféricos locais (altitude atual, deltaISA=zero)
        aux.air = sab_air_simple( -s.z, 0);
        
        % velocidade aerodinâmica, pressão dinâmica, Mach, AoA, ...
        aux.speed = sqrt(s.u^2 + s.w^2);
        aux.aoa = atan2(s.w, s.u);
        aux.dynpress = aux.air.density * aux.speed^2 / 2;
        aux.Mach = aux.speed / aux.air.sound_speed;
        
        %% cálculo do somatório das forças e torques externos
    
            % força peso
            Fweight.x = param.mass * aux.air.gravity * (-sin(s.theta));
            Fweight.z = param.mass * aux.air.gravity * (+cos(s.theta));

            %% esforços aerodinâmicos
                % cálculo dos coeficientes aerondâmicos da bomba na condição atual
                %   dadas as "derivadas de estabilidade aerodinâmicas"
                Xref = 1.3;
                Lref = 0.36;
                Sref = 0.1;

                if abs(aux.aoa) > 30*(pi/180)
                    warning("Ângulo de ataque fora do modelo aerodinâmico.");
                end

                damping = Lref / (2 * aux.speed);

                % força axial (arrasto) no eixo do corpo (x)
                Cx = interp1([    0     0.6     0.8     0.9    1.0    1.1    1.2    1.5    2.0   ], ...
                             [-0.12   -0.14   -0.17   -0.29  -0.49  -0.55  -0.52  -0.37  -0.33   ], ...
                             aux.Mach, 'pchip', 'extrap');

                % força vertical no eixo do corpo (z)
                Cz0 = 0;
                Cza = -5.6;
                Czq = -42;
                Cz = Cz0   +   aux.aoa*Cza   +   s.q * Czq * damping;

                % momento (torque) em pitch
                Cm0 = 0;
                Cma = -6.8;
                Cmq = -177;
                Cm = Cm0   +   aux.aoa*Cma   +   s.q * Cmq * damping;            
            
            % aplicando esforços (forças e momentos) aerodinâmicos
            Faer.x = aux.dynpress * Sref * Cx;
            Faer.z = aux.dynpress * Sref * Cz;
            Maer.y = aux.dynpress * Sref * Lref * Cm;
            
            % torque aerodinâmico adicional devido à ponto de ref. aer.
            Maer.y = Maer.y + (Xref - param.Xcm) * Faer.z;
            
        %% somatório das forças e torques no corpo
        Fx = Fweight.x + Faer.x;
        Fz = Fweight.z + Faer.z;
        My = Maer.y;

        %% equações de dinâmica
        ds.x = +s.u * cos(s.theta)  +  s.w * sin(s.theta);
        ds.z = -s.u * sin(s.theta)  +  s.w * cos(s.theta);
    
        ds.u = Fx/param.mass - s.q * s.w;
        ds.w = Fz/param.mass + s.q * s.u;
    
        ds.theta = s.q;
        ds.q = My/param.Jyy;
    end

    %% critério de parada do simulador
    function stop = fstop(t, s, ds, aux, param)        
        
        % interrompa a simulação se a bomba colidir com o solo (MSL h = 0)
        if s.z >= 0 
            stop = true;
        else
            stop = false;
        end
    end
end
