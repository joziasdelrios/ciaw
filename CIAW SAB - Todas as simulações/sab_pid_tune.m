% Laboratório de ajuste (tuning) de ganhos de controlador PID.
%   PID = Proporcional-Integral-Derivativo
%
% Estudo dos movimentos harmônicos (oscilações) amortecidas de um sistema 
%    Massa-Mola-Amortecedor com controlador PID.
%
% Implementa os seguintes recursos:
%  - Saturação do integrador (anti-windup)
%  - Saturação do controlador (limite de comando)
%  - Atraso do controlador (constante de tempo de primeira ordem)
%
% CIAW/MB - Centro de Instrução Almirante Wandenkolk
% ITA/FAB - Instituto Tecnológico de Aeronáutica (DCTA)
% IAE/FAB - Instituto de Aeronáutica e Espaço (DCTA)
% ASD/IAE - Divisão de Sistemas de Defesa
%
% Curso: SAB - Simulação e Controle de Artefatos Bélicos
%
% Instrutor: Jozias DEL RIOS - Cap Eng (ASD/IAE)
%            <delriosjdrvgs@fab.mil.br>
%            <joziasdelrios@gmail.com>
%            12-98177-9921
%
% Rev. 14-dez-2017
%
% Sintaxe:  sab_pid_tune(x_ref, x0,v0, m,b,k, kp,ki,kd, eimax, sat, tc)
%
%    x_ref = valor desejado de controle [m]
%    x0    = posição inicial [m]
%    v0    = velocidade inicial [m/s]
%    m     = massa do modelo [kg]
%    b     = coeficiente de amortecimento [kg/s]
%    k     = constante de mola [N/m]
%    kp    = ganho do controlador proporcional ao erro [N/m]
%    ki    = ganho do controlador integral ao erro [N/m/s]
%    kd    = ganho do controlador derivada do sinal [N.s/m]
%    eimax = valor máximo acumulado no controle integral [m.s]
%    sat   = saturação do comando de controle [N]:
%               F = Kp*(x_ref-x) + Ki*integral(x_ref-x,t) + Kd*derivada(x)
%    tc    = constante de tempo de atraso do sinal de controle [s]
%               use zero para desabilitar.

%{
  % Uso de exemplo: 
  %   Entrada degrau: x0 = v0 = 0;  x_ref = 1;
  %   Controle proporcional: kp = 100;
  %   Modelo com parâmetros: m = 1; b = 5; k = 10;
  %   Sem saturações e atrasos: sat = 100; tc = 0

  sab_pid_tune(1, 0, 0, 1,5,10, 100,0,0, 0, 0, 100, 0);

%}

function sab_pid_tune(x_ref, x0,v0, m,b,k, kp,ki,kd, eimax, sat, tc)
    addpath('sab');
    
    tmax = 5.0;
    
    %% Configuração da simulação
    % passo inicial da simulação e tempo limite
    sim.dt = 0.005;
    sim.tmax = tmax;

    % valores das condições iniciais dos estados a serem integrados
    sim.s.x = x0;
    sim.s.v = v0;
    sim.s.ei = 0;
    sim.s.F = 0;
    
    % funções de cálculo da derivada, de parada e de análise do erro
    sim.dsdt = @dsdt;
    sim.fstep = @fstep;
    
    % valores iniciais das variáveis de saída auxiliar (aux)
    sim.aux.a = 0;
    sim.aux.energy = 0;
    
    % parâmetros dos modelos do sistema massa-mola-amortecedor
    sim.param.m = m;     % massa do corpo de 2 kg
    sim.param.b = b;     % amortecedor: 1 kg/s
    sim.param.k = k;     % constante de mola 5 N/m
    
    sim.param.x_ref = x_ref;   % valor desejado 1 m
    sim.param.Kp = kp;       % ganho do controle proporcional
    sim.param.Ki = ki;       % ganho do controle integral
    sim.param.Kd = kd;       % ganho do controle integral
    sim.param.eimax = eimax; % saturação do integrador (anti-windup)
    sim.param.sat = sat;     % saturação da força de controle
    sim.param.delay = 0.0;   % atraso morto (paralisado)
    sim.param.tc = tc;       % constante de tempo de atuação
            
    %% cálculo dos parâmetros do sistema massa-mola-amortecedor
    fprintf('============================================================\n');
    fprintf('Mass-Spring-Damper Model\n');
    fprintf('  Inputs:\n');
    fprintf('     m     = %f [kg] (mass)\n', m);
    fprintf('     b     = %f [kg/s] (damper) \n', b);
    fprintf('     k     = %f [N/m] (spring)\n', k);
    fprintf('     x0    = %f [m] (initial position)\n', x0);
    fprintf('     x_ref = %f [m] (desired position)\n', x_ref);
    fprintf('     v0    = %f [m/s] (initial velocity)\n', v0);    
    fprintf('     Kp    = %f [N/m] (proportional gain)\n', kp);
    fprintf('     Ki    = %f [N/m/s] (integral gain), saturation = %f [m.s]\n', ki, eimax);
    fprintf('     Kd    = %f [N.s/m] (derivative gain)\n', kd);
    fprintf('     sat   = %f [N] (control signal saturation)\n', sat);
    fprintf('     tc    = %f [s] (control signal time constant)\n', sat);

    if m < 0 || (b+kd) < 0 || (k+kp) < 0
        fprintf('     Closed-loop Divergence.\n');
        return;
    end
    
    zeta = b / (2 * sqrt(k * m));
    wn = sqrt(k/m);
    wd = wn * sqrt(1 - zeta^2);
    phi = -atan( (zeta * wn * x0 + v0) / (x0 * wd) );
    A = x0 * sec(phi);
    
    fprintf('  Open-Loop Analytical Outputs: x(t) = A*exp(-zeta*wn*t)*cos(wd*t+phi)\n');
    fprintf('     zeta  = %f\n', zeta);
    if zeta > 0
        fprintf('     wn    = %f [rad/s] (natural frequency)  ---> Tn=%4.4f [s]\n', wn, 2*pi/wn );
        fprintf('     wd    = %f [rad/s] (dampened frequency) ---> Td=%4.4f [s]\n', wd, 2*pi/wd );
        fprintf('     A     = %f [m] (initial amplitude)\n', A);
        fprintf('     phi   = %f [deg]\n', 2*pi*phi);    
    end
    
    x_inf_pd = kp / (k + kp) * x_ref;
    zeta = (b+kd) / (2 * sqrt((k+kp) * m));
    wn = sqrt( (k+kp) / m);
    wd = wn * sqrt(1 - zeta^2);
    phi = atan( 1/wd * v0 / (x_inf_pd - x0) - zeta / sqrt(1 - zeta^2));
    A = (x0 - x_inf_pd) / cos(phi);
    
    fprintf('  PD Closed-Loop Analytical Outputs: x(t) = A*exp(-zeta*wn*t)*cos(wd*t+phi)\n');
    fprintf('     zeta  = %f\n', zeta);
    if zeta > 0
        fprintf('     wn    = %f [rad/s] (natural frequency)  ---> Tn=%4.4f [s]\n', wn, 2*pi/wn );
        fprintf('     wd    = %f [rad/s] (dampened frequency) ---> Td=%4.4f [s]\n', wd, 2*pi/wd );
        fprintf('     A     = %f [m] (initial amplitude)\n', A);
        fprintf('     phi   = %f [deg] (initial phase)\n', 2*pi*phi);    
        fprintf('     x_inf = %f [m] (steady state)\n', x_inf_pd); 
    end

    %% simula usando RK4
    sim.method = "RK4";
    [t,s,aux,out,tout] = sab_sim(sim);

    %% cálculo do desempenho
    x_inf = s(end).x;        
    td = sim.param.delay;
    xdir = sign(x_ref - x0);
    dx = abs(x_inf - x0);   
    
    fprintf('  Simulation Outputs:\n');
    fprintf('     x_inf = %f (steady value)\n', x_inf); 
    fprintf('     sse   = %f%% (steady state error)\n', 100 * abs(x_inf - x_ref) / (x_ref - x0)  ); 
    
    % calcular o overshoot e peak-time    
    if xdir > 0
        [xp, tpi] = max( [s.x] );
        Mp = (xp - x_inf)/dx * 100;
    else
        [xp, tpi] = min( [s.x] );
        Mp = (x_inf - xp)/dx * 100;
    end
            
    % calcular o rise-time 10% até 90% do valor final médio
    if xdir > 0
        tr10 = find( [s.x] > x0 + 0.1*dx, 1, 'first');
        tr90 = find( [s.x] > x0 + 0.9*dx, 1, 'first');
    else
        tr10 = find( [s.x] < x0 - 0.1*dx, 1, 'first');
        tr90 = find( [s.x] < x0 - 0.9*dx, 1, 'first');
    end        
    tr = t(tr90) - t(tr10);
    
    % calcular o settling time
    if xdir > 0
        tslo = find( [s.x] < x0 + 0.95*dx, 1, 'last');
        tshi = find( [s.x] > x0 + 1.05*dx, 1, 'last');
    else
        tslo = find( [s.x] > x0 - 0.95*dx, 1, 'last');
        tshi = find( [s.x] < x0 - 1.05*dx, 1, 'last');
    end        
    ts = t( max(tslo, tshi) );
    
    %% print do desempenho
    fprintf('     Mp    = %2.4f%% (Overshoot)\n', Mp );
    fprintf('     tp    = %2.4f [s] (Peak time)\n', t(tpi) );
    fprintf('     tr    = %2.4f [s] (Rise time)\n', tr);
    fprintf('     ts    = %2.4f [s] (Settling time 5%% band)\n', ts);        

    %% Gera os gráficos cinemáticos e energia
    close all;
    set(0,'DefaultFigureWindowStyle','docked');
    set(0,'DefaultLineLineWidth',2);

    figure('Name','Massa-Mola-Amortecedor simulado com RK4');
    plot(t, [s.x]); hold on;
    grid on; 
    title('Posição x - Massa-Mola-Amortecedor com controle Proporcional-Integral-Derivativo'); 
    xlabel('[s]'); ylabel('[m]'); 

    cte = ones(1, numel(t));
    plot(t, x_ref * cte, 'k', 'LineWidth', 2);
    plot(t, x_inf * cte, '--', 'Color', [0.6 0.0 0.6], 'LineWidth', 2);

    plot(t, (x0 + 0.95*(x_inf - x0)) * cte, 'r:', 'LineWidth', 2);
    plot(t, (x0 + 1.05*(x_inf - x0)) * cte, 'r:', 'LineWidth', 2);
    plot(t, [s.x], 'b'); 
    
    miny = min([s.x]);
    maxy = max([s.x]);
    
    plot(t, [aux.Fext], 'g:');
    ylim([miny, maxy]);
    


    %----------------------------------------------------------------------
    %% Função que calcula a derivada temporal dos estados:
    % t      tempo atual
    % s      estrutura de estados atual
    % dsp    estrutura de derivada temporal dos estados do passo anterior
    % auxp   estrutura de saídas auxiliares do passo anterior
    % param  estrutura de parâmetros constantes
    % ds     saída: estrutura de derivada dos estados
    % aux    saída: estrutura de saídas auxiliares
    function [ds, aux] = dsdt( t, dt, s, dsp, auxp, param, out )
        
        % cálculo do somatório das forças externas
        % ganho proporcional ao erro entre referencia e valor atual
        if t < sim.param.delay
            % antes de acabar o delay, não aplique forças externas
            aux.Fext = 0;
            ds.ei = 0;
        else
            % sinal de erro e(t)
            e = (param.x_ref - s.x);
                        
            % anti-windup: limitar o módulo do sinal de erro:
            if abs(s.ei) > param.eimax && sign(s.ei)==sign(e)
                % integrador saturado, não acumule mais
                ds.ei = 0;
            else
                % o sinal u(t) é a integral do sinal de erro e(t)
                % então o sinal de erro e(t) é a derivada do sinal u(t)
                ds.ei = e;
            end
            
            % a derivada do erro é o negativo da derivada da posição atual
            aux.Fext = param.Kp * e  +  param.Ki * s.ei  +  param.Kd * (-dsp.x);
            
            % saturação da força externa
            if aux.Fext > param.sat
                aux.Fext = param.sat;
                
            elseif aux.Fext < -param.sat
                aux.Fext = -param.sat;
            end
        end

        % Somatório das forças no corpo:
        %   Força externa (Fext)
        %   Força da mola: -k*x
        %   Força do amortecedor: -b*v
        F = aux.Fext - param.k * s.x - param.b * s.v;
        
        if param.tc > 0
            % aplica a dinâmica de primeira ordem
            ds.F = (1/param.tc) * (F - s.F);
        
            % cálculo da dinâmica de aceleração, 2a Lei de Newton    
            aux.a = s.F / param.m;
        else
            ds.F = 0;
            aux.a = F / param.m;
        end            

        % calcula a energia do sistema em cada instante de tempo
        %   Energia potencial na mola: k*x^2/2
        %   Energia cinética do corpo: m*v^2/2
        aux.energy = param.k * s.x^2/2  +  param.m * s.v^2/2;

        % derivada da velocidade é a aceleração (calculada)
        ds.v = aux.a;
        % derivada da posição é a velocidade (estado)
        ds.x = s.v;
    end

end
