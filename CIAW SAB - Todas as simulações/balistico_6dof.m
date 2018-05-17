% Simulação Balística de Bomba em 6DOF com Matriz DCM e Quaternions
% Jozias del Rios - 09.set.2016

% pbii = posição do corpo em relação ao referencial inercial, em coordenadas inerciais
% vbii = velocidades lineares do corpo em relação ao referencial inercial, em coordenadas inerciais
% euler = ângulos de Euler da orientação em relação ao referencial inercial, sendo: phi (roll), theta (pitch), psi (yaw)
% wbib = velocidades angulares do corpo em relação ao referencial inercial, em coordenadas do corpo
% weather = nível de 0 à 7 de severidade do clima (vento, rajada, turbulência)

%{
 [sim,t,s,aux,out,tout] = balistico_6dof( [0, 0, -10000*0.3048]', [300*0.5144, 0, 0]', [0,2,0]'*pi/180, [0,0,0]', 3);
 sab_plot_ballistic_6dof( sim,t,s,aux,out,tout );
%}


function [sim,t,s,aux,out,tout] = balistico_6dof( pbii0, vbii0, euler0, wbib0, weather )

    addpath('sab');
    tic;
        
    %% parâmetros do corpo da bomba
    sim.param.mass = 490;
    sim.param.inertia = diag([8.2 155 155]);            % matriz tensor de inércia (Jxx Jyy Jzz)
    sim.param.inv_inertia = inv(sim.param.inertia);     % inverso do tensor de inércia
    sim.param.cm = [1.300, 0.0, 0.0]';                  % posiçãodo centro de massa
    sim.param.vwii = sab_wind(weather);                 % vento constante aleatório
    sim.param.weather = weather;                        % intensidade de perturbação climática
    
    %% preparando processos discretos: 
    sim.freq = 100;                       % frequência de atualização dos processos discretos
    sim.out.turb = sab_turbulence_init;   % inicializa as variáveis de turbulência
    
    %% preparando os estados nas condições iniciais
    sim.s.pbii = pbii0;                  % posição inicial, no referencial inercial, coordenadas inerciais
    sim.s.qbi  = sab_euler2quat(euler0); % converter os ângulos de Euler iniciais para o quaternion de orientação do corpo em relação ao referencial inercial
    sim.s.wbib = wbib0;                  % velocidade angular inicial, no referencial inercial, coordenadas do corpo
    Rbi0 = sab_euler2dcm(euler0);        % converter os ângulos de Euler iniciais para a matriz de rotação de orientação do corpo em relação ao referencial inercial
    sim.s.vbib = Rbi0 * vbii0;           % converter a velocidade linear inicial no referencial inercial de coordenadas inerciais (vbii) para corpo (vbib)
    
    %% passo da simulação (10ms) e tempo máximo (5min)
    sim.dt   = 0.01;
    sim.tmax = 300;
    
    %% funções de operação e simulação
    sim.dsdt = @dsdt;
    sim.fstop = @fstop;
    sim.adjust = @fadjust;
    sim.fdiscrete = @fdiscrete;
    sim.method = "RK4";
    [t,s,aux,out,tout] = sab_sim(sim);
       
    %% resultado final da simulação
    r2d = 180/pi;
    fprintf('t=%5.2fs \timp=(%4.2f, %4.2f)  ', t(end), s(end).pbii(1), s(end).pbii(2));        
    fprintf('\tphi=%3.2f theta=%3.2f psi=%3.2f  ', aux(end).euler(1)*r2d, aux(end).euler(2)*r2d, aux(end).euler(3)*r2d  );
    fprintf('\tvwii=(%4.2f,%4.2f)\n', sim.param.vwii(1), sim.param.vwii(2));
    toc;
        
    %% função de cálculo derivada dos estados
    function [ds, aux] = dsdt( t, dt, s, dsp, auxp, param, out )
        %% cálculos iniciais de cinemática
        % converte de quaternion para matriz de rotação (DCM)
        Rbi = sab_quat2dcm( s.qbi );
        
        % converte a velocidade linear de ref. corpo para ref. inercial
        aux.vbii = Rbi' * s.vbib;
        
        % obtém os ângulos de Euler da orientação atual para os gráficos
        aux.euler = sab_quat2euler( s.qbi );

        %% cálculos atmosféricos
        % parâmetros atmosféricos locais (altitude atual, deltaISA=zero)
        h = -s.pbii(3);
        aux.air = sab_air_simple( h, 0);
        
        % velocidade aerodinâmica no sistema de coordenadas do corpo
        vwib = Rbi * (out.turb.vwii + param.vwii);
        aux.vaib = s.vbib - vwib;
        
        % velocidade aerodinâmica, pressão dinâmica, Mach, AoA, AoS, ...
        aux.velocity = norm( s.vbib );          % velocity é a velocidade absoluta do corpo (ref. inercial)
        aux.speed = norm( aux.vaib );           % speed é o TAS (True Air Speed): módulo da velocidade relativa ao ar.
        aux.dynpress = aux.air.density * aux.speed^2 / 2;
        aux.Mach = aux.speed / aux.air.sound_speed;
        aux.aoa = atan2(aux.vaib(3), aux.vaib(1));      % ângulo de ataque (w/u)
        aux.aos = atan2(aux.vaib(2), aux.vaib(1));      % ângulo de derrapagem (v/u)
        
        %% cálculo do somatório das forças e torques externos
        
        % força peso ocorre na direção do eixo +z (vertical inercial)
        %   é necessáro converter o eixo para as coordenadas do corpo
        weight.force = param.mass * aux.air.gravity * Rbi * [0 0 1]';

        % esforços aerodinâmicos: função de estimação das forças e momentos
        [aero] = sab_bomb_aero(aux.speed, aux.Mach, aux.dynpress, aux.aoa, aux.aos, s.wbib);

        %% somatório das forças e momentos totais
        forces = weight.force + aero.force;
        torques = aero.torque + cross(aero.ref - param.cm, aero.force);
        
        %% equações de dinâmica
        ds.pbii = Rbi' * s.vbib;
        ds.vbib = forces/param.mass - cross(s.wbib, s.vbib);
        ds.qbi  = (1/2) * sab_quatmult( s.qbi, [0 s.wbib'] );
        ds.wbib = param.inv_inertia * (torques - cross(s.wbib, param.inertia * s.wbib));
        
    end

    %% critério de parada do simulador
    function stop = fstop(t, s, ds, aux, param)        
        
        % interrompa a simulação se colidir com o solo
        h = -s.pbii(3);
        if  h < 0 
            stop = true;
        else
            stop = false;
        end
    end

    %% ajuste nos estados feito depois de cada passo:
    function [s] = fadjust(t, s, ds, aux, param)
        
        % apenas normalize o quaternion depois de cada passo
        s.qbi = sab_quatnormalize( s.qbi );        
    end

    %% processos discretos (não integrados)
    function [out] = fdiscrete( t, dt, s, ds, aux, param, outp, z )
        
        % segure os valores discretos anteriores que não forem alterados...
        out = outp(end);
        
        % atualize o processo de turbulência de Dryden
        out.turb = sab_turbulence( param.vwii, aux.vbii, t, dt, -s.pbii(3), param.weather, outp, z);
    end

end

