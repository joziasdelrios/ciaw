% Simulação de Míssil SRAAM (Short-Range Air-to-Air Missile) em 3DOF 
%   (Plano Horizontal) para estudo de Leis de Guiamento
%
%   Em 3DOF horizontal, os corpos tem:
%     Posição (x,y)
%     Velocidades (u,v) nos eixos do corpo: u=frente, v=direita. 
%     Orientação (ângulo psi), sendo que: 
%        psi=0     é na direção de x positivo
%        psi=pi/2  é na direção de y positivo
%
% Rev 17.maio.2018
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

% Entradas: estrutura sraam contendo:
%   mx,ym:      posição inicial do míssil
%   mu:         velocidade inicial do míssil (eixos do míssil)
%   tx,ty:      posição inicial do alvo
%   tpsi:       orientação inicial do alvo
%   tu:         velocidade constante do alvo
%   ta:         aceleração lateral do alvo (constante)
%
%   law:        lei de guiamento:
%                    "AP"    perseguição por atitude
%                    "VP"    perseguição por velocidade
%                    "TPN"   navegação proporcional verdadeira
%                    "PPN"   navegação proporcional pura
%   Kpn:    ganho da navegação proporcional (3 a 5 sugerido)
%   vc:     velocidade de aproximação míssil-alvo constante
%               se indefinido, será utilizado o valor real.

%{
    
    sraam.mx = 0; sraam.my = 0; 
    sraam.mu = 200;
    sraam.tx = 3000;
    sraam.ty = -1000;
    sraam.tpsi = 90*pi/180; 
    sraam.tu = 300;
    sraam.ta = 10;
    sraam.law = "AP";
    sraam.Kpn = 4;
    sraam_3dof( sraam );
    [sim,t,s,aux] = sraam_3dof( sraam ); 
    sab_plot_sraam_3dof( sim,t,s,aux );     

%}

function [sim,t,s,aux] = sraam_3dof( sraam )
    addpath('sab');
    tic;
    
    %% Passos e Tempos da simulação
    sim.dt   = 0.005;           % passo básico de 5 ms
    sim.tmax = 20;              % tempo máximo de simulação

    %% Estados iniciais do simulador
    % estados iniciais do míssil
    sim.s.mx = sraam.mx;       % posição do míssil
    sim.s.my = sraam.my;
    sim.s.mu = sraam.mu;       % velocidade do míssil
    sim.s.mv = 0;    
    sim.s.mpsi = 0;            % orientação do míssil
    sim.s.mr = 0;              % velocidade angular do míssil
    % estados do alvo
    sim.s.tx = sraam.tx;
    sim.s.ty = sraam.ty;
    sim.s.tvx = +sraam.tu * cos(sraam.tpsi); 
    sim.s.tvy = +sraam.tu * sin(sraam.tpsi); 
    % estados do atuador
    sim.s.delta = 0;
    
    %% Parâmetros dos modelos
    sim.aux.fuzed = 0;
    sim.aux.miss_distance = 1e10;
    sim.param.actuator_tc = 0.050;
    
    %% Ganhos dos controladores
    % autopiloto de atitude:
    sim.param.pilot_yaw_KP = 50;
    sim.param.pilot_yaw_KD = 1;
    % autopiloto de aceleração:   
    sim.param.pilot_latax_KP = 2 * 1e-5;
    sim.param.pilot_latax_KD = 2;
    
    %% Executa o simulador
    sim.dsdt = @dsdt;
    sim.fstop = @fstop;
    sim.method = "RK4";
    [t,s,aux,out,tout] = sab_sim(sim);
    
    %% resultado final da simulação
    [dm,~] = min([aux.miss_distance]);
    fprintf('\nSRAAM 3DOF MISSILE SIMULATION...');
    fprintf('\nt = %5.2fs', t(end) );
    fprintf('\nmxy=[%5.2f, %5.2f]m  mvuv=[%4.2f, %4.2f]m/s  mpsi=[%3.2f]deg', ...
        s(end).mx, s(end).my, aux(end).mvx, aux(end).mvy, s(end).mpsi*180/pi );
    fprintf('\ndmiss=%4.2fm\n', dm);
    fprintf('\n=====================================================\n');
    toc;    
    
    
    
    %% função de cálculo derivada dos estados ##################################################################
    function [ds, aux] = dsdt( t, dt, s, dsp, auxp, param, out )
        
        %% converte velocidade do ref corpo para inercial
        aux.mvx = +s.mu * cos(s.mpsi) - s.mv * sin(s.mpsi);
        aux.mvy = +s.mu * sin(s.mpsi) + s.mv * cos(s.mpsi);
        
        %% rumos do míssil e do alvo        
        aux.mgamma  = atan2( aux.mvy, aux.mvx );     % no ref inercial
        aux.tgamma  = atan2( s.tvy, s.tvx );         % no ref inercial        
        
        %% interpolação de massa, inércia e propulsão
        % míssil tem 90kg, 150mm de diâmetro, 2900mm de comprimento,
        % deve normalmente chegar em Mach 3.5 
        p = interp1( [    0    1.9    2.1    6.0   6.2    30]', ... % tempo
                     [ 30e3   30e3   11e3   10e3     0     0;   ... % propulsão
                         90     80     77     60    60    60;   ... % massa
                       1.50   1.30   1.28   1.02  1.00  1.00;   ... % cm
                         70     57     55     40    40    40]', ... % Jzz
                        t, 'pchip', 'extrap' );
        aux.thrust = p(1);
        aux.mass = p(2);
        aux.cm = -p(3);
        aux.Jzz = p(4);
        
        %% cálculos auxiliares
        % geometria míssil-alvo, miss-distance e ângulo LoS (Line of Sight)
        rtmx = s.tx - s.mx;
        rtmy = s.ty - s.my;        
        vtmx = s.tvx - aux.mvx;
        vtmy = s.tvy - aux.mvy;        
        aux.miss_distance = sqrt( rtmx^2 + rtmy^2 );
        aux.los = atan2(rtmy, rtmx);
        
        %% lógica da espoleta
        aux.fuzed = 1;
        if (auxp.fuzed == 0) && (auxp.miss_distance >= aux.miss_distance)
            aux.fuzed = 0;
        end
                
        %% obtém os dados atmosféricos locais, altitude fixa 10 kft
        air = sab_air_simple(10000*0.3048, 0);
        %   air.temperature     = temperatura local
        %   air.pressure        = pressão estática local
        %   air.density         = densidade do ar local
        %   air.sound_speed     = velocidade do som no local
        %   air.gravity         = aceleração da gravidade constante de 1901 (9.80665 m/s^2)
        
        %% velocidade aerodinâmica, pressão dinâmica, Mach, AoA, ...
        aux.speed = sqrt(s.mu^2 + s.mv^2);
        aux.aos = atan2(s.mv, s.mu);
        aux.dynpress = air.density * aux.speed^2 / 2;
        aux.Mach = aux.speed / air.sound_speed;
        
        %% NAVEGAÇÃO: sensores para medir e processadores para calcular:
        %aux.nav_axb = dsp.mu - s.mr * s.mv;    % acelerômetro no eixo X (forward)
        %aux.nav_ayb = dsp.mv + s.mr * s.mu;    % acelerômetro no eixo Y (right)        
        aux.nav_axb = dsp.mu;    % acelerômetro no eixo X (forward)
        aux.nav_ayb = dsp.mv;    % acelerômetro no eixo Y (right)        
        
        %% AUTODIRETOR: rastreando o alvo, faz as leituras: 
        % Leitura de azimute do autodiretor: direção do alvo no ref míssil
        aux.seeker_azimuth = aux.los - s.mpsi;
        % VAILV (LoS Rate) variação angular inercial do LoS obtida analiticamente
        aux.seeker_losrate = (vtmy*rtmx - vtmx*rtmy) / aux.miss_distance;
        % alcance atual até o alvo (seeker Radar)
        aux.seeker_range = aux.miss_distance;
        % velocidade de aproximação (seeker Radar)
        if ~isfield(sraam,'vc'), aux.seeker_rangerate = sqrt( vtmx^2 + vtmy^2 );
        else,                    aux.seeker_rangerate = sraam.vc;
        end
        
        %% GUIAMENTO: conforme a lei de guiamento escolhida
        if sraam.law == "AP"
            %% Guiamento por Perseguição de Atitude:
            % objetivo: zerar a leitura de azimute do autodiretor
            aux.guidance_yaw = aux.seeker_azimuth;
            
        elseif sraam.law == "VP"
            %% Guiamento por Perserguição de Velocidade
            % objetivo: alinhar a velocidade com o azimute do autodiretor
            aux.guidance_yaw = aux.seeker_azimuth - aux.aos;
        
        else
            accel = sraam.Kpn * aux.seeker_rangerate * aux.seeker_losrate;
            if sraam.law == "PPN"
                %% Guiamento por Navegação Proporcional Pura
                % objetivo: manobra LATAX na direção perpendicular à velocidade
                aux.guidance_latax = accel / cos(aux.aos);            

            elseif sraam.law == "TPN"
                %% Guiamento por Navegação Proporcional Verdadeira
                % objetivo: manobra LATAX na direção perpendicular ao LOS                
                aux.guidance_latax = accel / cos(aux.seeker_azimuth);
            end
        end
        
        %% CONTROLE: Autopiloto de atitude (yaw) OU de aceleração (latax)
        if sraam.law == "AP" || sraam.law == "VP"
            % Autopiloto de Atitude
            aux.delta = + param.pilot_yaw_KP * aux.guidance_yaw ...
                        - param.pilot_yaw_KD * s.mr;
        else %if sraam.law == "TPN" || sraam.law == "PPN"
            aux.delta = + param.pilot_latax_KP * (aux.guidance_latax - aux.nav_ayb) ...
                        - param.pilot_latax_KD * s.mr;
        end
        
        % Saturação dos comandos, geralmente em SKID-TO-TURN é de 10 graus
        aux.delta = 4*clamp( aux.delta, 10 );
        if t < 0.3
            aux.delta = 0;
        end
        
        %% cálculo do somatório das forças e torques externos
            % força peso nula (plano horizontal apenas)
    
            % esforços aerodinâmicos: estimação das forças e momentos
            aero.ref  = -1.360;
            aero.Lref = 0.152;
            aero.Sref = pi*(aero.Lref/2)^2;
            damping = aero.Lref / (2*aux.speed);
            
            % força no eixo longitudinal (sentido forward)
            CxM = interp1([    0.5      0.8      1.0      1.2     1.5    2.0     3.0      4.0   ], ...
                          [  -0.49    -0.47    -1.25    -1.35   -1.10   -1.1   -0.90    -0.79   ], ...
                          aux.Mach, 'pchip', 'extrap');
            Cxd = -0.087;
            Cx = CxM + Cxd*abs(s.delta);
            
            % força no eixo lateral (sentido right)
            Cyb = -36.6/180*pi*50;
            Cyr = 517/180*pi*50;
            Cyd = -0.518/4;
            Cy = aux.aos*Cyb  +  s.mr*Cyr*damping  +  Cyd*s.delta;
            
            % torque no eixo vertical (yawing moment)
            Cnb = 70.8/180*pi;
            Cnr = -5812/180*pi*800;
            Cnd = 3.84/10;
            Cn = aux.aos*Cnb  +  s.mr*Cnr*damping  +  Cnd*s.delta;
            
            % somando forças e momentos (torques) aerodinâmicos
            aero.Fx = aux.dynpress * aero.Sref * Cx;
            aero.Fy = aux.dynpress * aero.Sref * Cy;
            aero.Mz = aux.dynpress * aero.Sref * aero.Lref * Cn;
            
        %% forças e momentos totais (peso + aerodinâmico + propulsão)
        Fx = aero.Fx  +  aux.thrust;
        Fy = aero.Fy;
        Mz = aero.Mz + (aero.ref - aux.cm)*Fy;
        
        %% equações de dinâmica
        % do míssil
        ds.mx = +s.mu * cos(s.mpsi)  -  s.mv * sin(s.mpsi);
        ds.my = +s.mu * sin(s.mpsi)  +  s.mv * cos(s.mpsi);    
        ds.mu = Fx/aux.mass + s.mr * s.mv;
        ds.mv = Fy/aux.mass - s.mr * s.mu;    
        ds.mpsi = s.mr;
        ds.mr = Mz/aux.Jzz;
        
        % do atuador
        ds.delta = (1/param.actuator_tc) * (aux.delta - s.delta);
            
        % do alvo
        ds.tx = s.tvx;
        ds.ty = s.tvy;
        ds.tvx = sraam.ta * ( -sin(aux.tgamma) );
        ds.tvy = sraam.ta * ( +cos(aux.tgamma) );
    end    

    %% critério de parada do simulador ##################################################################
    function stop = fstop(t, s, ds, aux, param)        
        stop = false;
        
        % validade do modelo aerodinâmico para os ângulos de ataque
        if abs(aux.aos) > 45*pi/180
           stop = true;
           warning('Interrompendo a simulação: Modelo aerodinâmico extrapolado.');
        end
        
        % se a distância de passagem atender miss_distance
        if aux.fuzed == 1
            stop = true;
        end
    end

end
