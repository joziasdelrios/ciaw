% Simulação do Míssil MBDA Exocet MM-40 Block II 
%   Em 6DOF com Matriz DCM e Quaternions
% Rev.: 15/05/2018
%
%   Opções de guiamento:
%       Perseguição de atitude (Attitude Pursuit, AP)
%       Perseguição de velocidade (Velocity Pursuit, VP)
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
% Dados de Modelo do Exocet MM-40 Block II:   
%       Length: 5800 mm,
%       Tangent-ogive Nose: Length 460 mm, Blunted radius 5 mm
%       Diameter: 348 mm
%       Elevators: HEX,XLE=5.450,5.635,SSPAN=0.,0.280,CHORD=0.275,0.085
%       Weight: 855 kg (launch) 
%       Propulsion: boost-sustain solid propellant rocket
%           Boost:   dT = 3.00 s, dV = 310 m/s
%           Sustain: dT =  210 s
%       Range: 50-70 km
%       Speed: Mach 0.93 max
%       Altitude: Cruise: 30 m    Attack altitude: 8 m
%       Warhead: 165 kg HE SAP
%       Guidance: Way-points Sea-skimming,  Proportional Navigation
%       Navigation: Inertial navigation (INS), Radio-altimeter: up to 500 m
%       Seeker: Active X-band monopulse radar, 15 km range , 30deg FoV cone

%{
    % Selecione as linhas abaixo e aperte F9.
    % Depois da primeira vez, poderá utilizar F5.

    pbii0 = [0, 0, -5]';
    euler0 = [0, 20, 0]' * pi/180;
    vbib0 = 0.5144 * [5, -5, 0]';
    ptii0 = [4000, -3500, 0]';
    vtii = 0.5144 * [-100, 0, 0]';
    weather = 2;
    noise = 0;
    guidance_law = "AP";
    %guidance_law = "VP";
    [sim,t,s,aux,out,tout] = exocet_6dof( pbii0, vbib0, euler0, ptii0, vtii, weather, noise, guidance_law );

%}

function [sim,t,s,aux,out,tout] = exocet_6dof( pbii0, vbib0, euler0, ptii0, vtii, weather, noise, guidance_law )
    tic;
    rng(1001);
    addpath('sab');
    r2d = 180/pi;
    d2r = pi/180;
     
    %% passos e tempos da simulação
    sim.dt   = 0.005;           % passo básico de 5 ms
    sim.tmax = 20;              % tempo máximo de simulação: somente o necessário
    sim.freq = 100;             % frequência de atualização dos processos discretos    

    %% parâmetros do ambiente
    sim.param.vwii = sab_wind(weather);         % vento constante aleatório
    sim.param.weather = weather;                % intensidade de perturbação climática    
    
    %% parâmetros de subsistemas do míssil
    sim.param.seeker_max_range = 15000;         % distância máxima de detecção do autodiretor
    sim.param.seeker_max_look_angle = 15*d2r;   % meio-cone do FoV de detecção radar
    sim.param.fuze_detection = 10;              % distância de detecção de proximidade da espoleta
    sim.param.actuator_tc = 0.05;               % constante de tempo de resposta dos atuadores
    sim.param.actuator_unlock_delay = 0.8;      % atraso para liberação dos profundores
    sim.param.actuator_soft_limit = 20;         % batente de software do atuador
    sim.param.actuator_hard_limit = 21;         % batente mecânico do atuador
        
    %% parâmetros do alvo: velocidade inercial constante
    sim.param.vtii = vtii;
    
    %% preparando processos discretos: 
    sim.out.turb = sab_turbulence_init;       % inicializa as variáveis de turbulência
    
    %% preparando os estados nas condições iniciais
    % estados do corpo do míssil:
    wbib0 = 15*d2r * randn(3,1);              % velocidades angulares iniciais aleatórias (balanço)
    sim.s.pbii = pbii0;                       % posição inicial, no referencial inercial, coordenadas inerciais
    sim.s.qbi  = sab_euler2quat(euler0);      % converter os ângulos de Euler iniciais para o quaternion de orientação do corpo em relação ao referencial inercial
    sim.s.wbib = wbib0;                       % velocidade angular inicial, no referencial inercial, coordenadas do corpo
    Rbi0 = sab_euler2dcm(euler0);             % converter os ângulos de Euler iniciais para a matriz de rotação de orientação do corpo em relação ao referencial inercial
    sim.s.vbib = vbib0;                       % converter a velocidade linear inicial no referencial inercial de coordenadas inerciais (vbii) para corpo (vbib)
    
    % estados do corpo do alvo: posição
    sim.s.ptii = ptii0;
    
    % estados de deflexão angular dos atuadores dos profundores
    sim.s.actuator_delta = [0 0 0 0]';
    sim.out.actuator_delta_ref = [0 0 0 0]';
    
    % estados iniciais de navegação inercial: transferência de alinhamento
    sim.out.nav_qbi = sim.s.qbi;
    sim.out.nav_vbii = Rbi0' * vbib0;
    sim.out.nav_pbii = pbii0;

    %% erros e desvio-padrão dos modelos (run-to-run bias/noise)
    sim.param.thrust_error   = noise * 1.0   *d2r*randn(2,1);    % desvio do vetor de empuxo: 1.0 graus
    sim.param.imu_accel_bias = noise * 0.1   *9.81*randn(3,1);   % bias dos acelerômetros: 10 mg
    sim.param.imu_gyro_noise = noise * 0.1   *d2r;               % ruído dos girômetros: 0.1 dps
    sim.param.seeker_noise   = noise * 0.3   *d2r;               % ruído nas leituras do autodiretor (azimute/elevação): 0.3 graus
    sim.param.radalt_noise   = noise * 2.0;                      % ruído na leitura do radar-altímetro: 2 m
    
    
    %% Ganhos das malhas de Controle (AutoPilotos)
    % controle de rolamento (roll):
    sim.param.pilot_roll_KP  = 10;            % ganho proporcional 
    sim.param.pilot_roll_KD  = 1;             % ganho derivativo
    % controle de arfagem (pitch):
    sim.param.pilot_pitch_KP  = +0.100;        % ganho proporcional 
    sim.param.pilot_pitch_KI  = +0.030;        % ganho integral
    sim.param.pilot_pitch_KHD = -0.100;        % ganho derivativo de altitude: dz/dt
    sim.param.pilot_pitch_KWD = -1.000;        % ganho derivativo de velocidade angular: q
    sim.param.pilot_pitch_KI_sat = 30;         % saturação anti-windup do integrador
    sim.out.alt_err_integral = 0;              % inicialização do integrador
    % controle de guinada (yaw):
    sim.param.pilot_yaw_KP   = +45.0;          % ganho proporcional
    sim.param.pilot_yaw_KD   =   3.0;           % ganho derivativo
    
    %% funções de operação e simulação
    sim.dsdt = @dsdt;
    sim.fstop = @fstop;
    sim.adjust = @fadjust;
    sim.fdiscrete = @fdiscrete;
    sim.method = "RK4";
    [t,s,aux,out,tout] = sab_sim(sim);
    
    %% resultado final da simulação
    [dm,~] = min([aux.miss_distance]);
    fprintf('\n=====================================================');
    fprintf('\nt=%5.2fs --> pbii=[%5.2f, %5.2f, %5.2f]m  vbii=[%4.2f, %4.2f, %4.2f]m/s ', ...
        t(end), s(end).pbii(1), s(end).pbii(2), s(end).pbii(3), aux(end).vbii(1), aux(end).vbii(2), aux(end).vbii(3) );
    fprintf('\nphi=%3.2f theta=%3.2f psi=%3.2f', aux(end).euler(1)*r2d, aux(end).euler(2)*r2d, aux(end).euler(3)*r2d  );
    fprintf('\nvwii=(%3.2f,%3.2f)  dmiss=%4.2fm\n', sim.param.vwii(1), sim.param.vwii(2), dm);
    toc;
    sab_plot_exocet_6dof( sim,t,s,aux,out,tout );
    
    
    %% função de cálculo derivada dos estados ##################################################################
    function [ds, aux] = dsdt( t, dt, s, dsp, auxp, param, out )

        %% cálculos iniciais de cinemática do míssil        
        Rbi = sab_quat2dcm( s.qbi );                % converte de quaternion para matriz de rotação (DCM)
        aux.vbii = Rbi' * s.vbib;                   % converte a velocidade linear de ref. corpo para ref. inercial
        aux.euler = sab_quat2euler( s.qbi );        % obtém os ângulos de Euler da orientação atual para os gráficos
        aux.velocity = norm( s.vbib );              % velocity é a velocidade absoluta do corpo (ref. inercial)
        aux.gamma_yaw =   atan2(  aux.vbii(2), aux.vbii(1) );  % rumo no plano horizontal
        aux.gamma_pitch = atan2( -aux.vbii(3), aux.vbii(1) );  % rumo no plano vertical

        %% cálculos adicionais da geometria míssil-alvo
        % Distância de passagem: 
        aux.ptbi = s.ptii - s.pbii;
        aux.miss_distance = norm(aux.ptbi);
        
        %% cálculos interpolados de massa, inércia e propulsão (estimativa grosseira)         
        % Thrust = g * Isp * dm/dt,  Isp=230s,
        % Boost:   dV=310s, 10.5g boost => 3.0s boost, F=88.35kN, 
        %          mass_boost = 100kg, dm/dt = -100kg/3s = -33.3kg/s
        % Sustain: F=3514N => mass_sustain = 278kg, dm/dt=-1.324 kg/s
        p = interp1([    0    0.1    1.0    3.0   3.2   213   214   215 ]', ... % current time
                    [  855    853    822    755   754   477   477   477 ;   ... % p(1) = aux.mass
                      3.10   3.10   3.08   2.95  2.94  2.09  2.07  2.07 ;   ... % p(2) = -aux.cm.x
                         0  88000  88500  88300  3514  3514     0     0 ;   ... % p(3) = aux.thrust(1)
                      12.9   12.9   11.9   10.9  10.5   7.3   7.2   7.2 ;   ... % p(4) = Jxx
                      1924   1920   1800   1500  1490  1080  1073  1073 ]', ... % p(5) = Jyy=Jzz
                    t, 'pchip','extrap');                
        
        aux.mass = p(1);                   % massa atual
        aux.cm = [-p(2), 0, 0]';           % posição do Centro de Massa
        aux.thrust_force = p(3) * [1, param.thrust_error(1), param.thrust_error(2)]';
        aux.thrust_ref = [-5.800 0 0]';    % ponto de aplicação da força propulsão
        aux.Jxx = p(4); 
        aux.Jyy = p(5); 
        aux.Jzz = aux.Jyy;
        inertia = diag([aux.Jxx, aux.Jyy, aux.Jzz]);
        inv_inertia = diag([1/aux.Jxx, 1/aux.Jyy, 1/aux.Jzz]);        

        %% cálculos atmosféricos com vento e turbulência
        % parâmetros locais (altitude atual, deltaISA=zero)
        h = -s.pbii(3);
        aux.air = sab_air_simple(h, 0);
        %   aux.air.temperature     = temperatura local
        %   aux.air.pressure        = pressão estática local
        %   aux.air.density         = densidade do ar local
        %   aux.air.sound_speed     = velocidade do som no local
        %   aux.air.gravity         = aceleração da gravidade constante de 1901        
        vwib = Rbi * (out.turb.vwii + param.vwii);          % vwib: velocidade total do ar no sistema de coordenadas do corpo
        
        %% cálculos aerodinâmicos locais: velocidade, pressão dinâmica, Mach, AoA, AoS
        aux.vaib = s.vbib - vwib;                           % vaib: velocidade aerodinâmica no sistema de coordenadas do corpo
        aux.speed = norm( aux.vaib );                       % speed é o TAS (True Air Speed): módulo da velocidade relativa ao ar.
        aux.dynpress = aux.air.density * aux.speed^2 / 2;   % pressão dinâmica (q_inf)
        aux.Mach = aux.speed / aux.air.sound_speed;         % número de Mach
        aux.aoa = atan2(aux.vaib(3), aux.vaib(1));          % ângulo de ataque (plano vertical do corpo) (w/u)
        aux.aos = atan2(aux.vaib(2), aux.vaib(1));          % ângulo de derrapagem (plano horizontal do corpo) (v/u)                
        
        %% cálculo do somatório das forças e torques externos
        
        % força peso no eixo vertical inercial, converter para o ref corpo
        weight.force = aux.mass * aux.air.gravity * Rbi * [0 0 1]';

        %% esforços aerodinâmicos: função de estimação das forças e momentos
        % posição, comprimento e superfície de referência aerodinâmica
            aero.Lref=1.280;
            aero.Sref=0.554;
            aero.ref=[-2.500, 0, 0]';
            
            % monta deflexões de manobra (roll, pitch, yaw, brake)
            % a partir das deflexões atuais dos 4 profundores (1,2,3,4)
            aero.delta_brake = (1/4) * sum(abs(s.actuator_delta));
            aero.delta_pitch = (1/4) * [-1 -1 +1 +1] * s.actuator_delta;
            aero.delta_yaw   = (1/4) * [+1 -1 -1 +1] * s.actuator_delta;
            aero.delta_roll  = (1/4) * [-1 -1 -1 -1] * s.actuator_delta;
            
            % aplica aerodinâmica apenas para velocidade relevante
            if aux.speed < 10
                aero.damping = 0;  aux.aoa = 0; aux.aos = 0;
            else
                % amortecimento comum aos coeficientes de derivadas dinâmicas
                aero.damping = aero.Lref / (2 * aux.speed);
            end

            % força longitudinal (sentido forward) (arrasto) (Cx = -CA) 
            CxMach  = interp1([   0.20     0.40     0.60     0.70     0.80      0.90     0.95     1.00    1.10], ...
                              [-0.0482  -0.0430  -0.0400  -0.0390  -0.0382   -0.0495  -0.0635  -0.0784  -0.109], ...
                              aux.Mach, 'pchip', 'extrap');
                          
            Cxd = -(0.101 - 0.0495)/10;       % CA (BRAKE10) - CA (NOMINAL MACH 0.90)   [degrees]
            Cx = CxMach + Cxd * aero.delta_brake;

            % força lateral (sentido right) (Cy = CY)            
            Cyb = -4.99*2;        % =CNA @ alpha=10deg, Mach=0.90 
            Cyr = 3.42*2;         %  CYR @ alpha= 0deg, Mach=0.90
            Cyd = -0.2269/10;     %  CY  @ alpha= 0deg, Mach=0.90 (YAW10) (opposite) (degrees)
            Cy = aux.aos * Cyb   +   s.wbib(3) * Cyr * aero.damping  +  Cyd * aero.delta_yaw;

            % força no eixo vertical (sentido down) (Cz = -CN)
            Cza = Cyb; 
            Czq = -Cyr;          % signal convention
            Czd = Cyd;           % (degrees)
            Cz = aux.aoa * Cza   +   s.wbib(2) * Czq * aero.damping  +  Czd * aero.delta_pitch;

            % torque no eixo longitudinal (rolling moment)
            % (Cl=CLL, Clp=CLLP, Clb=CLLB)
            Clp = -0.65;         % CLLP @ alpha=10deg, Mach=0.90 (NOMINAL)
            Cld = 0.0385/10;     % CLL  @ alpha= 0deg, Mach=0.90 (ROLL10) (degrees)
            Cl = s.wbib(1) * Clp * aero.damping  +  Cld * aero.delta_roll;

            % torque no eixo lateral (pitching moment) 
            Cma = -2.836*2;      % CMA @ alpha=average, Mach=0.90 (NOMINAL) 
            Cmq = -19.34*4;      % CMQ @ alpha=average, Mach=0.90 (NOMINAL)
            Cmd = 0.5605/10*2;   % (degrees)
            Cm = aux.aoa * Cma  +  s.wbib(2) * Cmq * aero.damping  +  Cmd * aero.delta_pitch;

            % torque no eixo vertical (yawing moment)
            Cnb = -Cma;          % =CLNB
            Cnr = Cmq;
            Cnd = Cmd;           % (degrees)
            Cn = aux.aos * Cnb  +  s.wbib(3) * Cnr * aero.damping  +  Cnd * aero.delta_yaw; 

            % somando forças e momentos (torques) aerodinâmicos
            aero.force  = aux.dynpress * aero.Sref * [Cx Cy Cz]';
            aero.torque = aux.dynpress * aero.Sref * aero.Lref * [Cl Cm Cn]';
            
        %% forças e momentos totais (peso + aerodinâmico + propulsão)
        forces = weight.force + aero.force + aux.thrust_force;
        torques = aero.torque + ...
                  cross(aero.ref - aux.cm, aero.force) + ...        % trazer força aerodinâmica para o Centro de Massa
                  cross(aux.thrust_ref - aux.cm, aux.thrust_force); % torque causado pela propulsão
        
        %% Derivada dos estados: Equações de Dinâmica
        % Dinâmica do corpo do míssil
        ds.pbii = Rbi' * s.vbib;
        ds.vbib = forces/aux.mass - cross(s.wbib, s.vbib);
        ds.qbi  = (1/2) * sab_quatmult( s.qbi, [0 s.wbib'] );
        ds.wbib = inv_inertia * (torques - cross(s.wbib, inertia * s.wbib));
        aux.abib = ds.vbib;
        
        % Cinemática do corpo do alvo (velocidade constante => M.R.U.)
        ds.ptii = param.vtii;
        
        % Dinâmica de primeira ordem dos atuadores dos profundores
        ds.actuator_delta = (1/param.actuator_tc) * (out.actuator_delta_ref - s.actuator_delta);
    end

    %% critério de parada do simulador ##################################################################
    function stop = fstop(t, s, ds, aux, param)        
        stop = false;
        
        % se colidir com o solo
        if  s.pbii(3) > 0
            stop = true; 
        end
        
        % se a distância de passagem atender miss_distance (distância reta-ponto)
        if aux.miss_distance <= param.fuze_detection            
            stop = true;
        end
        
        % validade do modelo aerodinâmico para os ângulos de ataque
        if (aux.speed > 10) && (abs(aux.aoa) > 45*d2r || abs(aux.aos) > 45*d2r)
           stop = true;
           warning('Interrompendo a simulação: Modelo aerodinâmico extrapolado.');
        end

    end

    %% ajuste nos estados feito depois de cada passo ##################################################################
    function [s] = fadjust(t, s, ds, aux, param)
        
        % normalize o quaternion de orientação do míssil depois de cada passo.
        s.qbi = sab_quatnormalize( s.qbi );        
        
        % aplica o batente mecânico nos profundores.
        s.actuator_delta = clamp( s.actuator_delta, param.actuator_hard_limit );
        
    end

    %% processos discretos (não integrados) ##################################################################
    % nesta função consta também todo o software embarcado no míssil
    function [out] = fdiscrete( t, dt, s, ds, aux, param, outp, z )
        lastz = max(1, z-1);
        out = outp(lastz);     % inicializa o campo out
        
        % obtenha a matriz DCM do míssil, para utilizar depois
        Rbi = sab_quat2dcm( s.qbi );
        ptbb = Rbi * aux.ptbi;
        
        %% TURBULÊNCIA:
        % atualize o processo de turbulência atmosférica de Dryden
        out.turb = sab_turbulence( param.vwii, aux.vbii, t, dt, -s.pbii(3), param.weather, outp, z);
                        
        %% NAVEGAÇÃO: Leitura da IMU (sensores inerciais)
        % leituras dos sensores da IMU com inserção de ruídos, bias e gravidade
        out.imu_accel_abib = ds.vbib + param.imu_accel_bias + (aux.air.gravity * Rbi * [0 0 1]');
        out.imu_gyro_wbib  = s.wbib + param.imu_gyro_noise * randn(3,1);
                
        %% NAVEGAÇÃO: Integração do girômetro para obter a Orientação estimada
        % integre o navegador de orientação a partir do girômetro
        dqbi = (1/2) * sab_quatmult( out.nav_qbi, [0 out.imu_gyro_wbib'] );
        out.nav_qbi = sab_quatnormalize( out.nav_qbi + dt*dqbi );
        out.nav_euler = sab_quat2euler( out.nav_qbi );
        nav_Rbi = sab_quat2dcm( out.nav_qbi );
        
        %% NAVEGAÇÃO: Integração dupla do acelerômetro para obter a posição e velocidades estimadas
        % Estime a aceleração da gravidade local        
        out.nav_gbii = 9.80665 * [0 0 1]';
        out.nav_gbib = nav_Rbi * out.nav_gbii;
        
        % Desconte a gravidade e a força de Coriollis do acelerômetro
        out.nav_abib = out.imu_accel_abib ...
                       + cross(out.imu_gyro_wbib, nav_Rbi * out.nav_vbii) ...
                       - out.nav_gbib;
        out.nav_abii = Rbi' * out.nav_abib;

        % Integração linear (método de Euler, impreciso)
        out.nav_vbii = out.nav_vbii + dt*out.nav_abii;      % integra a aceleração
        out.nav_pbii = out.nav_pbii + dt*out.nav_vbii;      % integra a velocidade
        out.nav_vbib = nav_Rbi * out.nav_vbii;

        %% NAVEGAÇÃO: Leitura do Radar-Altímetro
        % sofre influência da inclinação (roll) atual do míssil, 
        %   pois inclina o feixe radar, desviando da vertical local.
        out.nav_radalt = -s.pbii(3) / cos(aux.euler(1)) + param.radalt_noise*randn;
        % TODO: trabalho futuro: Fusão de Dados: Filtro de Kalman RADALT <-> IMU
        
        %% AUTODIRETOR: Obtém elevação/azimute da direção do alvo rastreado        
        out.seeker_range = norm( ptbb );        
        out.seeker_azimuth =    atan2(  ptbb(2), sqrt(ptbb(1)^2 + ptbb(3)^2) ) + param.seeker_noise * randn;
        out.seeker_elevation =  atan2( -ptbb(3), ptbb(1) )                     + param.seeker_noise * randn;        
        % converte ângulos (azimuth,elevation) para (look,roll)
        out.seeker_look_angle =  acos( cos(out.seeker_azimuth) * cos(out.seeker_elevation) );
        out.seeker_roll_angle = atan2( sin(out.seeker_elevation), tan(out.seeker_azimuth) );

        % Obtém versor direção de apontamento do seeker, no referencial do corpo
        out.seeker_ptbb = [ cos(out.seeker_elevation) * cos(out.seeker_azimuth), ...
                            sin(out.seeker_azimuth), ...
                            sin(out.seeker_elevation) * cos(out.seeker_azimuth) ]';

        % Verifica se o alvo pode ser detectado (alcance e FoV)
        if (out.seeker_range > param.seeker_max_range) || (out.seeker_look_angle > param.seeker_max_look_angle)
            out.seeker_tracking = false;        
        else
            out.seeker_tracking = true;
        end

        %% GUIAMENTO: Calcula  proa ou rumo desejáveis
        %  Este guiamento supõe que o rolamento esteja controlado em zero,
        %       portanto, supõe que o eixo pitch comanda altitude apenas
        %                     e que o eixo de yaw comanda rumo ou proa.
        if guidance_law == "AP"
            % Lei de Guiamento por Perseguição de Atitude (Attitude Pursuit Guidance)
            % ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            % Comanda a guinada no sentido de zerar o ângulo de azimute do autodiretor
            % --> objetivo centralizar a mira no alvo --> azimute desejado nulo
            out.guidance_yaw = out.seeker_azimuth;
        
        elseif guidance_law == "VP"
            % Lei de Guiamento por Perseguição de Velocidade (Velocity Pursuit Guidance)
            % ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            % Calcula a direção da velocidade do míssil usando a navegação
            % Comanda a guinada para apontar a velocidade na direção da Linha 
            %   de Visada obtida pelo autodiretor
            out.guidance_gamma = atan2( out.nav_vbib(2), out.nav_vbib(1) );
            out.guidance_yaw = (out.seeker_azimuth - out.guidance_gamma);
        end
        
        %% AUTOPILOTO DE ATITUDE de ROLAMENTO: Controle PD de ROLL-HOLD-ZERO
        % mantém míssil com rolamento phi=0 
        % => eixo down do míssil controlado para ficar para baixo
        phi_ref = 0;
        out.pilot_roll = + param.pilot_roll_KP * (phi_ref - out.nav_euler(1) ) ...
                         - param.pilot_roll_KD * out.imu_gyro_wbib(1);
        
        %% AUTOPILOTO VERTICAL: ALTITUDE-LEVEL-HOLD
        %     mantém altura de referência, usando o radar-altímetro
        %     a referência é interpolada na distância restante
        
        altitude_hold = interp1( [ 0   30   150    500    1000    3000 ], ...
                                 [ 0  2.0     8     15      30      50 ], ...
                                 out.seeker_range, 'previous', 'extrap' );        
        
        % Calcular o erro atual de altitude
        out.alt_err = altitude_hold - out.nav_radalt;
        
        % aplica a lei de controle com ganhos
        out.pilot_pitch = param.pilot_pitch_KP  * out.alt_err + ...
                          param.pilot_pitch_KI  * out.alt_err_integral + ...
                          param.pilot_pitch_KHD * (-out.nav_vbii(3)) + ...
                          param.pilot_pitch_KWD * out.imu_gyro_wbib(2);
        
        % integra o erro de altitude e verifica saturação
        out.alt_err_integral = clamp( out.alt_err_integral + dt*out.alt_err, param.pilot_pitch_KI_sat );

        %% AUTOPILOTO HORIZONTAL: Controle PD com referência do Guiamento        
        out.pilot_yaw = + param.pilot_yaw_KP * (out.guidance_yaw) ...
                        - param.pilot_yaw_KD * out.imu_gyro_wbib(3);

                          
        %% AUTOPILOTOS: combinação dos comandos para deflexão de profundores
        % alocação e saturação dos comandos individuais (4 + 8 + 8 = 20 graus)
        out.pilot_roll =  clamp(out.pilot_roll,  4);
        out.pilot_pitch = clamp(out.pilot_pitch, 8);
        out.pilot_yaw =   clamp(out.pilot_yaw,   8);
                    
        % combinação de sinais dos comandos de manobras 
        out.actuator_delta_ref(1) = ( +out.pilot_yaw -out.pilot_pitch -out.pilot_roll );
        out.actuator_delta_ref(2) = ( -out.pilot_yaw -out.pilot_pitch -out.pilot_roll );
        out.actuator_delta_ref(3) = ( -out.pilot_yaw +out.pilot_pitch -out.pilot_roll );
        out.actuator_delta_ref(4) = ( +out.pilot_yaw +out.pilot_pitch -out.pilot_roll );
        
        % aplicação de batente de software nos profundores:
        out.actuator_delta_ref = clamp( out.actuator_delta_ref, param.actuator_soft_limit );
        
        % atraso de travamento de profundores
        if t < param.actuator_unlock_delay            
            out.actuator_delta_ref = [0,0,0,0]';     % travamento
            out.alt_err_integral = 0;                % não acumula erros
        end
    end
   
end
