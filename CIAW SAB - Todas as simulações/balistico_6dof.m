% Simula��o Bal�stica de Bomba em 6DOF com Matriz DCM e Quaternions
% Jozias del Rios - 09.set.2016

% pbii = posi��o do corpo em rela��o ao referencial inercial, em coordenadas inerciais
% vbii = velocidades lineares do corpo em rela��o ao referencial inercial, em coordenadas inerciais
% euler = �ngulos de Euler da orienta��o em rela��o ao referencial inercial, sendo: phi (roll), theta (pitch), psi (yaw)
% wbib = velocidades angulares do corpo em rela��o ao referencial inercial, em coordenadas do corpo
% weather = n�vel de 0 � 7 de severidade do clima (vento, rajada, turbul�ncia)

%{
 [sim,t,s,aux,out,tout] = balistico_6dof( [0, 0, -10000*0.3048]', [300*0.5144, 0, 0]', [0,2,0]'*pi/180, [0,0,0]', 3);
 sab_plot_ballistic_6dof( sim,t,s,aux,out,tout );
%}


function [sim,t,s,aux,out,tout] = balistico_6dof( pbii0, vbii0, euler0, wbib0, weather )

    addpath('sab');
    tic;
        
    %% par�metros do corpo da bomba
    sim.param.mass = 490;
    sim.param.inertia = diag([8.2 155 155]);            % matriz tensor de in�rcia (Jxx Jyy Jzz)
    sim.param.inv_inertia = inv(sim.param.inertia);     % inverso do tensor de in�rcia
    sim.param.cm = [1.300, 0.0, 0.0]';                  % posi��odo centro de massa
    sim.param.vwii = sab_wind(weather);                 % vento constante aleat�rio
    sim.param.weather = weather;                        % intensidade de perturba��o clim�tica
    
    %% preparando processos discretos: 
    sim.freq = 100;                       % frequ�ncia de atualiza��o dos processos discretos
    sim.out.turb = sab_turbulence_init;   % inicializa as vari�veis de turbul�ncia
    
    %% preparando os estados nas condi��es iniciais
    sim.s.pbii = pbii0;                  % posi��o inicial, no referencial inercial, coordenadas inerciais
    sim.s.qbi  = sab_euler2quat(euler0); % converter os �ngulos de Euler iniciais para o quaternion de orienta��o do corpo em rela��o ao referencial inercial
    sim.s.wbib = wbib0;                  % velocidade angular inicial, no referencial inercial, coordenadas do corpo
    Rbi0 = sab_euler2dcm(euler0);        % converter os �ngulos de Euler iniciais para a matriz de rota��o de orienta��o do corpo em rela��o ao referencial inercial
    sim.s.vbib = Rbi0 * vbii0;           % converter a velocidade linear inicial no referencial inercial de coordenadas inerciais (vbii) para corpo (vbib)
    
    %% passo da simula��o (10ms) e tempo m�ximo (5min)
    sim.dt   = 0.01;
    sim.tmax = 300;
    
    %% fun��es de opera��o e simula��o
    sim.dsdt = @dsdt;
    sim.fstop = @fstop;
    sim.adjust = @fadjust;
    sim.fdiscrete = @fdiscrete;
    sim.method = "RK4";
    [t,s,aux,out,tout] = sab_sim(sim);
       
    %% resultado final da simula��o
    r2d = 180/pi;
    fprintf('t=%5.2fs \timp=(%4.2f, %4.2f)  ', t(end), s(end).pbii(1), s(end).pbii(2));        
    fprintf('\tphi=%3.2f theta=%3.2f psi=%3.2f  ', aux(end).euler(1)*r2d, aux(end).euler(2)*r2d, aux(end).euler(3)*r2d  );
    fprintf('\tvwii=(%4.2f,%4.2f)\n', sim.param.vwii(1), sim.param.vwii(2));
    toc;
        
    %% fun��o de c�lculo derivada dos estados
    function [ds, aux] = dsdt( t, dt, s, dsp, auxp, param, out )
        %% c�lculos iniciais de cinem�tica
        % converte de quaternion para matriz de rota��o (DCM)
        Rbi = sab_quat2dcm( s.qbi );
        
        % converte a velocidade linear de ref. corpo para ref. inercial
        aux.vbii = Rbi' * s.vbib;
        
        % obt�m os �ngulos de Euler da orienta��o atual para os gr�ficos
        aux.euler = sab_quat2euler( s.qbi );

        %% c�lculos atmosf�ricos
        % par�metros atmosf�ricos locais (altitude atual, deltaISA=zero)
        h = -s.pbii(3);
        aux.air = sab_air_simple( h, 0);
        
        % velocidade aerodin�mica no sistema de coordenadas do corpo
        vwib = Rbi * (out.turb.vwii + param.vwii);
        aux.vaib = s.vbib - vwib;
        
        % velocidade aerodin�mica, press�o din�mica, Mach, AoA, AoS, ...
        aux.velocity = norm( s.vbib );          % velocity � a velocidade absoluta do corpo (ref. inercial)
        aux.speed = norm( aux.vaib );           % speed � o TAS (True Air Speed): m�dulo da velocidade relativa ao ar.
        aux.dynpress = aux.air.density * aux.speed^2 / 2;
        aux.Mach = aux.speed / aux.air.sound_speed;
        aux.aoa = atan2(aux.vaib(3), aux.vaib(1));      % �ngulo de ataque (w/u)
        aux.aos = atan2(aux.vaib(2), aux.vaib(1));      % �ngulo de derrapagem (v/u)
        
        %% c�lculo do somat�rio das for�as e torques externos
        
        % for�a peso ocorre na dire��o do eixo +z (vertical inercial)
        %   � necess�ro converter o eixo para as coordenadas do corpo
        weight.force = param.mass * aux.air.gravity * Rbi * [0 0 1]';

        % esfor�os aerodin�micos: fun��o de estima��o das for�as e momentos
        [aero] = sab_bomb_aero(aux.speed, aux.Mach, aux.dynpress, aux.aoa, aux.aos, s.wbib);

        %% somat�rio das for�as e momentos totais
        forces = weight.force + aero.force;
        torques = aero.torque + cross(aero.ref - param.cm, aero.force);
        
        %% equa��es de din�mica
        ds.pbii = Rbi' * s.vbib;
        ds.vbib = forces/param.mass - cross(s.wbib, s.vbib);
        ds.qbi  = (1/2) * sab_quatmult( s.qbi, [0 s.wbib'] );
        ds.wbib = param.inv_inertia * (torques - cross(s.wbib, param.inertia * s.wbib));
        
    end

    %% crit�rio de parada do simulador
    function stop = fstop(t, s, ds, aux, param)        
        
        % interrompa a simula��o se colidir com o solo
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

    %% processos discretos (n�o integrados)
    function [out] = fdiscrete( t, dt, s, ds, aux, param, outp, z )
        
        % segure os valores discretos anteriores que n�o forem alterados...
        out = outp(end);
        
        % atualize o processo de turbul�ncia de Dryden
        out.turb = sab_turbulence( param.vwii, aux.vbii, t, dt, -s.pbii(3), param.weather, outp, z);
    end

end

