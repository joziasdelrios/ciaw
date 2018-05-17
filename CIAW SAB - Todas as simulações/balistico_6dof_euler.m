% Simula��o Bal�stica de Bomba em 6DOF 
% Jozias del Rios - 11.jun.2016

%{
  [sim,t,s,aux] = balistico_6dof_euler( 0, 0, -10000*0.3048, 300*0.5144,0,0,  0,+2*pi/180,0, 0,0,0 );
  sab_plot_ballistic_6dof_euler(sim,t,s,aux);
%}

function [sim,t,s,aux] = balistico_6dof_euler( x0,y0,z0, vx0,vy0,vz0, psi0,theta0,phi0, p0,q0,r0 )
   
    tic;    
    addpath('sab');    
        
    %% par�metros do corpo da bomba
    sim.param.mass = 490;
    sim.param.Jxx = 8.2;
    sim.param.Jyy = 155;
    sim.param.Jzz = 155;
    sim.param.Xcg = 1.3;
       
    %% estados iniciais de posi��o
    sim.s.x = x0;
    sim.s.y = y0;
    sim.s.z = z0;

    %% convertendo velocidades iniciais para ref. corpo
    cr = cos(phi0);   sr = sin(phi0);
    cp = cos(theta0); sp = sin(theta0);
    cy = cos(psi0);   sy = sin(psi0);
    sim.s.u = vx0*( cy*cp )          + vy0*( sy*cp )          + vz0*(-sp );
    sim.s.v = vx0*( cy*sp*sr-sy*cr ) + vy0*( sy*cp*sr+cy*cr ) + vz0*( cp*sr );
    sim.s.w = vx0*( cy*sp*cr+sy*sr ) + vy0*( sy*sp*cr-cy*sr ) + vz0*( cp*cr );
    
    %% estados iniciais de orienta��o e velocidade angular
    sim.s.psi = psi0;
    sim.s.theta = theta0;
    sim.s.phi = phi0;
    sim.s.p = p0;
    sim.s.q = q0;
    sim.s.r = r0;
    
    %% tempo de simula��o: 5min
    sim.tmax = 300;   
    
    %% fun��es de simula��o
    sim.dsdt = @dsdt;
    sim.fstop = @fstop;
    
    %% simulador
    sim.method = "RK4";
    [t,s,aux,~,~] = sab_sim(sim);
    
    %% resultado final da simula��o
    r2d = 180/pi;
    fprintf('t=%5.2f. imp=(%4.2f, %4.2f). phi=%3.2f, theta=%3.2f, psi=%3.2f\n', ...
        t(end), s(end).x, s(end).y, s(end).phi*r2d, s(end).theta*r2d, s(end).psi*r2d );
    
    toc;

    %% fun��o de c�lculo derivada dos estados
    function [ds, aux] = dsdt( t, dt, s, dsp, auxp, param, out )

        %% obtendo cossenos e senos dos �ngulos de Euler
        cr = cos(s.phi);   sr = sin(s.phi);
        cp = cos(s.theta); sp = sin(s.theta);
        cy = cos(s.psi);   sy = sin(s.psi);

        %% par�metros atmosf�ricos locais (altitude atual, deltaISA=zero)
        aux.air = sab_air_simple( -s.z, 0);
        
        % velocidade aerodin�mica, press�o din�mica, Mach, AoA, AoS, ...
        aux.speed = sqrt(s.u^2 + s.w^2);
        aux.aoa = atan2(s.w, s.u);
        aux.aos = atan2(s.v, s.u);
        aux.dynpress = aux.air.density * aux.speed^2 / 2;
        aux.Mach = aux.speed / aux.air.sound_speed;    

        %% c�lculo do somat�rio das for�as e torques externos

        % for�a peso 
        Fweight.x = - param.mass * aux.air.gravity * sp;
        Fweight.y = + param.mass * aux.air.gravity * cp*sr;
        Fweight.z = + param.mass * aux.air.gravity * cp*cr;

        % esfor�os aerodin�micos: estima��o das for�as e momentos
        [aero] = sab_bomb_aero(aux.speed, aux.Mach, aux.dynpress, aux.aoa, aux.aos, [s.p, s.q, s.r]');
            
        % for�as totais sobre o corpo
        Fx = Fweight.x + aero.Fx;
        Fy = Fweight.y + aero.Fy;
        Fz = Fweight.z + aero.Fz;
        
        % torques totais sobre o corpo
        Mx = aero.Mx;
        My = aero.My + (aero.Xref - param.Xcg)*Fz;    % ajuste ref_aer -> cm (produto vetorial)
        Mz = aero.Mz + (aero.Xref - param.Xcg)*Fy;    % ajuste ref_aer -> cm (produto vetorial)

        %% equa��es de din�mica: integrando a velocidade linear
        ds.x = s.u*( cy*cp ) + s.v*( cy*sp*sr-sy*cr ) + s.w*( cy*sp*cr+sy*sr );
        ds.y = s.u*( sy*cp ) + s.v*( sy*cp*sr+cy*cr ) + s.w*( sy*sp*cr-cy*sr );
        ds.z = s.u*( -sp )   + s.v*( cp*sr )          + s.w*( cp*cr );

        %% equa��es de din�mica: integrando a acelera��o linear
        ds.u = Fx/param.mass + s.r*s.v - s.q*s.w;
        ds.v = Fy/param.mass + s.p*s.w - s.r*s.u;
        ds.w = Fz/param.mass + s.q*s.u - s.p*s.v;

        %% equa��es de din�mica: integrando a velocidade angular
        ds.phi   = s.p + (s.q * sr + s.r * cr)*sp/cp;
        ds.theta = s.q * cr - s.r * sr;
        ds.psi   = (s.q * sr + s.r * cr)/cp;  

        %% equa��es de din�mica: integrando a acelera��o angular
        ds.p = (Mx + s.q*s.r*(param.Jyy - param.Jzz)) / param.Jxx;
        ds.q = (My + s.p*s.r*(param.Jzz - param.Jxx)) / param.Jyy;
        ds.r = (Mz + s.p*s.q*(param.Jxx - param.Jyy)) / param.Jzz;

    end    
   
    %% crit�rio de parada do simulador
    function stop = fstop(t, s, ds, aux, param)        
        % interrompa a simula��o se colidir com o solo
        if s.z >= 0 
            stop = true;
        else
            stop = false;
        end
    end
    
end

