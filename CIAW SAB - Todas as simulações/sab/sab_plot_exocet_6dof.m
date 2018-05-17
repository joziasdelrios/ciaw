function sab_plot_exocet_6dof( sim,t,s,aux,out,tout )
    r2d = 180/pi;    
    c = lines;   
    
    close all;
    set(0,'DefaultFigureWindowStyle','docked');
    set(0,'DefaultLineLineWidth',2);
    
    %% PLOT TRAJETÓRIAS PLANOS
    figure('Name', 'Trajetórias nos Planos');
    
    [dm,idm] = min([aux.miss_distance]);
    
    subplot(3,1,1); 
    plot( sx([s.pbii]), -sz([s.pbii]), 'Color', c(1,:) ); hold on;
    plot( sx([s.ptii]), -sz([s.ptii]), 'Color', c(2,:) ); 
    plot( sx([out.nav_pbii]), -sz([out.nav_pbii]), ':', 'Color', c(3,:) );
    plot( s(idm).pbii(1), -s(idm).pbii(3), 'o', 'MarkerSize', 5, 'Color', c(1,:) );  
    plot( s(idm).ptii(1), -s(idm).ptii(3), 'o', 'MarkerSize', 5, 'Color', c(2,:) );    
    grid on; title('Trajetória vertical (Vista Lateral) (XZ)'); xlabel('x alcance [m]'); ylabel('h altura [m]');
    legend('Míssil', 'Alvo', 'Navegação');

    subplot(3,1,2); 
    plot( sx([s.pbii]), sy([s.pbii]), 'Color', c(1,:) ); hold on;
    plot( sx([s.ptii]), sy([s.ptii]), 'Color', c(2,:) ); 
    plot( sx([out.nav_pbii]), sy([out.nav_pbii]), ':', 'Color', c(3,:) );
    plot( s(idm).pbii(1), s(idm).pbii(2), 'o', 'MarkerSize', 5, 'Color', c(1,:));  
    plot( s(idm).ptii(1), s(idm).ptii(2), 'o', 'MarkerSize', 5, 'Color', c(2,:));    
    grid on; title('Trajetória horizontal (Vista de Topo) (XY)'); xlabel('x alcance [m]'); ylabel('y lateral [m]');
    set(gca,'YDir','Reverse'); legend('Míssil', 'Alvo', 'Navegação');
    
    subplot(3,1,3); 
    plot( sy([s.pbii]), -sz([s.pbii]), 'Color', c(1,:) );  hold on;
    plot( sy([s.ptii]), -sz([s.ptii]), 'Color', c(2,:) );
    plot( sy([out.nav_pbii]), -sz([out.nav_pbii]), ':', 'Color', c(3,:) );
    plot( s(idm).pbii(2), -s(idm).pbii(3), 'o', 'MarkerSize', 5,  'Color', c(1,:));  
    plot( s(idm).ptii(2), -s(idm).ptii(3), 'o', 'MarkerSize', 5,  'Color', c(2,:));    
    grid on; title('Trajetória (Vista Frontal) (YZ)'); xlabel('y lateral [m]'); ylabel('z altura [m]');
    legend('Míssil', 'Alvo', 'Navegação');
    
    %% PLOT 3D
    figure('Name', 'Trajetórias 3D');
    plot3( sx([s.pbii]), sy([s.pbii]), -sz([s.pbii]),  'Color', c(1,:), 'LineWidth',3 ); hold on;
    plot3( sx([s.ptii]), sy([s.ptii]), -sz([s.ptii]),  'Color', c(2,:), 'LineWidth',3 ); 
    plot3( s(idm).pbii(1), s(idm).pbii(2), -s(idm).pbii(3), 'o', 'MarkerSize', 5,  'Color', c(1,:));  
    plot3( s(idm).ptii(1), s(idm).ptii(2), -s(idm).ptii(3), 'o', 'MarkerSize', 5,  'Color', c(2,:)); 
    %zlim([-10000, +20]);
    title('Trajetória tridimensional'); xlabel('x alcance [m]'); ylabel('y lateral [m]'); zlabel('h altura [m]');
    legend('Míssil', 'Alvo');
    grid on; 
    axis vis3d; 
    axis equal;
    set(gca,'YDir','Reverse'); 
    set(gca,'Projection','perspective'); 
    set(gca,'DataAspectRatio', [1,1,1/20]);
    set(gca,'FontSize', 14);
    
    %% PLOT VELOCIDADES E ACELERAÇÕES LINEARES
    figure('Name', 'Velocidades e Acelerações');
    
    subplot(2,2,1); 
    plot( t, [aux.velocity], t,sx([aux.vbii]), t, sy([aux.vbii]), t, sz([aux.vbii]) ); 
    grid on; title('Velocidades');  xlabel('t [s]'); ylabel('[m/s]');
    legend('Velocity', 'vbii_x', 'vbii_y', 'vbii_z');

    subplot(2,2,2); plot(t, [aux.Mach],'k', 'LineWidth',3  ); 
    grid on; xlabel('t [s]'); ylabel('[Mach]');

    subplot(2,2,3); plot(t, sx([aux.thrust_force]).*0.001,'k', 'LineWidth',3  ); 
    grid on; title('Empuxo do motor-foguete'); xlabel('t [s]'); ylabel('[kN]');

    subplot(2,2,4); plot( t, sx([aux.abib]),'r', t, sy([aux.abib]),'g', t,sz([aux.abib]),'b' ); 
    grid on; title('Acelerações'); xlabel('t [s]'); ylabel('aceleração [m/s^2]');
    legend('Axial', 'LATAX right', 'LATAX down');
        
    %% PLOT ANGULARES
    figure('Name', 'Orientação e Velocidades Angulares');
    
    subplot(3,2,1); 
    plot(t, sx([aux.euler]).*r2d, 'Color', c(1,:) ); hold on;
    plot(tout, sx([out.nav_euler]).*r2d, ':', 'Color', c(5,:)); 
    grid on; title('\phi Ângulo de Rolamento (roll angle)'); xlabel('t [s]'); ylabel('[\circ]'); 
    legend('\phi roll', '\phi roll nav');
    
    subplot(3,2,2); 
    plot(t, sx([s.wbib]).*r2d, 'Color', c(1,:) ); 
    grid on; title('p: Velocidade de Rolamento (roll speed)'); xlabel('t [s]'); ylabel('[\circ/s]'); 

    subplot(3,2,3);
    plot(t, sy([aux.euler]).*r2d, 'Color', c(2,:)); hold on;
    plot(t, [aux.gamma_pitch].*r2d, 'Color', c(4,:));
    plot(tout, sy([out.nav_euler]).*r2d, ':', 'Color', c(5,:)); 
    grid on; title('\theta: Ângulo de Arfagem (pitch angle)'); xlabel('t [s]'); ylabel('[\circ]');  
    legend('\theta pitch', '\gamma rumo vert.','\theta nav');
    subplot(3,2,4); 
    plot(t, sy([s.wbib]).*r2d, 'Color', c(2,:)); 
    grid on; title('q: Velocidade de Arfagem (pitch speed)'); xlabel('t [s]'); ylabel('[\circ/s]');
    
    subplot(3,2,5); 
    plot(t, sz([aux.euler]).*r2d, 'Color', c(3,:) ); hold on;
    plot(t, [aux.gamma_yaw].*r2d, 'Color', c(4,:)); 
    plot(tout, sz([out.nav_euler]).*r2d, ':', 'Color', c(5,:)); 
    grid on; title('\psi: Ângulo de Guinada (yaw angle) e de Rumo'); xlabel('t [s]'); ylabel('[\circ]'); 
    legend('\psi yaw proa', '\gamma rumo horiz.', '\psi yaw nav');
    subplot(3,2,6);
    plot(t, sz([s.wbib]).*r2d, 'Color', c(3,:)); 
    grid on; title('r: Velocidade de Guinada (yaw speed)'); xlabel('t [s]'); ylabel('[\circ/s]');    
    
    %% NAVEGAÇÃO
    figure('Name','Guiamento, Navegação e Controle');
    
    subplot(3,2,1); 
    plot(t,[aux.aoa].*r2d,  t,[aux.aos].*r2d ); 
    grid on; title('Ângulos de Ataque'); xlabel('t [s]'); ylabel('[\circ]'); 
    legend('\alpha AoA','\beta AoS');
    
    subplot(3,2,2);     
    plot(t,si([s.actuator_delta],1)); hold on;
    plot(t,si([s.actuator_delta],2));
    plot(t,si([s.actuator_delta],3));
    plot(t,si([s.actuator_delta],4));
    grid on; title('Comandos de manobra nas empenas'); xlabel('t [s]'); ylabel('[\circ]'); 
    legend('\delta_1', '\delta_2', '\delta_3', '\delta_4');
    
    subplot(3,2,3); 
    plot(tout, [out.seeker_look_angle].*r2d, 'Color', c(3,:) ); hold on;
    %plot(tout, [out.seeker_roll_angle].*r2d, 'Color', c(4,:) ); 
    plot(tout, [out.seeker_azimuth].*r2d, 'Color', c(1,:) ); 
    plot(tout, [out.seeker_elevation].*r2d, 'Color', c(2,:) ); 
    grid on; title('Ângulos de Azimute, Elevação, Look e Rolamento do Autodiretor (Seeker)'); 
    xlabel('t [s]'); ylabel('[\circ]'); 
    legend('\lambda_{AD} Look-Angle', '\psi_{AD} Azimuth','\theta_{AD} Elevation' );
    ylim([-90, +90]);
    
    subplot(3,2,4);
    plot(tout,[out.pilot_roll] ); hold on;
    plot(tout,[out.pilot_pitch] );
    plot(tout,[out.pilot_yaw] ); 
    grid on; title('Comandos de manobra do Autopiloto'); xlabel('t [s]'); ylabel('[\circ]'); 
    legend('\phi_{CMD} roll', '\theta_{CMD} pitch', '\psi_{CMD} yaw');
    
    subplot(3,2,[5 6]);
    plot(t, -sz([s.pbii]) ); hold on;    
    plot(tout, [out.alt_err_integral] );
    grid on; title('Controle de altitude'); xlabel('t [s]'); ylabel('[m]'); 
    legend('altitude(t)', 'integral do erro de altitude' );

    %% Turbulência
    figure('Name', 'Turbulência');
    turb = [out.turb];
    plot(tout, sx([turb.vwii]),'r'); hold on;
    plot(tout, sy([turb.vwii]),'g');
    plot(tout, sz([turb.vwii]),'b');
    grid on; title('Turbulência'); xlabel('t [s]'); ylabel('[m/s]');

end