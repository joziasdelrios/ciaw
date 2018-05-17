function sab_plot_ballistic_3dof( sim,t,s,aux )
    r2d = 180/pi;    
    
    set(0,'DefaultFigureWindowStyle','docked');
    set(0,'DefaultLineLineWidth',2);
    
    figure;
    
    subplot(2,2,1); plot( [s.x], -[s.z] ); 
    grid on; title('Trajetória'); xlabel('x [m]'); ylabel('h [m]'); ylim([0 200-s(1).z]); daspect([1 1 1]);
    
    subplot(2,2,2); plot(t, [aux.vx],'g', t,[aux.vz],'b', t,[aux.speed],'r' ); 
    grid on; title('Velocidades');  xlabel('t [s]'); ylabel('[m/s]');
    legend('vx', 'vz', 'va');
        
    subplot(2,2,3); plot(t, [s.theta]*r2d, 'b', t,[aux.aoa]*r2d,'r' ); 
    grid on; xlabel('t [s]'); ylabel('[\circ]');
    legend('\theta Ângulo de Arfagem (pitch)', '\alpha Ângulo de ataque (AoA)' );
    
    subplot(2,2,4); plot(t, [s.q]*r2d ); 
    grid on; title('Velocidade Angular'); xlabel('t [s]'); ylabel('\omega [\circ/s]'); 

end