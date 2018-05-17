function sab_plot_ballistic_6dof_euler(sim,t,s,aux)
    r2d = 180/pi;

    close all;
    set(0,'DefaultFigureWindowStyle','docked');
    set(0,'DefaultLineLineWidth',2);

    figure('Name', 'Trajetórias e Velocidades Lineares');

    subplot(3,1,1); plot( [s.x], -[s.z] );
    grid on; title('Trajetória vertical (XZ)'); xlabel('x alcance [m]'); ylabel('h altura [m]');

    subplot(3,1,2); plot( [s.x], [s.y] ); 
    grid on; title('Trajetória horizontal (XY)'); xlabel('x alcance [m]'); ylabel('y lateral [m]');

    subplot(3,1,3); plot(t,[aux.speed],'r' ); 
    grid on; title('Velocidade total');  xlabel('t [s]'); ylabel('[m/s]');

    figure('Name', 'Orientação e Velocidades Angulares');        

    subplot(3,2,1); plot(t, [s.phi]*r2d,'r'); 
    grid on; title('\phi Ângulo de Rolamento (roll)'); xlabel('t [s]'); ylabel('[\circ]'); 

    subplot(3,2,2); plot(t, [s.p]*r2d,'r' ); 
    grid on; title('p: Velocidade de Rolamento (roll)'); xlabel('t [s]'); ylabel('[\circ/s]'); 

    subplot(3,2,3); plot(t, [s.theta]*r2d,'g' ); 
    grid on; title('\theta: Ângulo de Arfagem (pitch)'); xlabel('t [s]'); ylabel('[\circ]');  
    subplot(3,2,4); plot(t, [s.q]*r2d,'g' ); 
    grid on; title('q: Velocidade de Arfagem (pitch)'); xlabel('t [s]'); ylabel('[\circ/s]');

    subplot(3,2,5); plot(t, [s.psi]*r2d,'b'); 
    grid on; title('\psi: Ângulo de Guinada (yaw)'); xlabel('t [s]'); ylabel('[\circ]'); 

    subplot(3,2,6); plot(t, [s.r]*r2d,'b' ); 
    grid on; title('r: Velocidade de Guinada (yaw)'); xlabel('t [s]'); ylabel('[\circ/s]');
end