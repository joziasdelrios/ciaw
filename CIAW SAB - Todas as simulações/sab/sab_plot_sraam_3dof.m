function sab_plot_sraam_3dof( sim,t,s,aux )
    r2d = 180/pi;    
    c = colormap('gray');
    
    close all;
    set(0,'DefaultFigureWindowStyle','docked');
    set(0,'DefaultLineLineWidth',2);
    
    %% informação de miss-distance
    [dm,idm] = min([aux.miss_distance]);
    mx = [s.mx];
    my = [s.my];
    tx = [s.tx];
    ty = [s.ty];
    psi = [s.mpsi];
    sp = sin(psi);
    cp = cos(psi);
        
    %% PLOT TRAJETÓRIA
    figure('Name', 'Trajetória');
    rt = 2;
    
    plot( my(idm), mx(idm), 'b*', 'MarkerSize', 8); hold on;
    plot( ty(idm), tx(idm), 'r*', 'MarkerSize', 8);
    
    for i = 1:rt*floor(t(end))
        ti = i/sim.dt/rt;
        line( [my(ti), ty(ti)], [mx(ti), tx(ti)], 'Color',[0.8 0.8 0.8], 'LineStyle','-', 'LineWidth', 1 );        
    end
    
    % line([x1,x2],[y1,y2])
    line( [my(idm), ty(idm)], [mx(idm), tx(idm)], 'Color','black', 'LineWidth', 3 );

    
    sz = 100;
    for i = 1:rt*floor(t(end))
        ti = i/sim.dt/rt;
        line( [my(ti)-sz*sp(ti), my(ti)+sz*sp(ti)], [mx(ti)-sz*cp(ti), mx(ti)+sz*cp(ti)], 'Color',[0.5 0.5 0.5], 'LineWidth', 2 );        
    end        
    
    plot( my, mx, 'b', 'LineWidth', 3 ); 
    plot( ty, tx, 'r', 'LineWidth', 3 );    

    for i = 1:rt*floor(t(end))
        ti = i/sim.dt/rt;
        plot( my(ti), mx(ti), 'b*', 'MarkerSize', 8);
        plot( ty(ti), tx(ti), 'r*', 'MarkerSize', 8);
    end
    
    grid on; title('Trajetória'); xlabel('y [m]'); ylabel('x [m]');
    daspect([1 1 1]);
    
    %% PLOT VELOCIDADES E ACELERAÇÕES LINEARES
    figure('Name', 'Velocidades e Acelerações');

    ax(1) = subplot(2,2,1); 
    plot( t, [aux.speed], t,[aux.mvx], t, [aux.mvy] ); 
    grid on; title('Velocidades');  xlabel('t [s]'); ylabel('[m/s]');
    legend('Velocity', 'vbii_x', 'vbii_y');

    ax(2) = subplot(2,2,2);
    plot(t, [aux.Mach],'k', 'LineWidth',3  ); 
    grid on; xlabel('t [s]'); ylabel('[Mach]');

    ax(3) = subplot(2,2,3);
    plot(t, [aux.thrust].*0.001,'k', 'LineWidth',3  ); 
    grid on; title('Empuxo do motor-foguete'); xlabel('t [s]'); ylabel('[kN]');

    ax(4) = subplot(2,2,4);
    plot( t, [aux.nav_axb],'r', t, [aux.nav_ayb],'b' ); 
    grid on; title('Acelerações'); xlabel('t [s]'); ylabel('aceleração [m/s^2]');
    legend('Axial', 'LATAX right' );
    linkaxes(ax(1:4),'x');
    
    %% PLOT ANGULARES
    figure('Name', 'Orientação e Velocidades Angulares');
    
    ax(1) = subplot(3,1,1);
    plot(t, [s.mpsi].*r2d ); hold on;
    plot(t, [aux.mgamma].*r2d );
    plot(t, [aux.los].*r2d );    
    grid on; title('Ângulos planares no referencial inercial'); xlabel('t [s]'); ylabel('[\circ]');  
    legend('\psi proa', '\gamma rumo','\lambda LoS' );
    
    ax(2) = subplot(3,1,2);
    plot(t, [s.mr].*r2d ); 
    grid on; title('Velocidade angular (yaw speed)'); xlabel('t [s]'); ylabel('[\circ]/s');     
    
    ax(3) = subplot(3,1,3);
    plot(t, [aux.aos].*r2d ); hold on;
    plot(t, [aux.seeker_azimuth].*r2d );
    plot(t, [aux.delta] );
    grid on; title('Ângulos no referencial corpo'); xlabel('t [s]'); ylabel('[\circ]');  
    legend('\beta AoS', '\psi_{AD} azimute', '\delta atuador');
    
    linkaxes(ax(1:3),'x');

end