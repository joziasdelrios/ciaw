

function p = nez_sraam_3dof
    lethality_criteria = 20;

    sraam.mx = 0; sraam.my = 0; 
    sraam.mu = 200;
    sraam.law = "TPN";
    sraam.Kpn = 4;
    
    xsteps = 9;
    ysteps = 6;
    p = zeros(xsteps,ysteps);
    
    for y = 1:ysteps
        for x = 1:xsteps
            sraam.tx = x * 1000;
            sraam.ty = (y-1) * 1000;    
            sraam.tpsi = 180*pi/180;    % head-on --> escaping
            sraam.tu = 300;
            sraam.ta = -10;

            [sim,t,s,aux] = sraam_3dof( sraam );
            dm = min([aux.miss_distance]);            
            fprintf('Target @ (x = %5.2f, y = %5.2f) [m] --> dmiss = %5.2f [m] \n', sraam.tx, sraam.ty, dm);
            
            p(x,y) = dm;
        end
    end

    figure;
    for y = 1:ysteps
        for x = 1:xsteps
            tx = x * 1000;
            ty = (y-1) * 1000;
                       
            if p(x,y) > lethality_criteria                
                % miss
                plot( +ty, +tx, 'rs', 'MarkerSize', 10 ); hold on;
                plot( -ty, +tx, 'rs', 'MarkerSize', 10 );
                t = sprintf('%2.1f', p(x,y) );
                text( +ty-300, +tx+250, t, 'Color','red','FontSize',14, 'Color','red'  );
                text( -ty-300, +tx+250, t, 'Color','red','FontSize',14, 'Color','red'  );
            else
                % hit
                plot( -ty, +tx, 'bx', 'MarkerSize', 10 ); hold on;
                plot( +ty, +tx, 'bx', 'MarkerSize', 10 );
                t = sprintf('%2.1f', p(x,y) );
                text( +ty-300, +tx+250, t, 'Color','blue','FontSize',14, 'Color','blue' );
                text( -ty-300, +tx+250, t, 'Color','blue','FontSize',14, 'Color','blue' );
            end            
        end
    end
    
    ylim([0, (xsteps+1)*1000]);
    xlim([-(ysteps+1)*1000, +(ysteps+1)*1000]);
end

            

                
            

    
    