% Turbulence Dryden discrete spectrum model 
% according to MIL-F-8785C (Flying Qualities of Piloted Aircraft)
%
% returns only turb.vwii (linear air speed turbulence)
%
% initialize:
%    sim.out.turb = sab_turbulence_init;
%
%  usage as a discrete random process:
%    function [out] = fdiscrete( t, dt, s, ds, aux, param, outp, z )
%       out.turb = sab_turbulence( param.vwii, aux.vbii, t, dt, -s.pbii(3), param.weather, outp, z);
% 

function [turb] = sab_turbulence(vwii, vbii, t, dt, altitude, weather, outp, z)       
    % hold previous values...
    zlast = max(1, z-1);
    turb = outp( zlast ).turb;
    
    % no turbulence? null wind
    if weather == 0
        turb.vwii = [0 0 0]';
        return;
    end    
    
    % 10 Hz / 0.1sec turbulence wind update
    [zi,zm] = sab_discrete_subrate(z, dt, 0.1);
    if zi > 0
        if z > zm
            % interpolate 
            vwiip = outp(z-zi-1).turb.vwiif;
            turb.vwii = vwiip + (turb.vwiif - vwiip) * zi / zm;
        else
            % just hold            
        end
    else
        % update turbulence wind intensities
    
        % convert height to feet and clip above 10ft
        h = max(10, altitude / 0.3048);

        % get windspeed in feet/sec (fps) and wingspan in feet
        windspeed_at_20ft = norm(vwii) / 0.3048;

        % Probability of Exceedance (MIL-HDBK-1797 fig 262 / MIL-F-9490D Table V)
        poe_h = [ 500  1750  3750  7500 15000 25000 35000 45000 55000 65000 75000 80000]';
        poe   = [ 3.2,  2.2,  1.5,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0; ...
                  4.2,  3.6,  3.3,  1.6,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0; ...
                  6.6,  6.9,  7.4,  6.7,  4.6,  2.7,  0.4,  0.0,  0.0,  0.0,  0.0,  0.0; ...
                  8.6,  9.6, 10.6, 10.1,  8.0,  6.6,  5.0,  4.2,  2.7,  0.0,  0.0,  0.0; ...
                 11.8, 13.0, 16.0, 15.1, 11.6,  9.7,  8.1,  8.2,  7.9,  4.9,  3.2,  2.1; ...
                 15.6, 17.6, 23.0, 23.6, 22.1, 20.0, 16.0, 15.1, 12.1,  7.9,  6.2,  5.1; ...
                 18.7, 21.5, 28.4, 30.2, 30.7, 31.0, 25.2, 23.1, 17.5, 10.7,  8.4,  7.2];

        if h <= 1000
            L_u = h / (0.177 + 0.000823*h)^1.2;         % MIL-F-8785C, fig 10
            L_w = h;
            sig_w = 0.1 * windspeed_at_20ft;            % MIL-F-8785C, para 3.7.3.4 (Low-altitude disturbance model)
            sig_u = sig_w / (0.177 + 0.000823*h)^0.4;   % MIL-F-8785C, fig 11

        elseif h <= 2000
            % linear interpolation between low altitude and high altitude models
            L_u = 1000 + (h - 1000)/1000*750;
            L_w = L_u;
            gust = interp2( poe_h, 1:7, poe, h, weather);
            sig_u = 0.1*windspeed_at_20ft + (h-1000)/1000*(gust - 0.1*windspeed_at_20ft);
            sig_w = sig_u;
        else
            % MIL-STD-1797, para Appendix A - 4.9.2 / 2.a
            L_u = 1750;
            L_w = L_u/2;
            sig_u = interp2( poe_h, 1:7, poe, h, weather);
            sig_w = sig_u;
        end

        va = norm( vwii - vbii ) / 0.3048;
        tovtau_u = dt*va/L_u;
        tovtau_w = dt*va/L_w;

        turb.uvw = [0 0 0]';
        puvw = outp(zlast).turb.uvw;
        turb.uvw(1) = (1 -   tovtau_u)*puvw(1) + sig_u*sqrt(2*tovtau_u)*randn;
        turb.uvw(2) = (1 - 2*tovtau_u)*puvw(2) + sig_u*sqrt(4*tovtau_u)*randn;
        turb.uvw(3) = (1 - 2*tovtau_w)*puvw(3) + sig_w*sqrt(4*tovtau_w)*randn;

        % psiw (Wind heading) is the direction the wind is blowing towards
        psiw = atan2( vwii(2), vwii(1) );

        % turn into wind direction
        Rwi = [cos(psiw), sin(psiw), 0; -sin(psiw), cos(psiw), 0; 0, 0, 1];       
        turb.vwiif = 0.3048 * Rwi*turb.uvw;
        
        % current wind turbulence is previous computed turbulence
        turb.vwii = outp(zlast).turb.vwiif;
    end
end
