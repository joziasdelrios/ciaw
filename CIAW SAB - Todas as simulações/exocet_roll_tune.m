function exocet_roll_tune(roll_kp, roll_kd)

    Xcm = -2.940;
    Jxx = 10.5;
    Jyy = 1490;
    Jzz = 1490;
    Va = 310;

    [air] = sab_air_simple(0, 0);
    qinf = 0.5 * air.density * Va^2;

    Lref = 1.28;
    Sref = 0.554;
    Xref = -2.500;
    
    Cyb = -4.99*2;        % =CNA @ alpha=10deg, Mach=0.90 
    Cyr = 3.42*2;         %  CYR @ alpha= 0deg, Mach=0.90
    Cyd = -0.2269/10;     %  CY  @ alpha= 0deg, Mach=0.90 (YAW10) (opposite) (degrees)
    Cza = Cyb; 
    Czq = -Cyr;          % signal convention
    Czd = Cyd;           % (degrees)
    Cld = 0.0385/10;
    Clp = -0.65;
    Cma = -2.836*2;      % CMA @ alpha=average, Mach=0.90 (NOMINAL) 
    Cmq = -19.34*4;      % CMQ @ alpha=average, Mach=0.90 (NOMINAL)
    Cmd = 0.5605/10*2;   % (degrees)
    Cnb = -Cma;          % =CLNB
    Cnr = Cmq;
    Cnd = Cmd;           % (degrees)

    % parâmetros do modelo Massa-Mola-Amortecedor
    m = Jxx/(qinf*Sref*Lref*Cld);
    b = -(Clp/Cld)*Lref/(2*Va);
    k = 0;

    % posição angular de referência, inicial e velocidade angular inicial
    x_ref = 0;
    x0 = 30*pi/180;
    v0 = 60*pi/180;

    % sem controle intergral
    imax = 0;
    ki = 0;

    % ganho Kp e Kd
    kp = roll_kp;
    kd = roll_kd;
    
    % saturação do controlador: 4 graus (não converte para radianos)
    sat = 4.0;
    
    % constante de tempo de atuação:
    tc = 0.050;
    
    sab_pid_tune(x_ref, x0, v0, m,b,k, kp,ki,kd, imax, sat, tc);
