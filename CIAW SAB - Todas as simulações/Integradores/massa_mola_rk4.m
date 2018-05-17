% Simulação Massa-Mola por Método de Runge-Kutta ordem 4
% Jozias del Rios - 31.mai.2016

function massa_mola_rk4

clear all;

% valores dos parâmetros
global m b k;
m = 2;
b = 0;
k = 5;

% valores das condições iniciais
v0 = 0;
x0 = 1;

% Simulação por 10 segundos, passo fixo 10ms
t_inicial = 0;
t_final = 10;
dt = 0.01;

% preenchendo condições iniciais
i = 1;
x = [x0];
v = [v0];
t = [t_inicial];

% parâmetros do Movimento Harmônico real
w0 = sqrt(k/m);
d = b/(2 * sqrt(k*m));
w = w0 * sqrt(1 - d^2);

xr = [x0];

while t(i) < t_final
        
    % integrador Runge-Kutta de ordem 4
    kv1 = dt * a( v(i), x(i), t(i) );    
    kp1 = dt * v(i);
    
    kv2 = dt * a( v(i) + kv1/2, x(i) + kp1/2, t(i) + dt/2 );
    kp2 = dt * (v(i) + kv1/2);
    
    kv3 = dt * a( v(i) + kv2/2, x(i) + kp2/2, t(i) + dt/2 );
    kp3 = dt * (v(i) + kv2/2);
    
    kv4 = dt * a( v(i) + kv3, x(i) + kp3, t(i) + dt );
    kp4 = dt * (v(i) + kv3);
    
    % ponderando velocidades e posições
    kv = (kv1 + 2*kv2 + 2*kv3 + kv4) / 6;
    kp = (kp1 + 2*kp2 + 2*kp3 + kp4) / 6;
    
    % integrando velocidade e posição
    v(i + 1) = v(i) + kv;
    x(i + 1) = x(i) + kp;

    % curva real: cálculo da posição real, sem integração
    xr(i + 1) = x0 * exp(-t(i)*w0*d) * cos(w*t(i));    
    
    % próximo passo
    t(i + 1) = t(i) + dt; 
    i = i + 1;
end

% energia total = energia potencial da mola somado com energia cinética
energy = k*x.^2/2 + m*v.^2/2;

figure; 
subplot(3,1,1); plot(t, x); hold all; plot(t, xr); grid on; title('Posição'); 
subplot(3,1,2); plot(t, v); grid on; title('Velocidade');
subplot(3,1,3); plot(t, energy); grid on; title('Energia');

end


% Função local que calcula a aceleração em uma posição, velocidade e tempo
function a = a(v, x, t)
    global m b k;
    
    % cálculo do somatório das forças
    Fext = 0;
    F = -k*x -b*v + Fext;

    % cálculo de aceleração 
    a = F/m;
end

