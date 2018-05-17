% Simulação Massa-Mola por Método de Euler
% Jozias del Rios - 31.mai.2016

clear all;

% valores dos parâmetros
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
    
    % cálculo do somatório das forças
    Fext = 0;
    F = -k*x(i) -b*v(i) + Fext;
    
    % cálculo de aceleração
    a = F/m;
   
    % integração por Euler obtem a velocidade
    v(i + 1) = v(i) + dt * a;
    
    % integração por Euler obtem a posição
    x(i + 1) = x(i) + dt * v(i);
       
    % curva real: cálculo da posição real, sem integração
    xr(i + 1) = x0 * exp(-t(i)*w0*d) * cos(w*t(i));    
    
    % próximo passo de tempo
    t(i + 1) = t(i) + dt; 
    i = i + 1;
end

% energia total = energia potencial da mola somado com energia cinética
energy = k*x.^2/2 + m*v.^2/2;

figure; 
subplot(3,1,1); plot(t, x); hold all; plot(t, xr); grid on; title('Posição'); 
subplot(3,1,2); plot(t, v); grid on; title('Velocidade');
subplot(3,1,3); plot(t, energy); grid on; title('Energia');

