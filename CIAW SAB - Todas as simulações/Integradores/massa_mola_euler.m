% Simula��o Massa-Mola por M�todo de Euler
% Jozias del Rios - 31.mai.2016

clear all;

% valores dos par�metros
m = 2;
b = 0;
k = 5;

% valores das condi��es iniciais
v0 = 0;
x0 = 1;

% Simula��o por 10 segundos, passo fixo 10ms
t_inicial = 0;
t_final = 10;
dt = 0.01;

% preenchendo condi��es iniciais
i = 1;
x = [x0];
v = [v0];
t = [t_inicial];

% par�metros do Movimento Harm�nico real
w0 = sqrt(k/m);
d = b/(2 * sqrt(k*m));
w = w0 * sqrt(1 - d^2);

xr = [x0];

while t(i) < t_final
    
    % c�lculo do somat�rio das for�as
    Fext = 0;
    F = -k*x(i) -b*v(i) + Fext;
    
    % c�lculo de acelera��o
    a = F/m;
   
    % integra��o por Euler obtem a velocidade
    v(i + 1) = v(i) + dt * a;
    
    % integra��o por Euler obtem a posi��o
    x(i + 1) = x(i) + dt * v(i);
       
    % curva real: c�lculo da posi��o real, sem integra��o
    xr(i + 1) = x0 * exp(-t(i)*w0*d) * cos(w*t(i));    
    
    % pr�ximo passo de tempo
    t(i + 1) = t(i) + dt; 
    i = i + 1;
end

% energia total = energia potencial da mola somado com energia cin�tica
energy = k*x.^2/2 + m*v.^2/2;

figure; 
subplot(3,1,1); plot(t, x); hold all; plot(t, xr); grid on; title('Posi��o'); 
subplot(3,1,2); plot(t, v); grid on; title('Velocidade');
subplot(3,1,3); plot(t, energy); grid on; title('Energia');

