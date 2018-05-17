% Simulação Massa-Mola por ODE45
% Jozias del Rios - 31.mai.2016

function massa_mola_ode45

% valores dos parâmetros
global m b k;
m = 2;
b = 0;
k = 5;

% valores das condições iniciais
v0 = 0; 
x0 = 1;

% Simulação por 10 segundos
t_inicial = 0;
t_final = 10;

% Configurando o ODE45 para não fazer passos maiores que 500ms
options = odeset('MaxStep', 0.5 );

% Simulando usando a função f e espaço de estados inicial [x0, v0]
[t, s] = ode45(@f, [t_inicial t_final], [x0 v0], options);

energy = k*s(:,1).^2/2 + m*s(:,2).^2/2;

% parâmetros do Movimento Harmônico real
w0 = sqrt(k/m);
d = b/(2 * sqrt(k*m));
w = w0 * sqrt(1 - d^2);

% calculando a curva real, sem simulação
for i = 1:numel(t)
    xr(i) =  x0 * exp(-t(i)*w0*d) * cos(w*t(i));
end

figure; 
subplot(3,1,1); plot(t, s(:,1)); hold all; plot(t, xr); grid on; title('Posição'); 
subplot(3,1,2); plot(t, s(:,2)); grid on; title('Velocidade');
subplot(3,1,3); plot(t, energy); grid on; title('Energia');
size(s)
end

% Função que calcula a derivada dos estados
function dsdt = f(t, s)
    
    % expande a variável de estados
    x = s(1);
    v = s(2);
    
    % cálculo do somatório das forças
    global m b k;
    Fext = 0;
    F = -k*x -b*v + Fext;
    
    % cálculo de aceleração
    a = F/m;
    
    % preenche a variação da variável de estados
    dsdt = [v; a];
end
  
