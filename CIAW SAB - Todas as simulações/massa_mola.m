% Esta função estuda movimentos harmônicos (oscilações) simples e 
% amortecidas através de um simulação massa-mola-amortecedor sem força externa
% 
% O objetivo é a resolução de equações diferenciais ordinárias (EDO) com 
% aplicação física (cinemática e dinâmica) utilizando várias opções de
% integradores: RKF45 (passo adaptativo), RK4 (quarta ordem) e Euler.
%
% Autor: Jozias del Rios - Rev 1 - 04.dez.2017

function massa_mola
    
    addpath('sab');
    tic;
    
    %% passo inicial da simulação e tempo limite
    sim.dt = 0.01;      % passo de tempo 10ms
    sim.tmax = 25;      % tempo máximo 25 segundos

    %% valores das condições iniciais dos estados a serem integrados
    sim.s.x = 1;        % posição
    sim.s.v = 0;        % velocidade
    
    %% funções de cálculo da derivada e de parada do simulador
    sim.dsdt = @dsdt;
    sim.fstop = @fstop;
       
    %% parâmetros dos modelos do sistema massa-mola-amortecedor
    sim.param.m = 2;    % massa do corpo de 2 kg
    sim.param.b = 0;    % amortecedor nulo: 0 kg/s
    sim.param.k = 5;    % constante de mola 5 N/m
    
    %% simular com algum método: "Euler", "RK4" ou "RKF45"
    sim.method = "RK4";
    [t,s,aux,~,~] = sab_sim(sim);
    
    %% plot dos gráficos massa-mola
    set(0,'DefaultFigureWindowStyle','docked');
    set(groot,'defaultLineLineWidth',2);

    figure;
    subplot(3,1,1); plot(t, [s.x]); grid on; title('Posição x');    xlabel('[s]'); ylabel('[m]');
    subplot(3,1,2); plot(t, [s.v]); grid on; title('Velocidade v'); xlabel('[s]'); ylabel('[m/s]');
    subplot(3,1,3); plot(t, [aux.energy]); grid on; title('Energia [J]'); xlabel('[s]'); ylabel('[J]');    
    
    toc;
    
    %----------------------------------------------------------------------
    % Função que calcula a derivada temporal dos estados:
    % t      tempo atual
    % dt     passo de tempo atual
    % s      estrutura de estados atual
    % dsp    estrutura de derivada temporal dos estados do passo anterior
    % auxp   estrutura de saídas auxiliares do passo anterior
    % param  estrutura de parâmetros constantes
    % out    estrutura com as saídas discretas atuais
    % ds     saída: estrutura de derivada dos estados
    % aux    saída: estrutura de saídas auxiliares
    function [ds, aux] = dsdt(t, dt, s, dsp, auxp, param, out)
        
        % Somatório das forças no corpo:
        %   Força da mola: -k*x
        %   Força do amortecedor: -b*v
        F = - param.k * s.x - param.b * s.v;

        % cálculo da dinâmica de aceleração, 2a Lei de Newton
        aux.a = F / param.m;

        % calcula a energia do sistema em cada instante de tempo
        %   Energia potencial na mola: k*x^2/2
        %   Energia cinética do corpo: m*v^2/2
        aux.energy = param.k * s.x^2/2   +   param.m * s.v^2/2;

        % derivada da velocidade é a aceleração (calculada)
        ds.v = aux.a;
        % derivada da posição é a velocidade (estado)
        ds.x = s.v;
    end

    %----------------------------------------------------------------------
    % Função que analisa o critério de parada da simulação
    % t      tempo atual
    % s      estrutura de estados atual
    % ds     estrutura de derivada dos estados calculada pela função dsdt
    % aux    estrutura de saídas auxiliares calculada pela função dsdt
    % param  estrutura de parâmetros constantes
    function stop = fstop(t, s, ds, aux, param)
        
        % interrompe a simulação se a energia do sistema for muito baixa
        if aux.energy < 0.01
            stop = true;
        else
            stop = false;
        end
    end

end
