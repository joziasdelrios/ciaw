% Esta fun��o estuda movimentos harm�nicos (oscila��es) simples e 
% amortecidas atrav�s de um simula��o massa-mola-amortecedor sem for�a externa
% 
% O objetivo � a resolu��o de equa��es diferenciais ordin�rias (EDO) com 
% aplica��o f�sica (cinem�tica e din�mica) utilizando v�rias op��es de
% integradores: RKF45 (passo adaptativo), RK4 (quarta ordem) e Euler.
%
% Autor: Jozias del Rios - Rev 1 - 04.dez.2017

function massa_mola
    
    addpath('sab');
    tic;
    
    %% passo inicial da simula��o e tempo limite
    sim.dt = 0.01;      % passo de tempo 10ms
    sim.tmax = 25;      % tempo m�ximo 25 segundos

    %% valores das condi��es iniciais dos estados a serem integrados
    sim.s.x = 1;        % posi��o
    sim.s.v = 0;        % velocidade
    
    %% fun��es de c�lculo da derivada e de parada do simulador
    sim.dsdt = @dsdt;
    sim.fstop = @fstop;
       
    %% par�metros dos modelos do sistema massa-mola-amortecedor
    sim.param.m = 2;    % massa do corpo de 2 kg
    sim.param.b = 0;    % amortecedor nulo: 0 kg/s
    sim.param.k = 5;    % constante de mola 5 N/m
    
    %% simular com algum m�todo: "Euler", "RK4" ou "RKF45"
    sim.method = "RK4";
    [t,s,aux,~,~] = sab_sim(sim);
    
    %% plot dos gr�ficos massa-mola
    set(0,'DefaultFigureWindowStyle','docked');
    set(groot,'defaultLineLineWidth',2);

    figure;
    subplot(3,1,1); plot(t, [s.x]); grid on; title('Posi��o x');    xlabel('[s]'); ylabel('[m]');
    subplot(3,1,2); plot(t, [s.v]); grid on; title('Velocidade v'); xlabel('[s]'); ylabel('[m/s]');
    subplot(3,1,3); plot(t, [aux.energy]); grid on; title('Energia [J]'); xlabel('[s]'); ylabel('[J]');    
    
    toc;
    
    %----------------------------------------------------------------------
    % Fun��o que calcula a derivada temporal dos estados:
    % t      tempo atual
    % dt     passo de tempo atual
    % s      estrutura de estados atual
    % dsp    estrutura de derivada temporal dos estados do passo anterior
    % auxp   estrutura de sa�das auxiliares do passo anterior
    % param  estrutura de par�metros constantes
    % out    estrutura com as sa�das discretas atuais
    % ds     sa�da: estrutura de derivada dos estados
    % aux    sa�da: estrutura de sa�das auxiliares
    function [ds, aux] = dsdt(t, dt, s, dsp, auxp, param, out)
        
        % Somat�rio das for�as no corpo:
        %   For�a da mola: -k*x
        %   For�a do amortecedor: -b*v
        F = - param.k * s.x - param.b * s.v;

        % c�lculo da din�mica de acelera��o, 2a Lei de Newton
        aux.a = F / param.m;

        % calcula a energia do sistema em cada instante de tempo
        %   Energia potencial na mola: k*x^2/2
        %   Energia cin�tica do corpo: m*v^2/2
        aux.energy = param.k * s.x^2/2   +   param.m * s.v^2/2;

        % derivada da velocidade � a acelera��o (calculada)
        ds.v = aux.a;
        % derivada da posi��o � a velocidade (estado)
        ds.x = s.v;
    end

    %----------------------------------------------------------------------
    % Fun��o que analisa o crit�rio de parada da simula��o
    % t      tempo atual
    % s      estrutura de estados atual
    % ds     estrutura de derivada dos estados calculada pela fun��o dsdt
    % aux    estrutura de sa�das auxiliares calculada pela fun��o dsdt
    % param  estrutura de par�metros constantes
    function stop = fstop(t, s, ds, aux, param)
        
        % interrompe a simula��o se a energia do sistema for muito baixa
        if aux.energy < 0.01
            stop = true;
        else
            stop = false;
        end
    end

end
