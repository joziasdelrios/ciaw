% Simulação física da disciplina Simulação e Controle de Artefato Bélico
% ITA - Instituto Tecnológico de Aeronáutica
% IAE - Instituto de Aeronáutica e Espaço
% CIAW - Centro de Instrução Almirante Wandenkolk
%
% Autor: Jozias del Rios - Cap Eng <joziasdelrios@gmail.com>
% Rev. 2 [21 nov 2017]
%
% Chamada:
%       [t, s, aux, out, tout] = sab_sim(sim)
%
% Saídas:
%       t       vetor dos tempos de cada passo (contínuo)
%       s       vetor com estrutura de estados
%       aux     vetor com estrutura de saídas auxiliares
%       out     vetor com estrutura de saídas discretas
%       tout    vetor com o tempo das saídas discretas
%
% Entradas: fornecer estrutura 'sim' preenchida com:
%   sim.method  opções de integradores:
%       "Euler"     Euler de primeira ordem passo fixo
%       "RK4"       Runge-Kutta de quarta ordem passo fixo (default)
%       "RKF45"     Runge-Kutta-Fehlberg adaptativo
%
%   sim.tmax    tempo máximo de simulação (default 10sec)
%   sim.dt      passo de tempo fixo da simulação (default 10ms)
%   sim.dtmax   máximo passo de tempo (default 100ms)
%   sim.dtmin   mínimo passo de tempo (default 1ps)
%   sim.aux     (saída opcional) estrutura inicial de saídas auxiliares
%   sim.param   (entrada opcional) estrutura de parâmetros constantes
%   sim.freq    (opcional) frequência de chamada da função discreta
%   sim.out     (saída opcional) estrutura inicial de saídas discretas
%
%   sim.dsdt = @dsdt    função para calcular a derivada dos estados:
%       [ds, aux] = dsdt( t, dt, s, dsp, auxp, param, out ) 
%           t       tempo atual
%           s       estrutura do estado atual
%           dsp     estrutura da derivada dos estados do passo anterior
%           auxp    estrutura das saídas auxiliares do passo anterior
%           param   estrutura dos parâmetros constantes
%           aux     estrutura a ser preenchida com as saídas auxiliares
%           ds      estrutura a ser preenchida com a derivada dos estados
%           out     estrutura com as saídas discretas atuais
%
%   sim.fstop = @fstop      (opcional) testar o critério de parada:
%       stop = fstop(t, s, ds, aux, param)
%           stop    retorne true ou false quanto ao desejo de interromper
%
%   sim.adjust = @adjust    (opcional) ajustar os estados após cada passo
%       [s] = adjust(t, s, ds, aux, param)
%
%   sim.fstep = @fstep      (opcional) função para ajuste do passo
%       [scale,repeat] = fstep(s, serr)
%           serr    estrutura com o erro absoluto estimado em cada estado
%           scale   fator que multiplica o passo de tempo
%           repeat  variável booleana que indica se deve repetir o passo
%
%   sim.fdiscrete = @fdiscrete          (opcional) função discreta periódica
%       [out] = fdiscrete( t, dt, s, ds, aux, param, outp, z )
%           outp    estrutura das saídas discretas anteriores
%           z       número do passo discreto atual
%
% Autor: Jozias del Rios - 21-nov-2017

function [t,s,aux,out,tout] = sab_sim(sim)
    
    % lista dos nomes dos campos da estrutura de estados
    sfn = fieldnames(sim.s);
    nfn = 1:numel(sfn);
       
    % preencher os campos não utilizados com valores default
    if ~isfield(sim, 'method');     sim.method='RK4';   end
    if ~isfield(sim, 'aux');        sim.aux=0;          end
    if ~isfield(sim, 'param');      sim.param=0;        end
    if ~isfield(sim, 'out');        sim.out=0;          end
    if ~isfield(sim, 'dt');         sim.dt = 0.01;      end
    if ~isfield(sim, 'dtmin');      sim.dtmin = 1e-12;  end
    if ~isfield(sim, 'dtmax');      sim.dtmax = 0.1;    end
    if ~isfield(sim, 'tmax');       sim.tmax = 10;      end
    
    % passo de tempo discreto
    if ~isfield(sim, 'freq');
        zdt = Inf;
    else
        zdt = 1/sim.freq;
    end
    

    % vetor de estados, com prealocação:
    sprealloc = 1 + ceil(sim.tmax / sim.dt);
    zprealloc = 1 + ceil(sim.tmax / zdt);
    
    %   s(i)   é o vetor de estrutura de estados
    %  ds(i)   é o vetor de estrutura de derivadas dos estados
    % aux(i)   é o vetor de estrutura de saídas auxiliares
    %     dt   é o valor atual do passo
    s(1:sprealloc) = sim.s;    
    [ds(1), aux(1)] = sim.dsdt( 0, sim.dt, sim.s, sim.s, sim.aux, sim.param, sim.out );
    aux(2:sprealloc) = aux(1);
    t(1:sprealloc) = 0;
        
    %out(1:zprealloc) = sim.out;
    out(1) = sim.out;
    tout(1:zprealloc) = 0;
    %tout(1) = 0;
    lastout = sim.out;
    i = 1;
    z = 1;

    % RKF45 Butcher tableau: maths.cnam.fr/IMG/pdf/RungeKuttaFehlbergProof.pdf
    a = [1/4, 3/32, 9/32, 1932/2197, -7200/2197, 7296/2197, 439/216, ...
         -8, 3680/513, -845/4104, -8/27, 2, -3544/2565, 1859/4104, -11/40];
    b = [1/4, 3/8, 12/13, 1, 1/2];
    c = [25/216, 0, 1408/2565, 2197/4101, -1/5];
    d = [16/135, 0, 6656/12825, 28561/56430, -9/50, 2/55];    
    
    
    % valor proposto do passo de tempo
    dtp = min(sim.dt, zdt);
    
    stopping = false;
    
    % inicia o loop de simulação global...
    while t(i) < sim.tmax
        
        % valor atual do passo de tempo, a ser usado na integração
        dt = dtp;
                
        % fazer um passo de simulação contínua...
        if sim.method == "RK4"

            % integrador de Runge-Kutta ordem 4 passo fixo
            ht = dt/2;
            [ds1,        ~] = sim.dsdt( t(i),    dt, s(i), ds, aux(i), sim.param, lastout ); for j=nfn, n=sfn{j}; su.(n) = s(i).(n) + ht*ds1.(n); end
            [ds2,        ~] = sim.dsdt( t(i)+ht, dt, su,   ds, aux(i), sim.param, lastout ); for j=nfn, n=sfn{j}; su.(n) = s(i).(n) + ht*ds2.(n); end
            [ds3,        ~] = sim.dsdt( t(i)+ht, dt, su,   ds, aux(i), sim.param, lastout ); for j=nfn, n=sfn{j}; su.(n) = s(i).(n) + dt*ds3.(n); end
            [ds4, aux(i+1)] = sim.dsdt( t(i)+dt, dt, su,   ds, aux(i), sim.param, lastout );

            % novo vetor de estados
            for j = 1:numel(sfn)
                n = sfn{j}; 
                ds.(n) = 1/6 * ( ds1.(n) + 2*ds2.(n) + 2*ds3.(n) + ds4.(n) );
                s(i+1).(n) = s(i).(n) + dt*ds.(n);
            end
            
        elseif sim.method == "Euler"
            % integrador de Euler
            [ds, aux(i+1)] = sim.dsdt( t(i), dt, s(i), ds, aux(i), sim.param, lastout );
            
            % novo vetor de estados
            for j = 1:numel(sfn)
                n = sfn{j}; 
                s(i+1).(n) = s(i).(n) + dt*ds.(n);
            end
            
        elseif sim.method == "RKF45"
            % integrador de Runge-Kutta-Fehlberg ordem 4 com erro ordem 5
            [ds1,        ~] = sim.dsdt( t(i),         dt, s(i), ds, aux(i), sim.param, lastout ); for j=nfn, n=sfn{j}; su.(n) = s(i).(n) + dt*( a(1)*ds1.(n) ); end
            [ds2,        ~] = sim.dsdt( t(i)+dt*b(1), dt, su,   ds, aux(i), sim.param, lastout ); for j=nfn, n=sfn{j}; su.(n) = s(i).(n) + dt*( a(2)*ds1.(n) + a(3)*ds2.(n) ); end
            [ds3,        ~] = sim.dsdt( t(i)+dt*b(2), dt, su,   ds, aux(i), sim.param, lastout ); for j=nfn, n=sfn{j}; su.(n) = s(i).(n) + dt*( a(4)*ds1.(n) + a(5)*ds2.(n) +  a(6)*ds3.(n)); end
            [ds4,        ~] = sim.dsdt( t(i)+dt*b(3), dt, su,   ds, aux(i), sim.param, lastout ); for j=nfn, n=sfn{j}; su.(n) = s(i).(n) + dt*( a(7)*ds1.(n) + a(8)*ds2.(n) +  a(9)*ds3.(n) + a(10)*ds4.(n)); end
            [ds5,        ~] = sim.dsdt( t(i)+dt*b(4), dt, su,   ds, aux(i), sim.param, lastout ); for j=nfn, n=sfn{j}; su.(n) = s(i).(n) + dt*( a(11)*ds1.(n) + a(12)*ds2.(n) +  a(13)*ds3.(n) + a(14)*ds4.(n) + a(15)*ds5.(n)); end
            [ds6, aux(i+1)] = sim.dsdt( t(i)+dt*b(5), dt, su,   ds, aux(i), sim.param, lastout ); 
 
            for j=nfn
                n = sfn{j}; 
                ds4th.(n) = c(1)*ds1.(n) + c(2)*ds2.(n) + c(3)*ds3.(n) + c(4)*ds4.(n) + c(5)*ds5.(n);
                ds.(n) =    d(1)*ds1.(n) + d(2)*ds2.(n) + d(3)*ds3.(n) + d(4)*ds4.(n) + d(5)*ds5.(n) + d(6)*ds6.(n);

                serr.(n) = dt * abs( ds4th.(n) - ds.(n) );
                s(i+1).(n) = s(i).(n) + dt*ds.(n);
            end
            
        else
            warning('SABSIM: Method name error');
        end
        
        % atualiza o tempo do próximo passo, usando o passo da integração
        t(i+1) = t(i) + dt;
                
        % testa erro para ajustar e aceitar ou recusar o passo dado
        if sim.method == "RKF45" && stopping == false
            if isfield(sim, 'fstep')
                % função customizada do usuário para análise do passo
                [scale, repeat] = sim.fstep(s(i+1), serr);
                dtp = dt * scale;
                if dtp < sim.dtmin; dtp = sim.dtmin; end
                if dtp > sim.dtmax; dtp = sim.dtmax; end
                if dtp > zdt;       dtp = zdt; end
                if scale < 1.0 && repeat == true; continue; end
            else
                 % ajuste automático do passo com base no erro absoluto
                 %{
                 err = 0;
                 for j=nfn
                     n = sfn{j};
                     err = err + abs( ds4th.(n) - ds.(n) );
                 end
                 if err > eps(err)
                     est = ( eps(err) * dt / err
                 eps
                 %}
            end 
        end
        
        % ajuste da estrutura de estados entre os passos
        if isfield(sim, 'adjust')
            s(i+1) = sim.adjust( t(i+1), s(i+1), ds, aux(i+1), sim.param );
        end
        
        % se entrar condição de parada...
        if isfield(sim, 'fstop')
            if sim.fstop( t(i), s(i), ds, aux(i), sim.param ) == true
                stopping = true;
                dtp = dtp/2;          % reduz o passo pela metade
                if dtp > sim.dtmin
                    % refaz o passo se o passo não for inferior ao mínimo
                    continue        
                else
                    % se reduziu demais, interrompe a simulação
                    i = i+1;
                    break;          
                end
            end
        end
        
        % passo aceito, verifica o tempo da simulação discreta
        if (t(i) + dt) > ( (z-1) * zdt)
            if isfield(sim,'fdiscrete')
                % chamada da função discreta
                lastout = sim.fdiscrete( t(i), zdt, s(i), ds, aux(i+1), sim.param, out, z ); 
                if z > 1
                    out(z) = lastout;
                else
                    out = lastout;
                    out(2:zprealloc) = lastout;
                end
                tout(z) = t(i);
                z = z + 1;
            end
        end
        
        % passo aceito, avança o índice
        i = i+1;    
    end
    
    if i < sprealloc
        s = s(1:i);
        t = t(1:i);
        aux = aux(1:i);
    end    
    
    if z < zprealloc
        zlast = max(1, z-1);
        out = out(1:zlast);
        tout = tout(1:zlast);
    end
    
end

    