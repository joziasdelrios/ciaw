% Simulação MonteCarlo Balística de Bomba em 6DOF
% Jozias del Rios - 30.jul.2016

%{
    montecarlo_balistico_6dof(10);
%} 

function montecarlo_balistico_6dof( histories )
    % semente de números aleatória fixada:
    rng(1000);

    % posição fixa do alvo em relação ao ponto nominal de lançamento:
    xtarget = 3780;
    ytarget = 0;

    % preparação da figura de trajetórias
    set(0,'DefaultFigureWindowStyle','docked');
    set(0,'DefaultLineLineWidth',2);    
    figure('Name', 'Trajetórias Monte-Carlo Balístico 6DOF');
    trajvert = subplot(2,1,1);
    trajhoriz = subplot(2,1,2);

    %% LOOP de simulação Monte-Carlo:
    fprintf('\nIniciando a simulação Monte-Carlo...\n');
    for i = 1:histories
        fprintf('História %d/%d ============================================\n',i,histories);

        % condições iniciais de lançamento da bomba:
        psi = 0.5 + 0.3*randn;
        vel = 300 + randn*5;
        fprintf('Lançamento: psi=%3.2fdeg  velocidade=%4.2f m/s\n', psi, vel);
        
        % Executa a simulação da bomba, sem plotar
        [~,~,s,~,~,~] = balistico_6dof( [0, 0, -10000*0.3048]', [vel*0.5144, 0, 0]', [0,2,psi]'*(pi/180), [0,0,0]', 3);
        
        % plota trajetória vertical deste lançamento, com cor aleatória
        c = rand(1,3);
        subplot(trajvert); hold on; plot( sx([s.pbii]), -sz([s.pbii]), 'color', c); 
        grid on; title('alcance x vertical'); xlabel('x alcance [m]'); ylabel('h altura [m]');

        % plota trajetória horizontal deste lançamento
        subplot(trajhoriz); hold on;  plot( sx([s.pbii]), sy([s.pbii]), 'color', c); 
        grid on; title('alcance x lateral'); xlabel('x alcance [m]'); ylabel('y lateral [m]');
        
        % registra este ponto de impacto no solo
        x(i) = s(end).pbii(1);
        y(i) = s(end).pbii(2);        
    end
    fprintf('RESULTADOS ====================================\n');

    %% calculando e apresentando o ponto MIP
    xmip = mean( x );
    ymip = mean( y );    
    fprintf('MIP = (%4.2f, %4.2f)\n', xmip, ymip);
    
    %% calculando o CEP em torno do MIP e em torno do Alvo
    
    for i=1:histories
        dmip(i)    = sqrt( (x(i) - xmip)^2     +   (y(i) - ymip)^2    );
        dtarget(i) = sqrt( (x(i) - xtarget)^2  +   (y(i) - ytarget)^2 );
    end
    
    % a mediana inclui 50% dos valores
    mip_cep    = median( dmip );
    target_cep = median( dtarget );
    
    fprintf('CEP no MIP  = %4.2f\nCEP no Alvo = %4.2f\n', mip_cep, target_cep);
    
    %% gráfico dos pontos de impacto, Alvo, MIP e CEP
    figure('Name', 'Pontos de Impacto, MIP e CEP');
    plot( x, y,             'r*', 'MarkerSize', 5, 'DisplayName','Ponto de impacto' ); hold on;    
    plot( xtarget, ytarget, 'bo', 'MarkerSize',15, 'DisplayName','ALVO' );
    plot( xmip, ymip,       'rs', 'MarkerSize',15, 'DisplayName','MIP' );
    legend('show');
    axis equal;
    
    % Círculo do CEP em torno do MIP
    pos = [xmip-mip_cep, ymip-mip_cep, mip_cep*2, mip_cep*2];
    rectangle('Position', pos, 'Curvature',[1 1], 'EdgeColor', 'r');
    
    % Círculo do CEP em torno do Alvo
    pos = [xtarget-target_cep, ytarget-target_cep, target_cep*2, target_cep*2];
    rectangle('Position', pos, 'Curvature',[1 1], 'EdgeColor', 'b');
    
end
        
