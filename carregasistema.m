function [DBAR, DLIN] = carregasistema(sistema, Pb)
%CARREGASISTEMA realiza a leitura dos dados
%Dados de barra (DBAR) e dados de linha (DLIN)

%% Inicialização dos dados
% O programa suporta um sistema de até 100 barras e 300 linhas
DBAR = zeros(100, 11);  % Matriz de dados de barra
DLIN = zeros(300, 11);  % Matriz de dados de linha

%% Leitura
% Abre arquivo de leitura
file = fopen(sistema);

% Lê primeiras linhas até a string DBAR ser encontrada
tline = 0;
while ~strcmp(tline,'DBAR')
    tline = fgetl(file);
end

% Leitura dos dados de barra
% Leitura da linha de cabeçalho
fgetl(file);

cont_lidas = 0;  % Contador de barras lidas
while 1
    % Atualiza contador de barras lidas
    cont_lidas = cont_lidas+1;
    
    % Verifica número máximo de barras permitido pelo programa
    if cont_lidas > 100
        disp('oops...')
        disp('Número máximo de barras excedido!!!!')
    end
    
    % Leitura das linhas
    tline = fgetl(file);
    
    % Critério de parada
    % A leitura dos dados de barras é finalizada quando a string
    % '99999' é encontrada
    if strcmp(tline,'99999')
        cont_lidas = cont_lidas-1;
        break
    end
    
    % Número da barra
    DBAR(cont_lidas, 1) = str2double(tline(1:5));
    
    % Tipo de barra
    DBAR(cont_lidas, 2) = str2double(tline(8));
    
    % Tensão (Dividida por 1000)
    DBAR(cont_lidas, 3) = str2double(tline(25:28))/1000;
    
    % Ângulo (transformado de graus para radianos)
    DBAR(cont_lidas, 4) = degtorad(str2double(tline(29:32)));
    %(em radiano dividido por 100)
%    DBAR(cont_lidas, 4) = str2double(tline(29:32))/100;
    % Pg (transformado em pu)
    DBAR(cont_lidas, 5) = str2double(tline(33:37))/Pb;
    
    % Qg (transformado em pu)
    DBAR(cont_lidas, 6) = str2double(tline(38:42))/Pb;
    
    % Qn (transformado em pu)
    DBAR(cont_lidas, 11) = str2double(tline(43:47))/Pb;
    
    % Qm (transformado em pu)
    DBAR(cont_lidas, 7) = str2double(tline(48:52))/Pb;
    
    % Pl (transformado em pu)
    DBAR(cont_lidas, 8) = str2double(tline(59:63))/Pb;
    
    % Ql (transformado em pu)
    DBAR(cont_lidas, 9) = str2double(tline(64:68))/Pb;
    
    % Qsh (transformado em pu)
    DBAR(cont_lidas, 10) = str2double(tline(69:73))/Pb;
    % O valor é invertido pois o dado está em MVAr e deve ser
    % convertido para uma reatância em pu.
end

% Remove NaN's da matriz
% Os NaN's são gerados pelos elementos vazios dos dados de entrada.
DBAR(cont_lidas+1:end, :) = [];
% Remove linhas zeradas (capacidade máxima de barras do programa)
DBAR(isnan(DBAR))=0;

% Leitura dos dados de linha
while ~strcmp(tline,'DLIN')
    tline = fgetl(file);
end
% Leitura da linha de cabeçalho
fgetl(file);


cont_lidas = 0;  % Contador de linhas lidas
while 1
    
    % Atualiza contador de barras lidas
    cont_lidas = cont_lidas+1;
    
    % Verifica número máximo de linhas permitido pelo programa
    if cont_lidas > 300
        disp('oops...')
        disp('Número máximo de linhas excedido!!!!')
    end
    
    % Leitura das linhas
    tline = fgetl(file);
    
    % Critério de parada
    if strcmp(tline,'99999')
        cont_lidas = cont_lidas-1;
        break
    end
    
    
    
    % Barra de
    DLIN(cont_lidas, 1) = str2double(tline(1:5));
    
    % Barra para
    DLIN(cont_lidas, 2) = str2double(tline(11:15));
    
    % r%
    DLIN(cont_lidas, 3) = str2double(tline(21:26))/100;
    
    % x%
    DLIN(cont_lidas, 4) = str2double(tline(27:32))/100;
    
    % 0.5*x_sh%
    DLIN(cont_lidas, 5) = str2double(tline(33:38))/(2*Pb);
    % O valor é invertido pois o dado está em MVAr e deve ser
    % convertido para uma reatância em pu.
    % Dividiu-se o valor da potência em MVAr por dois pois o programa
    % trabalha com o modelo pi da linha de transmissão.
    
    % Tap (Invertido pois o ANAREDE trabalha de t:1)
    DLIN(cont_lidas, 6) = 1/str2double(tline(39:43));
    
    % TapMAX (Invertido pois o ANAREDE trabalha de t:1)
    DLIN(cont_lidas, 8) = 1/str2double(tline(44:48));
    
    % TapMIN (Invertido pois o ANAREDE trabalha de t:1)
    DLIN(cont_lidas, 7) = 1/str2double(tline(49:53));
    
    % Elimina Tap de LT's
    if isnan(DLIN(cont_lidas, 6))
        DLIN(cont_lidas, 6:8) = ones(3,1); 
    end
    
    % Phs (transformado de graus para radianos)
    DLIN(cont_lidas, 9) = degtorad(str2double(tline(54:58)));
    
    % Cap (transformado em pu)
    DLIN(cont_lidas, 10) = str2double(tline(65:68))/Pb;
    
    % Bc (Barra controlada - TAP)
    DLIN(cont_lidas, 11) = str2double(tline(59:64));
end

% Remove NaN's da matriz e linhas zeradas
DLIN(cont_lidas+1:end, :) = [];
% Remove linhas zeradas (capacidade máxima de linhas do programa)
DLIN(isnan(DLIN))=0;
%
%
fclose('all');
end

