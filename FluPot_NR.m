%% An�lise de Redes El�tricas I
% Fluxo de pot�ncia - Newton-Raphson
% Arthur Neves de Paula

clc
clear
close all
%
%% Par�metros de simula��o
PAR = struct(...
    'MAX_ITR', 100,... % M�ximo de itera��es
    'EPS', 0.001,...  % Crit�rio de converg�ncia
    'BIGN', 1e5);  % Big Number
%
%% Leitura do arquivo de dados
% Escolha do sistema
%sistema = '4barrasJLUIZ.pwf';
sistema = '6barrasJLUIZ.pwf';
%sistema = '14barrasIEEE.pwf';
%sistema = 'IEEE 24 Barras Bc3.pwf';
%sistema = 'IEEE 24 Barras.pwf';
%sistema = 'IEEE 24 Barras-MODIFICADO.pwf';
% sistema = 'Teste3barra.pwf';
%sistema = 'Teste2barra.pwf';

% Leitura dos dados (acessa a fun��o que l� o arquivo de dados)
disp('Leitura dos dados...')
Pb = 100;  % Pot�ncia base
[DBAR, DLIN] = carregasistema(['Casos\', sistema], Pb);

% Calcula n�mero de barras e linhas do sistema
[NBAR, ~] = size(DBAR);
[NLIN, ~] = size(DLIN);

% Dados de linha
SB = DLIN(:, 1); % Vetor de barras emissoras
EB = DLIN(:, 2); % Vetor de barras receptoras
R_E = DLIN(:, 3); % Vetor de resist�ncias s�rie de linha
X_E = DLIN(:, 4); % Vetor de reat�ncias s�rie de linha
Xsh_E = DLIN(:, 5); % Vetor de reat�ncia shunt de linha
tap = DLIN(:, 6); % Vetor de taps de transformadores
tapN = DLIN(:, 7); % Vetor de taps m�nimos de transformadores
tapM = DLIN(:, 8); % Vetor de taps m�ximos de transformadores
phs = DLIN(:, 9); % Vetor de defasagem de defasadores
Bc = DLIN(:, 11); % Vetor de Barras controladas

% Dados de barra
numB = DBAR(:, 1); % Vetor de n�mero de barras
tipB = DBAR(:, 2); % Vetor de tipo de barras 0-PQ 1-PV 2-SW
VB = DBAR(:, 3); % Vetor de tens�es de barras
TB = DBAR(:, 4); % Vetor de �ngulos nodais
PG = DBAR(:, 5); % Vetor de pot�ncias ativas geradas
QG = DBAR(:, 6); % Vetor de pot�ncias reativas geradas
QN = DBAR(:, 11); % Vetor gera��es m�nimas de pot�ncias reativas
QM = DBAR(:, 7); % Vetor gera��es m�ximas de pot�ncias reativas
PD = DBAR(:, 8); % Vetor de pot�ncias ativas demandadas
QD = DBAR(:, 9); % Vetor de pot�ncias reativas demandadas
Xsh = DBAR(:, 10); % Vetor de suscept�ncias shunt de barra

% tipos de barra iniciais
tipB0 = tipB;

% Pot�ncias especificadas
Pesp = PG - PD;
Qesp = QG - QD;

% Inicializa converg�ncia
conv = 0;
%
%% Montagem da Ybus
disp('Montagem da Ybus...')
Ybus = MontaYbus(NBAR, NLIN, SB, EB, R_E, X_E, Xsh_E, tap, phs, Xsh);

%% C�lculo do Fluxo de Pot�ncia - Newton Raphson
% Inicializando FP
V = VB; % Tens�es iniciais
T = TB; % �ngulo nodais

% Submatrizes zeradas para composi��o da Jacobiana
Za = zeros(NBAR, NBAR);
disp('M�todo de Newton-Raphson...')
for itr = 1:PAR.MAX_ITR
    
    % Condut�ncia e suscept�ncia s�rie
    G = real(Ybus);
    B = imag(Ybus);
    
    % Inicializa as submatrizes
    H = zeros(NBAR, NBAR);
    N = zeros(NBAR, NBAR);
    M = zeros(NBAR, NBAR);
    L = zeros(NBAR, NBAR);
    
    % Inicializa vetores de pot�ncias injetadas e res�duos
    Pcalc = zeros(NBAR,1);
    Qcalc = zeros(NBAR,1);
    deltaP = zeros(NBAR,1);
    deltaQ = zeros(NBAR,1);
    deltaV = zeros(NBAR,1);
    
    % C�lculo das pot�ncias
    for kbar = 1:NBAR
        for mbar = 1:NBAR
            if Ybus(kbar, mbar) == 0
                continue
            end
            Gkm = G(kbar, mbar);
            Bkm = B(kbar, mbar);
            Tkm = T(kbar)-T(mbar);
            
            Pcalc(kbar)=Pcalc(kbar)+V(kbar)*V(mbar)*(Gkm*cos(Tkm)+Bkm*sin(Tkm));
            Qcalc(kbar)=Qcalc(kbar)+V(kbar)*V(mbar)*(Gkm*sin(Tkm)-Bkm*cos(Tkm));
            
            % Submatrizes do Jacobiano (fora da diagonal)
            if kbar ~= mbar
                H(kbar, mbar) = +V(kbar)*V(mbar)*(Gkm*sin(Tkm)-Bkm*cos(Tkm));
                N(kbar, mbar) = +V(kbar)*(Gkm*cos(Tkm)+Bkm*sin(Tkm));
                M(kbar, mbar) = -V(kbar)*V(mbar)*(Gkm*cos(Tkm)+Bkm*sin(Tkm));
                L(kbar, mbar) = +V(kbar)*(Gkm*sin(Tkm)-Bkm*cos(Tkm));
            end
        end
            
        %Submatrizes do Jacobiano (diagonal)
        H(kbar, kbar) = -Qcalc(kbar)-V(kbar)^2*B(kbar,kbar);
        N(kbar, kbar) = (Pcalc(kbar)+V(kbar)^2*G(kbar,kbar))/V(kbar);
        M(kbar, kbar) = +Pcalc(kbar)-V(kbar)^2*G(kbar,kbar);
        L(kbar, kbar) = (Qcalc(kbar)-V(kbar)^2*B(kbar,kbar))/V(kbar);
        
        % Res�duos
        deltaP(kbar)=Pesp(kbar)-Pcalc(kbar);
        deltaQ(kbar)=Qesp(kbar)-Qcalc(kbar);
    end
    
    % Verifica converg�ncia
    if conv == 1
        disp('Convergiu!')
        break
    end
    
    %Construindo o vetor de res�duos
    gX=[deltaP;deltaQ;deltaV];
    
    % Montagem da Jacobiana
    J = [H  N; 
         M  L];
    
    % Big number
    for kbar = 1:NBAR
        if tipB(kbar) == 0
            J(2*NBAR+kbar, 2*NBAR+kbar) = PAR.BIGN;
        elseif tipB(kbar) == 1
            J(NBAR+kbar, NBAR+kbar) = PAR.BIGN;
            J(2*NBAR+kbar, 2*NBAR+kbar) = PAR.BIGN;
        elseif tipB(kbar) == 2
            J(kbar, kbar) = PAR.BIGN;
            J(NBAR+kbar, NBAR+kbar) = PAR.BIGN;
            J(2*NBAR+kbar, 2*NBAR+kbar) = PAR.BIGN;
        end
    end
    
    % C�lculo do deltaX
    deltaX = J\gX;
    
    % Atualiza tens�es nodais
    T = T+deltaX(1:NBAR);
    V = V+deltaX(NBAR+1:2*NBAR);

    % Crit�rio de parada
    if max(abs(deltaX)) < PAR.EPS
        conv = 1;
    end
end

%Diverg�ncia do m�todo
if itr == PAR.MAX_ITR+1
    disp('Oops... o n�mero m�ximo de itera��es foi atingido');
    return
end

%Subsistema 2
disp('Gerando relat�rios...')
[Pger, Qger, Pflu, Qflu, Qsh, Pperdas, Qperdas, QLIM__] =...
    CalculaSS2(V, T, Ybus, Xsh, Xsh_E, NLIN, Pcalc, Qcalc,...
    SB, EB, PD, QD, NBAR, tipB, tipB0);

%Imprimindo dados de barra em txt
NB_ = num2str(numB,'%03d');
V_____ = num2str(V,'%5.4f');
T_____ = num2str(radtodeg(T),'%+06.1f');
Pger__ = num2str(Pger*Pb,'%+06.1f');
Qger__ = num2str(Qger*Pb,'%+06.1f');
MVArSH = num2str(Qsh*Pb,'%+06.1f');

% Removendo valores nulos
for kbar = 1:NBAR
    if strcmp(Pger__(kbar, :),'+000.0') || strcmp(Pger__(kbar, :),'-000.0')
        Pger__(kbar, :) = '      ';
    end
    if strcmp(Qger__(kbar, :),'+000.0') || strcmp(Qger__(kbar, :),'-000.0')
        Qger__(kbar, :) = '      ';
    end
    if strcmp(MVArSH(kbar, :),'+000.0') || strcmp(MVArSH(kbar, :),'-000.0')
        MVArSH(kbar, :) = '      ';
    end
end

t1 = table(NB_,V_____,T_____,Pger__,Qger__,QLIM__,MVArSH);
writetable(t1,'Resultados\RBAR.txt','Delimiter','tab')

%Imprimindo dados de linha em txt
SB____ = num2str(SB,'%06d');
EB____ = num2str(EB,'%06d');
PES___ = num2str(Pflu(:,1)*Pb,'%+06.1f');
QES___ = num2str(Qflu(:,1)*Pb,'%+06.1f');
PSE___ = num2str(Pflu(:,2)*Pb,'%+06.1f');
QSE___ = num2str(Qflu(:,2)*Pb,'%+06.1f');
Pperda = num2str(Pperdas*Pb,'%+06.1f');
Qperda = num2str(Qperdas*Pb,'%+06.1f');

t2 = table(SB____,EB____,PES___,QES___,PSE___,QSE___,Pperda,Qperda);
writetable(t2,'Resultados\RLIN.txt','Delimiter','tab')

disp('Fim do programa.')