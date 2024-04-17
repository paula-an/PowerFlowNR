function Ybus = MontaYbus(NBAR, NLIN, SB, EB, R_E, X_E, Xsh_E, tap, phs, Xsh)
%MONTAYBUS constr�i a Ybus

% Inicializa��o da Ybus - todos os elementos nulos
Ybus = zeros(NBAR, NBAR);
for ilin = 1:NLIN
    de = SB(ilin); % Barra 'de' da linha i 
    pa = EB(ilin); % Barra 'para' da linha i
    a = tap(ilin); % Tap do transformador i
    fase = phs(ilin); % Defasagem do defasador i 
    ykm = 1/(R_E(ilin)+1j*X_E(ilin)); % C�lculo da admit�ncia s�rie da linha i
    % Zera admit�ncias infinitas
    % A n�o exist�ncia das linhas geram uma imped�ncia infinita que 
    % � traduzido como uma admit�ncia nula
    if ykm == inf
        ykm = 0;
    end
    % Calcula a suscept�ncia shunt da linha
    b = 1j*Xsh_E(ilin);
    
    % Montagem por inspe��o da Ybus
    % Elementos em kk recebem o somat�rio de todas
    % as admit�ncias conectadas � barra k
    Ybus(de, de) = Ybus(de, de) + a^2*ykm + b;
    Ybus(pa, pa) = Ybus(pa, pa) + ykm + b;
    % Elementos em km recebem o negativo da admit�ncia 
    % conectada entre k e m
    Ybus(de, pa) = Ybus(de, pa) - ykm*a*exp(-1j*fase);
    Ybus(pa, de) = Ybus(pa, de) - ykm*a*exp(+1j*fase);
end

% La�o para adicionar elementos shunt na barra
for kbar = 1:NBAR
    % C�lculo da suscept�ncia shunt na barra
    b = 1j*Xsh(kbar);
    % Elementos em kk recebem o somat�rio de todas
    % as admit�ncias conectadas � barra k
    Ybus(kbar, kbar) = Ybus(kbar, kbar) + b;
end
%
%
end

