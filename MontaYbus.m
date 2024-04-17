function Ybus = MontaYbus(NBAR, NLIN, SB, EB, R_E, X_E, Xsh_E, tap, phs, Xsh)
%MONTAYBUS constrói a Ybus

% Inicialização da Ybus - todos os elementos nulos
Ybus = zeros(NBAR, NBAR);
for ilin = 1:NLIN
    de = SB(ilin); % Barra 'de' da linha i 
    pa = EB(ilin); % Barra 'para' da linha i
    a = tap(ilin); % Tap do transformador i
    fase = phs(ilin); % Defasagem do defasador i 
    ykm = 1/(R_E(ilin)+1j*X_E(ilin)); % Cálculo da admitância série da linha i
    % Zera admitâncias infinitas
    % A não existência das linhas geram uma impedância infinita que 
    % é traduzido como uma admitância nula
    if ykm == inf
        ykm = 0;
    end
    % Calcula a susceptância shunt da linha
    b = 1j*Xsh_E(ilin);
    
    % Montagem por inspeção da Ybus
    % Elementos em kk recebem o somatório de todas
    % as admitâncias conectadas à barra k
    Ybus(de, de) = Ybus(de, de) + a^2*ykm + b;
    Ybus(pa, pa) = Ybus(pa, pa) + ykm + b;
    % Elementos em km recebem o negativo da admitância 
    % conectada entre k e m
    Ybus(de, pa) = Ybus(de, pa) - ykm*a*exp(-1j*fase);
    Ybus(pa, de) = Ybus(pa, de) - ykm*a*exp(+1j*fase);
end

% Laço para adicionar elementos shunt na barra
for kbar = 1:NBAR
    % Cálculo da susceptância shunt na barra
    b = 1j*Xsh(kbar);
    % Elementos em kk recebem o somatório de todas
    % as admitâncias conectadas à barra k
    Ybus(kbar, kbar) = Ybus(kbar, kbar) + b;
end
%
%
end

