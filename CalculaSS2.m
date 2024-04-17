function [Pger, Qger, Pflu, Qflu, Qsh, Pperdas, Qperdas, QLIM__] = CalculaSS2(V, T, Ybus, Xsh, Xsh_E, NLIN, Pcalc, Qcalc, SB, EB, PD, QD, NBAR, tipB, tipB0)
%CALCULASS2 Calcula o subsistema 2 do método gerando os relatórios.

% Condutância e susceptância série
G = real(Ybus);
B = imag(Ybus);

% Cálculo das potências geradas
Pger = Pcalc+PD;
Qger = Qcalc+QD;

% Potência reativa injetada por reatores
Qsh = Xsh.*V.^2;

% Fluxos nas linhas
Pflu = zeros(NLIN, 2); % Fluxos ativos nas linhas
Qflu = zeros(NLIN, 2); % Fluxos reativos nas linhas

% Perdas nas linhas
Pperdas = zeros(NLIN, 1); % Perdas ativas nas linhas
Qperdas = zeros(NLIN, 1); % Perdas reativas nas linhas

% Cálculo dos fluxos e perdas
for ilin=1:NLIN
    kbar = SB(ilin); % Barra 'de' da linha ilin 
    mbar = EB(ilin); % Barra 'para' da linha ilin
    
    Gkm = G(kbar, mbar);
    Bkm = B(kbar, mbar);
    Tkm = T(kbar)-T(mbar);
    Gmk = G(mbar, kbar);
    Bmk = B(mbar, kbar);
    Bsh = Xsh_E(ilin);
    
    Pflu(ilin, 1) = -Gkm*V(kbar)^2+V(kbar)*V(mbar)*(Gkm*cos(Tkm)+Bkm*sin(Tkm));
    Qflu(ilin, 1) = (Bkm-Bsh)*V(kbar)^2+V(kbar)*V(mbar)*(Gkm*sin(Tkm)-Bkm*cos(Tkm));
    Pflu(ilin, 2) = -Gmk*V(mbar)^2+V(kbar)*V(mbar)*(Gmk*cos(Tkm)-Bmk*sin(Tkm));
    Qflu(ilin, 2) = (Bkm-Bsh)*V(mbar)^2-V(kbar)*V(mbar)*(Gmk*sin(Tkm)+Bmk*cos(Tkm));
    
    
    %Cálculo das perdas
    Pperdas(ilin,1) = Pflu(ilin, 1) + Pflu(ilin, 2);
    Qperdas(ilin,1) = Qflu(ilin, 1) + Qflu(ilin, 2);
end

% Verifica barras PV que atingiram limite de reativo
QLIM__ = num2str(zeros(NBAR, 1),'%06d');
for kbar = 1:NBAR
    if tipB(kbar) == 0 && tipB0(kbar) == 1
        QLIM__(kbar,:) = 's     ';
    else
        QLIM__(kbar,:) = '      ';
    end
end
%
%
end

