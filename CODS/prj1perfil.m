clc, clear

% ===================================================
% PROJETO – IMPLEMENTAÇÃO DE UM SISTEMA DE CONTROLE 
% EM UM REATOR CSTR COM RESPOSTA INVERSA. 
% ===================================================


% VARIAVEIS DE ESTADO: CA, CB, T.


% CHUTE INICIAL (BASEADO NA INF. DO TEXTO.)
% CAee=1.25  CBee=0.90 Tee=407.15
x0 = [1.,0.5,400]';


% ATRIBUICAO DE VALORES: VARIAVEIS EXTERNAS
QH  = -451.509 ;
uop1   = 19.52   ;  %F/V (taxa de diluicao)
CA0 = 5.0     ;
T0  = 403.15  ;


% CALCULANDO OS ZEROS COND. ESPECIFICADAS;
[xop1,fval] = fsolve('prj1modelo',x0,[],QH,CA0,T0,uop1);



% CALCULANDO PERFIL DE CONCENTRACOES EM RELACAO A F/V
u = [4:.5:250];
xee = []

for i=1:length(u),
    i,
    [xcalc,fval] = fsolve('prj1modelo',x0,[],QH,CA0,T0,u(i));
    xee = [xee xcalc];
end

plot(u,xee(2,:))
hold on
plot(uop1,xop1(2),'r*')
xlabel('F/V (h^{-1})')
ylabel('CB (mol/l)')
grid()








