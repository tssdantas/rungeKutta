clc, clear

% ===================================================
% PROJETO – IMPLEMENTAÇÃO DE UM SISTEMA DE CONTROLE 
% EM UM REATOR CSTR COM RESPOSTA INVERSA. 
% ===================================================


% VARIAVEIS DE ESTADO: CA, CB, T.


% CHUTE INICIAL (BASEADO NA INF. DO TEXTO.)
% CAee=1.25  CBee=0.90 Tee=407.15
x0 = [1.,0.9,400]';


% ATRIBUICAO DE VALORES: VARIAVEIS EXTERNAS
QH  = -451.509 ;
u   = 19.52   ;  %F/V (taxa de diluicao)
CA0 = 5.0     ;
T0  = 403.15  ;


% CALCULANDO OS ZEROS COND. ESPECIFICADAS;
[xee,fval] = fsolve('prj1modelo',x0,[],QH,CA0,T0,u);

disp('Ponto de operacao: ')
xee







