clc, clear

% ====================================================================
% ====================================================================
% ATRIBUICAO DE VALORES: 
%
% VARIAVEIS EXTERNAS
QH  = -451.509 ;
u   = 19.52   ;  %F/V (taxa de diluicao)
CA0 = 5.0     ;
T0  = 403.15  ;

% ====================================================================
% ====================================================================

% CALCULO DO PONTO DE OPERACAO

% CHUTE INICIAL (BASEADO NA TAB.)
x0 = [1.,0.9,400]';

% CALCULANDO OS ZEROS PRIMEIRO PONTO DE OPERACAO
[xop1,fval] = fsolve('prj1modelo',x0,[],QH,CA0,T0,u);


% APLICANDO PERTUBACOES
%upert = u+25;
CA0pert = CA0+1;
[xop2,fval] = fsolve('prj1modelo',x0,[],QH,CA0pert,T0,u);


disp('Ponto de operacao 1: ')
xop1

disp('Ponto de operacao 2: ')
xop2


% ====================================================================
% ====================================================================
% SIMULACAO PARTINDO DO PONTO DE OPERACAO 2 PARA O PONTO DE OPERACAO 1
% (APLICANDO PERTUBACOES NAS VELOCIDADES DAS BOMBAS)
% OP2: [56,44]  >> OP1: [50,50]

par = [QH,CA0pert,T0,u];


% PASSO DE INTEGRACAO
Ts = 1e-3; % (passo pequeno pois esta em horas)

% VETOR VAR. INDEPENDENTE (TEMPO)
t = 0:Ts:0.25;

% SIMULACAO
[t,x] = ode23tb('prj1modelodin',t,xop1,[],par); 


% ====================================================================
% ====================================================================
% SAIDA

figure()

t = t*60;

subplot(3,1,1)
plot(t,x(:,1)), hold on

plot(t(1),xop1(1),'ro')
plot(t(end),xop2(1),'ro')

ylabel('CA (mol/l)'), xlabel('t (min)')

subplot(3,1,2)
plot(t,x(:,2))

hold on

plot(t(1),xop1(2),'ro')
plot(t(end),xop2(2),'ro')

ylabel('CB (mol/l)'), xlabel('t (min)')


subplot(3,1,3)
plot(t,x(:,3))

hold on

plot(t(1),  xop1(3),'ro')
plot(t(end),xop2(3),'ro')

ylabel('T (K)'), xlabel('t (min)')



