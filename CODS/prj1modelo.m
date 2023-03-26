function sis = prj1modelo(y,QH,CA0,T0,u)

    % PARAMETROS FIXOS

    k10 = 1.287E12 ;
    k20 = 1.287E12 ;
    k30 = 9.043E09 ;

    Er1 = 9758.3 ;
    Er2 = 9758.3 ;
    Er3 = 8560.  ;

    DeltaH1 = 4.2    ;
    DeltaH2 = -11    ;
    DeltaH3 = -41.85 ; 

    rho =  0.9342 ;
    cp  =  3.01   ;
   
    % PROJETO � IMPLEMENTA��O DE UM SISTEMA DE CONTROLE 
    % EM UM CONJUNTO DE TANQUES INTERATIVOS 
    
    % VARIAVEIS DE ESTADO
    CA  = y(1)                    ;
    CB  = y(2)                    ;
    T   = y(3)                    ;
    
    % RELACOES ALGEBRICAS
    
    k1 = k10*exp(-Er1/T);
    k2 = k20*exp(-Er2/T);
    k3 = k30*exp(-Er3/T);  
                 
    % FUNCOES DO MODELO NO E.E. 
    y1dot = -k1*CA - k3*CA^2 + (CA0 - CA)*u                                                  ;
    y2dot = k1*CA - k2*CB - CB*u                                                             ;
    y3dot = (T0-T)*u + ( (-DeltaH1*k1*CA) + (-DeltaH2*k2*CB) + (-DeltaH3*k3*CA^2) + QH )/(rho*cp);

    % SISTEMAS
    sis = [y1dot; y2dot; y3dot];
    
 
