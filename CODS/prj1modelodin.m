function sis = prj1modelodin(t,y,flag,par)

    % ===================================================
    % PROJETO – IMPLEMENTAÇÃO DE UM SISTEMA DE CONTROLE 
    % EM UM REATOR CSTR COM RESPOSTA INVERSA. 
    % ===================================================
    
    % MODELO DINAMICO DESENVOLVIDO A PARTIR
    % DO MODELO EE JA IMPLEMENTADO.

    QH  = par(1);
    CA0 = par(2);
    T0  = par(3);
    u   = par(4);
    
   sis = prj1modelo(y,QH,CA0,T0,u);
    
end