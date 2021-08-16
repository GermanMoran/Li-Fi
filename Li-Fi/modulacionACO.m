function [y,data,ceros_fluj]=modulacionACO(N_inf,M,TM)

    data=randi([0,1], 1, 1000);                             %Generamos datos aleatorios 
    bin = data;
    
    S=log2(M);                                              %Numero de bits por simbolo
    Sym_rell = (N_inf*S)-mod(length(bin), N_inf*S);         %Calculo el número de bits faltantes  
    ceros_fluj=randi([0 1], 1, Sym_rell);                   %Creo el flujo de ceros faltantes 
    bin_tx = [bin ceros_fluj];                              %Agrego los ceros faltantes al bloque de información
    tam_bin_tx=size(bin_tx);
    
  
    %Convierto los bits serie a paralelo segun M
    bin_tx=reshape(bin_tx, N_inf*S, length(bin_tx)/(N_inf*S));          %Agrupo bits grupos 64
    
    %Selecciono el esquema de modulación
    if TM==1
        y = qammod(bin_tx, M, 'InputType', 'bit');                      %Simbolos  M-QAM
    else
        y = pskmod(bin_tx, M);                                          %Simbolos  M-PSK
    end
    
    tam_y=size(y);
    per  = reshape(y,1,tam_y(1)*tam_y(2));
    norm = abs(real(max(per)));
    per=per/sqrt(2);
    y = reshape(per,tam_y(1),tam_y(2));
    
end



   
 
    

   