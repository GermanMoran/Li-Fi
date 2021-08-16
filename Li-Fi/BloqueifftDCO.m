function ofdm_hermi=BloqueifftDCO(y)

    %Conformo el complemento de los simbolos para que cumplan la condicion de simetria hermitica
    sim_conj=[];
    sim_comp=[];
    tam_sim=size(y);                            %Variable para guardar el numero de columnas necesarias para invertir 
    for colum=1:tam_sim(2)
        sim_conj=conj(y(:,colum));              %Calculo el conjugado a cada simbolo
        sim_comp(:,colum)=flipud(sim_conj);     %Invierto columna para conformar simbolos hermiticos   
    end
  

    
    colum_zeros=size(y);
    sdco=[zeros(1,colum_zeros(2));y; zeros(1,colum_zeros(2)); sim_comp];
   
    % Calculo la IFFT a cada simbolo OFDM presente en las columnas del
    % sim_hermi y los guardo en el nuevo vector ofdm_hermi
    
    ofdm_hermi=ifft(sdco);
    
      
end