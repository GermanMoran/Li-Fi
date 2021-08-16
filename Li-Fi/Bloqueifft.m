%Este bloque permite calcular la IFFT

function ofdm_hermi=Bloqueifft(y)

    %Conformo el complemento de los simbolos para que cumplan la condicion de simetria hermitica
    sim_conj=[];

    tam_sim=size(y);                            %Variable para guardar el numero de columnas necesarias para invertir 
    
    for colum=1:tam_sim(2)
        sim_conj=conj(y(:,colum));              %Calculo el conjugado a cada simbolo
        sim_comp(:,colum)=flipud(sim_conj);     %Invierto columna para conformar simbolos hermiticos   
    end
  
    %Comformacion de la estructura ACO-OFDM
    sim_hermi=[y;sim_comp];
    tam_inf=size(sim_hermi);  
    saco=zeros(1,tam_inf(2));
       
    for columnas=1:tam_inf(2)
        for filas=1:tam_inf(1)
            saco(2*filas,columnas) = sim_hermi(filas,columnas); 
        end
    end
    
    % Calculo la IFFT a cada simbolo OFDM presente en las columnas del
    % sim_hermi y los guardo en el nuevo vector ofdm_hermi

    ofdm_hermi=ifft(saco);
end