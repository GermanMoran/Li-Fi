 
function sim_serie=recortecero(ofdm_hermi)
 
 sim_CP=ofdm_hermi;
 tam_ofdm=size(sim_CP);                             %Tamaño total de un simbolo OFDM

    for columnas_sim=1:tam_ofdm(2)                  %Recorro cada columna de la matriz de información
        for filas_sim=1:tam_ofdm(1)                 %Recorro cada fila de la matriz de información
            if sim_CP(filas_sim,columnas_sim) < 0   %Ubico las portadoras negativas
               sim_CP(filas_sim,columnas_sim) = 0;  %Hago cero las portadoras negativas
            end
        end
    end

    
    %Conversion de paraleo a serie de la matriz de información

    sim_serie=reshape(sim_CP,[1,tam_ofdm(1)*tam_ofdm(2)]);
    
end