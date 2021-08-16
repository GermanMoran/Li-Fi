function [sim_serie tam_ofdm]=recortecerodco(ofdm_hermi)
 
%Adiciono la componente DC y el recorte a cero a cada simbolo OFDM


 sim_CP=ofdm_hermi;
 DC=4*std(sim_CP(:,:));  % Calculo la desviacion estandar de la señal
 sim_DC=sim_CP+DC;       % Simbolo con componente DC
 %Recorte a cero
 tam_ofdm=size(sim_DC);                             %Tamaño total de un simbolo OFDM

    for columnas_sim=1:tam_ofdm(2)                  %Recorro cada columna de la matriz de información
        for filas_sim=1:tam_ofdm(1)                 %Recorro cada fila de la matriz de información
            if sim_DC(filas_sim,columnas_sim) < 0   %Ubico las portadoras negativas
               sim_DC(filas_sim,columnas_sim) = 0;  %Hago cero las portadoras negativas
            end
        end
    end

    %Modificamos algo
    
    %Conversion de paraleo a serie de la matriz de información

    sim_serie=reshape(sim_DC,[1,tam_ofdm(1)*tam_ofdm(2)]);
    
end