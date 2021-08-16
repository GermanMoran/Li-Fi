function sim_rx_serie= demodulacion(sim_rx_par_des,y,M,ceros_fluj) 
    
    tam_y=size(y);
    sim_rx_par_dem=reshape(sim_rx_par_des, tam_y(1), tam_y(2));
    sim_dem = qamdemod(sim_rx_par_dem, M, 'OutputType', 'bit');      % Demodulación.
    tam_sim_dem=size(sim_dem);
    sim_dem = reshape(sim_dem, 1, tam_sim_dem(1)*tam_sim_dem(2));
    sim_dem= sim_dem(1:length(sim_dem)-length(ceros_fluj));        %Quito el relleno de cerosm utilizadas para conformar bloques exactos


    %% Conformar los pixeles en formato decimal uint 8

    tam_sim_rx=size(sim_dem);            % variable para calcular el tamaño serial de la informacion rx
    sim_rx_serie = reshape(sim_dem, 1, tam_sim_rx(1)*tam_sim_rx(2));
   
end