
function sim_rx_par= Bloquefft(signal_rx_discreta,N_inf,ofdm_hermi,y)


  tam_ofdm=size(ofdm_hermi);
  tam_sim=size(y); 
  sim_par=reshape(signal_rx_discreta, tam_ofdm(1),tam_ofdm(2)); 
  sim_sin_cp=sim_par;
  tam_sin_cp = size(sim_sin_cp);
  sim_rx=fft(sim_sin_cp);
  sim_equ=sim_rx;

    %% Extraigo portadoras impares de información

        tam_sim_rx = size(sim_equ);

        for columnas=1:tam_sim_rx(2)
            for filas=1:(tam_sim_rx(1))/2
                sim_hermi_rx(filas,columnas) = sim_equ(2*filas,columnas); 
            end
        end


    sim_sin_hermi=sim_hermi_rx(1:N_inf,:);
    sim_rx_par=reshape(sim_sin_hermi,[tam_sim(1)*tam_sim(2),1]);
end