
function sim_rx_par = Bloquefftdco(signal_rx_discreta,N_inf,tam_ofdm,y)

%tam_ofdm=size(tam_ofdm);
tam_sim=size(y); 
%Conversión serial a paralelo y remuevo el CP.
sim_par=reshape(signal_rx_discreta, tam_ofdm(1),tam_ofdm(2));  % Conversion serial a paralelo

%sim_sin_cp=sim_par(num_cp+1:end,:);                   %Elimino el prefijo ciclico
sim_sin_cp=sim_par;

%Convierto la señal al dominio de la frecuencia para su posterior demapeo.
sim_rx=fft(sim_sin_cp);
%scatterplot(reshape(sim_rx,[130*8,1]))
%Estimación del canal (Ecualización)
sim_equ=sim_rx;

%Elimino portadoras DC y la simetria hermitica.
sim_hermi_rx=sim_equ(2:N_inf+1,:);

%Conversión paralelo a serie, demodulación de simbolos, descisor
sim_rx_par=reshape(sim_hermi_rx,[tam_sim(1)*tam_sim(2),1]);

end
