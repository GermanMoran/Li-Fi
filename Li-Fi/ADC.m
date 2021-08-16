
function signal_rx_discreta=ADC(rx_ruido,Rate)

    
    %Parametros del filtro
    
    rollof=0.3;     %Factor de roll-off entre 0 y 1 -> 1 mas angosto y pronunciado
    N_T = 4;        %Spam del filtro
    rca=rcosfir(rollof,N_T,Rate,1,'sqrt');
    
    
    %% Recepto OFDM 
 
    signal_rx_real = rx_ruido;
    %signal_rx_imag = imag(rx_ruido);
    
    % filtrado
    
    signal_filter_real=2.*filter(rca, 1, signal_rx_real);
    %signal_filter_imag=filter(rca, 1, signal_rx_imag);
    
    % dowsample
    
    signal_down_real = downsample(signal_filter_real,Rate);
    %signal_down_imag = downsample(signal_filter_imag,Rate);
    
    signal_rx = signal_down_real(9:length(signal_down_real));% + i*signal_down_imag(9:length(signal_down_real));
    
    % Señal total recibida.
    
    signal_rx_discreta = signal_rx;
    
end


