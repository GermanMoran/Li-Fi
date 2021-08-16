%Modelado de Ruido AWGN_DCO-OFDM
function rx_ruido=RuidoAwgnDco(snr,pot_rx,Smod_tx_DC,M,Es)
      
        ebno_v=10^(snr/10);
        sigma=sqrt(Es/(2*log2(M)*ebno_v));
        ruido=sigma*randn(1,length(Smod_tx_DC));% + i*sigma*randn(1,length(Smod_tx_canal));
        rx_ruido = pot_rx+ruido;

end