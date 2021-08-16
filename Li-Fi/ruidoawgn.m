
function rx_ruido=ruidoawgn(snr,pot_rx,Smod_tx_DC,M)

        if M == 4
            Es = 0.2345/256;
        elseif M==16
            Es = 1.1724/256;
        elseif M == 64
            Es = 4.9240/256;
        end
        
        ebno_v=10^(snr/10);
        sigma=sqrt(Es/(2*log2(M)*ebno_v));
        ruido=sigma*randn(1,length(Smod_tx_DC));% + i*sigma*randn(1,length(Smod_tx_canal));
        rx_ruido = pot_rx+ruido;

end