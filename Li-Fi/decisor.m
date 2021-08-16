
function sim_rx_par_des=decisor(sim_rx_par,M)

        sim_rx_par_des=[];
        
        sim_rx_real = real(sim_rx_par);
        sim_rx_imag = imag(sim_rx_par); 
        
        signo_real = sign(sim_rx_real);
        signo_imag = sign(sim_rx_imag);
        
        norma_sim_real = abs(sim_rx_real);
        norma_sim_imag = abs(sim_rx_imag);

%% Descisor 4-QAM
if M == 4
        for ini=1:length(sim_rx_par)
           if real(sim_rx_par(ini)) >= 0 && imag(sim_rx_par(ini)) >= 0
               sim_rx_par_des(ini) = 1+1i;
           elseif real(sim_rx_par(ini)) >= 0 && imag(sim_rx_par(ini)) < 0
               sim_rx_par_des(ini) = 1-i;
           elseif real(sim_rx_par(ini)) <= 0 && imag(sim_rx_par(ini)) >= 0
               sim_rx_par_des(ini) = -1+1i;
           elseif real(sim_rx_par(ini)) <= 0 && imag(sim_rx_par(ini)) < 0
               sim_rx_par_des(ini) = -1-1i;
           end
        end
end

% Descisor 16-QAM

    if M == 16
        
        sim_rx_par_des=[];
        
        %sim_rx_par=round(sim_rx_par,1);
        
        sim_rx_real = real(sim_rx_par);
        sim_rx_imag = imag(sim_rx_par); 
        
        signo_real = sign(sim_rx_real);
        signo_imag = sign(sim_rx_imag);
        
        norma_sim_real = abs(sim_rx_real);
        norma_sim_imag = abs(sim_rx_imag);
        
% Descisor 4-QAM

        for ini=1:length(norma_sim_real)
           if norma_sim_real(ini) >= 0 && norma_sim_real(ini) <= 0.3424*2 && norma_sim_imag(ini) >= 0 && norma_sim_imag(ini) <= 0.3424*2
               sim_rx_par_des(ini) = 1+1i;
           elseif norma_sim_real(ini) >= 0 && norma_sim_real(ini) <= 0.3424*2 && norma_sim_imag(ini) > 0.3424*2
               sim_rx_par_des(ini) = 1+3i;
           elseif norma_sim_real(ini) >= 0.3424*2 && norma_sim_imag(ini) >= 0 && norma_sim_imag(ini) <= 0.3424*2
               sim_rx_par_des(ini) = 3+1i;
           elseif norma_sim_real(ini) >= 0.3424*2 && norma_sim_imag(ini) > 0.3424*2
               sim_rx_par_des(ini) = 3+3i;
           end
        end
    sim_rx_par_des = (signo_real'.*real(sim_rx_par_des)) + i*(signo_imag'.*imag(sim_rx_par_des));
    end
    
    
    if M == 64
        
        sim_rx_par_des=[];
        
        %sim_rx_par=round(sim_rx_par,1);

        sim_rx_real = real(sim_rx_par);
        sim_rx_imag = imag(sim_rx_par); 

        signo_real = sign(sim_rx_real);
        signo_imag = sign(sim_rx_imag);

        norma_sim_real = abs(sim_rx_real);
        norma_sim_imag = abs(sim_rx_imag);
    
    % Descisor 16-QAM

        for ini=1:length(norma_sim_real)
           if norma_sim_real(ini) >= 0 && norma_sim_real(ini) <= 0.3424*2 && norma_sim_imag(ini) >= 0 && norma_sim_imag(ini) <= 0.3424*2
               sim_rx_par_des(ini) = 1+1i;
           elseif norma_sim_real(ini) >= 0 && norma_sim_real(ini) <= 0.3424*2 && norma_sim_imag(ini) > 0.3424*2 && norma_sim_imag(ini) <= 0.3424*4
              sim_rx_par_des(ini) = 1+3i;
           elseif norma_sim_real(ini) >= 0 && norma_sim_real(ini) <= 0.3424*2 && norma_sim_imag(ini) > 0.3424*4 && norma_sim_imag(ini) <= 0.3424*6
              sim_rx_par_des(ini) = 1+5i;
           elseif norma_sim_real(ini) >= 0 && norma_sim_real(ini) <= 0.3424*2 && norma_sim_imag(ini) > 0.3424*6
               sim_rx_par_des(ini) = 1+7i;
           elseif norma_sim_real(ini) > 0.3424*2 && norma_sim_real(ini) <= 0.3424*4 && norma_sim_imag(ini) >= 0 && norma_sim_imag(ini) <= 0.3424*2
               sim_rx_par_des(ini) = 3+1i;
           elseif norma_sim_real(ini) > 0.3424*2 && norma_sim_real(ini) <= 0.3424*4 && norma_sim_imag(ini) > 0.3424*2 && norma_sim_imag(ini) <= 0.3424*4
                sim_rx_par_des(ini) = 3+3i;
           elseif norma_sim_real(ini) > 0.3424*2 && norma_sim_real(ini) <= 0.3424*4 && norma_sim_imag(ini) > 0.3424*4 && norma_sim_imag(ini) <= 0.3424*6
               sim_rx_par_des(ini) = 3+5i;
           elseif norma_sim_real(ini) > 0.3424*2 && norma_sim_real(ini) <= 0.3424*4 && norma_sim_imag(ini) > 0.3424*6
                sim_rx_par_des(ini) = 3+7i;
           elseif norma_sim_real(ini) > 0.3424*4 && norma_sim_real(ini) <= 0.3424*6 && norma_sim_imag(ini) >= 0 && norma_sim_imag(ini) <= 0.3424*2
               sim_rx_par_des(ini) = 5+1i;
           elseif norma_sim_real(ini) > 0.3424*4 && norma_sim_real(ini) <= 0.3424*6 && norma_sim_imag(ini) > 0.3424*2 && norma_sim_imag(ini) <= 0.3424*4
                sim_rx_par_des(ini) = 5+3i;
           elseif norma_sim_real(ini) > 0.3424*4 && norma_sim_real(ini) <= 0.3424*6 && norma_sim_imag(ini) > 0.3424*4 && norma_sim_imag(ini) <= 0.3424*6
               sim_rx_par_des(ini) = 5+5i;
           elseif norma_sim_real(ini) > 0.3424*4 && norma_sim_real(ini) <= 0.3424*6 && norma_sim_imag(ini) > 0.3424*6
              sim_rx_par_des(ini) = 5+7i;
           elseif norma_sim_real(ini) > 0.3424*6 && norma_sim_imag(ini) >= 0 && norma_sim_imag(ini) <= 0.3424*2
              sim_rx_par_des(ini) = 7+1i;
           elseif norma_sim_real(ini) > 0.3424*6 && norma_sim_imag(ini) > 0.3424*2 && norma_sim_imag(ini) <= 0.3424*4
              sim_rx_par_des(ini) = 7+3i;
           elseif norma_sim_real(ini) > 0.3424*6 && norma_sim_imag(ini) > 0.3424*4 && norma_sim_imag(ini) <= 0.3424*6
               sim_rx_par_des(ini) = 7+5i;
           elseif norma_sim_real(ini) > 0.3424*6 && norma_sim_imag(ini) > 0.3424*6
               sim_rx_par_des(ini) = 7+7i;           
           end
        end
    sim_rx_par_des = (signo_real'.*real(sim_rx_par_des)) + i*(signo_imag'.*imag(sim_rx_par_des));
    end
end

