%Función conversion Analoga/Digital

function [Smod_tx_DC,val_pot,Es]=ConversionDAC(bandera,sim_serie,N_sub,Rate)


        if bandera==1

            Vco=6;
            rollof=0.3; % Factor de roll-off entre 0 y 1 -> 1 mas angosto y pronunciado
            N_T = 4;    % spam del filtro
       

            Amp=[sim_serie zeros(1,8)];   % Amplitud, relleno al final con ceros.
            rca=rcosfir(rollof,N_T,Rate,1,'sqrt');
            signal_filter=filter(rca,1,upsample(Amp,Rate)); %los pulsos ya estan conformados

            Smod_tx_DC = signal_filter + Vco;
            
            
            % 1. calculamos el vector de voltaje y corriente, de acuerdo a la grafica.
            vol = [5.375 5.600 5.800 5.960 6.125 6.300 6.450 6.625 6.750 6.900 7.050 7.175 7.325 7.450 7.600 7.750];
            I = [0.025 0.050 0.075 0.100 0.125 0.150 0.175 0.200 0.225 0.250 0.275 0.300 0.325 0.350 0.375 0.400];

            % Calculo de la función Voltaje corriente.
            coef = polyfit(vol,I,5);        % Coeficientes del polinomio de la señal voltaje corriente.
            xvi = linspace(5.375, 7.75, 100);  % Rango funcion voltaje-corriente

            % Calculo de la función lineal voltaje-corriente.

            voll = [6.8 5.8]; 
            Il = [0.235 0.075]; 
            coef_lin = polyfit(voll,Il,1); 
            xvil = linspace(5.8, 6.8, 100);

            val_vil=polyval(coef_lin, xvil);      % Valores voltaje-corriente en el rango

            %Calculo del polinomio para la función corriente--voltaje 

            xiv = linspace(0.025, 0.4, 100);   % Rango función corriente-voltaje

            coef_iv = polyfit(I, vol, 5);   % Coeficientes del polinomio corriente-voltaje

            val_pre=polyval(coef_iv, val_vil); % Evaluo la función lineal en la función corriente-voltaje

            % Generar polinomio de predistorción.

            coef_pred = polyfit(xvil, val_pre, 4); 

            %2. calculamos el vector de corriente a flujo, de acuerdo a la grafica
            %relativo a 350mA (71lm)

            I2 = [0 25 50 75 100 125 150 175 200 225 250 275 300 325 350 375 400]*1/1000;
            flujo = [0 0.2 0.35 0.55 0.7 0.85 1 1.15 1.3 1.4 1.55 1.65 1.8 1.9 2.05 2.15 2.3];
            poto = [0 0.2 0.35 0.55 0.7 0.85 1 1.15 1.3 1.4 1.55 1.65 1.8 1.9 2.05 2.15 2.3]*(144/300); % W

            % Calculo del coeficiente de los polinomios.

            coef_if = polyfit(I2,flujo,6);  % Coeficientes corriente flujo luminoso

            x = linspace(0,0.4,100);

            fig_pol=polyval(coef_if,x);

            % 3. calculamos el vector de de corriente a potencia optica, de acuerdo a la grafica
            % relativo a 300lm/W

            % Calculo del coeficiente de los polinomios 

            coef_ip = polyfit(I2,poto,6);  % Coeficientes corriente potencia optica
            x = linspace(0,0.4,100);
            fig_pol=polyval(coef_ip,x);

            % Paso los valores de voltaje de la señal de predistorsión por la grafica
            % V_I

            val_pred=polyval(coef_pred, Smod_tx_DC);    % Evaluo la señal de voltaje de predistorsión en la función voltaje-corriente
            val_I=polyval(coef, val_pred);    % Evaluo la señal de voltaje de predistorsión en la función voltaje-corriente
            val_F=polyval(coef_if, val_I);    % Evaluo la señal de corriente en la función de corriente-flujo luminoso
            val_pot=polyval(coef_ip, val_I);  % Evaluo la señal de corriente en la funcion de corriente-potencia optica.
            
            %% Calculo de la energia de simbolo

            val_pred_ES=polyval(coef_pred, 6);    % Evaluo la señal de voltaje de predistorsión en la función voltaje-corriente
            val_I_ES=polyval(coef, val_pred_ES);    % Evaluo la señal de voltaje de predistorsión en la función voltaje-corriente
            val_pot_ES=polyval(coef_ip, val_I_ES);  % Evaluo la señal de corriente en la funcion de corriente-potencia optica.

            signal_filter_ES=filter(rca, 1, val_pot-val_pot_ES);

            signal_down_ES = downsample(signal_filter_ES,Rate);

            signal_rx_ES = signal_down_ES(9:length(signal_down_ES));% + i*signal_down_imag(9:length(signal_down_real));

            %sim_sin_cp=sim_par(num_cp+1:end,:);             % Elimino el prefijo ciclico

            Es = mean(real(signal_rx_ES).^2+imag(signal_rx_ES).^2)*2;
        end
end