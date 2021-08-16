     % Codigo realizado en matlab con la simulación de un canal
     % multitrayecto y con ruido, sin considerar el prefijo ciclico.
    clear all
    close all
    clc
    
    tic                     % Comando utilizado para medir el tiempo de simulación
    N_inf=64;               % Numero de portadoras de informacion
    DP=1/4;                 % Duración del prefijo ciclico.
    num_cp=round(DP*N_inf); % Número de portadora de CP
    Vco=6;                  % Valor de voltaje de conducción del LED
    bin=[];                 % Vector para guardar los datos binarios
    data=randi([0,1], 1, 10000);    % Generamos datos aleatorios
    num_datos = length(data);       % Longitud de datos creados
    pos = 1;                % Variable para guardar valores de BER
    pasos_ebno=0.5;         % Variaciones de Eb/No
    ebno_max=26;            % Maxima ebno alcanzable
    Ber=[];                 % Vector para guardar los diferentes valores de BER
    
    aux1=0;
    aux2=1;                 % Auxiliares utilizados para guardar diferentes 
    aux3=0;                 % valores de BER de diferentes modulaciones
    aux=0;

    bin = data;             % Cargo data en el vector bin 
    
    %% Relleno de ceros para generar un flujo de bits multiplo del esquema de modulación

    for cambio=2:2:6

    M = 2^cambio;       % Modulación utilizada
    S=log2(M);          % Número de bits por simbolo QAM
    Rb = 2;             % Tasa de datos en Gbps
    T_bit=1/Rb;         % Tiempo de simbolo en nano segundos
    
    Sym_rell = (N_inf*S)-mod(length(bin), N_inf*S);   % Calculo el número de bits faltantes
    
    ceros_fluj=randi([0 1], 1, Sym_rell);             % Creo el flujo de ceros faltantes
    
    bin_tx = [bin ceros_fluj];                  % Agrego los ceros faltantes al bloque de información
    
%% Divido los simbolos en grupos de 64, para asignarlos cada uno a una subportadora

    bin_tx=reshape(bin_tx, N_inf*S, length(bin_tx)/(N_inf*S)); 
   
    y = qammod(bin_tx, M, 'InputType', 'bit');  % Modulacion QAM    
    tam_y=size(y);
    per  = reshape(y,1,tam_y(1)*tam_y(2));
    norm = abs(real(max(per)));
    per=per/sqrt(2);
    y = reshape(per,tam_y(1),tam_y(2));
    % Determinar potencia pico y promedio del orden de la simulación.
    meanPower = mean (abs(y).^2);
    peakPower = max (abs(y).^2);
    %constelacion qam
    scatterplot(reshape(y,1,tam_y(1)*tam_y(2)));
    % Conformo el complemento de los simbolos para que cumplan la condicion de simetria hermitica
    sim_conj=[];

    tam_sim=size(y);                  % Variable para guardar el numero de columnas necesarias para invertir 
    
    sim_conj=[];
    sim_comp=[];
    
    for colum=1:tam_sim(2)
        sim_conj=conj(y(:,colum));    % Calculo el conjugado a cada simbolo
        sim_comp(:,colum)=flipud(sim_conj); % Invierto columna para conformar simbolos hermiticos   
    end

    % Conformo los simbolos hermiticos y calculo la IFFT a cada uno de ellos
    %Estructura de la señal ACO-OFDM
    
    sim_hermi=[y;sim_comp];
    
    tam_inf=size(sim_hermi);  % Variable para calcular el tamaño de la informacion en filas y columnsa luego del bloque de simetria hermitica
    
    saco=zeros(1,tam_inf(2));
       
    for columnas=1:tam_inf(2)
        for filas=1:tam_inf(1)
            saco(2*filas,columnas) = sim_hermi(filas,columnas); 
        end
    end
    
    % Calculo la IFFT a cada simbolo OFDM presente en las columnas del
    % sim_hermi y los guardo en el nuevo vector ofdm_hermi

    %  Recorte a cero a cada simbolo OFDM


    %% Conversion de paraleo a serie de la matriz de información

    sim_ofdm=ifft(saco);
    
    sim_CP=sim_ofdm;

     % Recorte a cero
      tam_ofdm=size(saco);      % Tamaño total de un simbolo OFDM
% % 
    for columnas_sim=1:tam_ofdm(2)   % Recorro cada columna de la matriz de información
        for filas_sim=1:tam_ofdm(1)  % Recorro cada fila de la matriz de información
            if sim_CP(filas_sim,columnas_sim) < 0   % ubico las portadoras negativas
               sim_CP(filas_sim,columnas_sim) = 0;  % hago cero las portadoras negativas
            end
        end
    end
    
    
    sim_serie=reshape(sim_CP,[1,tam_ofdm(1)*tam_ofdm(2)]);
      
    
    %% Conversion digital analogica de la señal de información
    
    rollof=0.5; % Factor de roll-off entre 0 y 1 -> 1 mas angosto y pronunciado
    N_T = 4;    % spam del filtro
    Rate=16*N_inf;   % Factor de sobre muestreo-> numero de muestras por cada pulso 

    Amp=[sim_serie zeros(1,8)];   % Amplitud, relleno al final con ceros.
    rca=rcosfir(rollof,N_T,Rate,1,'sqrt');
    signal_filter=filter(rca,1,upsample(Amp,Rate)); %los pulsos ya estan conformados

    Smod_tx_canal = signal_filter;
    
    Smod_tx_DC = signal_filter + Vco;
    
   
%% Calculo de polinomios para realizar la aproximación de V---I y Corriente --- potencia optica para modelar el LED

% 1. calculamos el vector de voltaje y corriente, de acuerdo a la grafica.
vol = [5.375 5.600 5.800 5.960 6.125 6.300 6.450 6.625 6.750 6.900 7.050 7.175 7.325 7.450 7.600 7.750];
I = [0.025 0.050 0.075 0.100 0.125 0.150 0.175 0.200 0.225 0.250 0.275 0.300 0.325 0.350 0.375 0.400];

% Calculo de la función Voltaje corriente.
coef = polyfit(vol,I,5);        % Coeficientes del polinomio de la señal voltaje corriente.

% Calculo de la función lineal voltaje-corriente.

voll = [6.8 5.8]; 
Il = [0.235 0.075]; 
coef_lin = polyfit(voll,Il,1); 
xvil = linspace(5.8, 6.8, 100);

val_vil=polyval(coef_lin, xvil);      % Valores voltaje-corriente en el rango

%Calculo del polinomio para la función corriente--voltaje 

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
coef_ip = polyfit(I2,poto,6);  % Coeficientes corriente potencia optica

% Paso los valores de voltaje de la señal de predistorsión por la grafica
% V_I

val_pred=polyval(coef_pred, Smod_tx_DC);    % Evaluo la señal de voltaje de predistorsión en la función voltaje-corriente
val_I=polyval(coef, val_pred);    % Evaluo la señal de voltaje de predistorsión en la función voltaje-corriente
val_F=polyval(coef_if, val_I);    % Evaluo la señal de corriente en la función de corriente-flujo luminoso
val_pot=polyval(coef_ip, val_I);  % Evaluo la señal de corriente en la funcion de corriente-potencia optica.
%val_pot=signal_filter;

%% Respuesta al impulso en diferentes puntos de rx 

%Programa para calcular la multitrayectoria(Caso ideal de 4 paredes)
%Tiempo en ns (Depuracion del programa
C=3e8*1e-9;
%Parametros iniciales (angulo media potencia, indice de refleccion paredes)
theta=60;
m=-log10(2)/log10(cosd(theta));
Aef=7.02e-3;
p=0.5;
Ts=1;
%index=1.5;
FOV=60;
%G_Con=(index^2)/(sind(FOV).^2);

hpr = 0.85;

%Dimensiones del cuarto
ancho=5; largo=5; alto=3;

h_tx_rx=alto-hpr;

Nx=ancho*10; Ny=largo*10; Nz=round(50);
%dA=17.4e-3;
dA=h_tx_rx*largo/(Ny*Nz);

x=0:ancho/Nx:ancho;
y=0:largo/Ny:largo;
z=hpr:alto/Nz:alto;

[X,Y,Z]=meshgrid(x,y,z);

%Posición del transmisor
rt=[2.5 2.5 3];
%Vector normal del transmisor
nt=[0 0 -1];
%Vector normal del receptor
nr=[0 0 1];

%Vectores normales de c/d pared
VN1=[1 0 0];
VN2=[-1 0 0];
VN3=[0 -1 0];
VN4=[0 1 0];

%Tiempo de resolucion
%delta_t=1/44.36;

delta_t=1/8;

%% codigo principal respuesta al impulso.

     r = [4.9 2.5 hpr];
     %r = [2.5 2.5 hpr];
     %r = [0.1 0.1 hpr];

     t_vector=0:25/(delta_t); % vector de tiempo en ns
     h_vector=zeros(1,length(t_vector));
     % receiver position vector
     %LOS channel gain
     d=sqrt(dot(rt-r,rt-r));
     cos_phi=dot(nt,(r-rt))/d;
     tau0=d/C;
     index=find(round(tau0/delta_t)==t_vector);
         if abs(acosd(cos_phi))<=FOV
         h_vector(index)=Aef*(m+1).*(cos_phi.^m).*(cos_phi)./(2*pi.*d.^2);
         end
         
         
        %Reflexión de la primera pared
        count=1;

            for k=1:length(y)
                for l=1:length(z)
                wr=[0,y(k),z(l)] ;           %Posicion del diferencial de area pared1
                d1=sqrt(dot(rt-wr, rt-wr));
                cos_phi=dot(nt,wr-rt)/d1;
                cos_alpha=dot(VN1,rt-wr)/d1;
                d2=sqrt(dot(wr-r,wr-r));
                cos_beta=dot(VN1,r-wr)/d2;
                cos_tetha1=dot(nr,wr-r)/d2; 
                
                tau1=(d1+d2)/C;
                index=find(round(tau1/delta_t)==t_vector);
                     if abs(acosd(cos_tetha1))<=FOV
                     h_vector(index)=h_vector(index)+(m+1)*p*Aef*dA*(cos_phi^m)*(cos_alpha)*(cos_beta)*(cos_tetha1)/(2*pi^2*d1^2*d2^2);
                     end
                count=count+1;
                 end
            end
            %Reflexion de la segunda Pared
            count=1;
            for k=1:length(y)
                 for l=1:length(z)

                  wr2=[ancho, y(k),z(l)];
                  d1=sqrt(dot(rt-wr2,rt-wr2));
                  cos_phi=dot(nt,wr2-rt)/d1;
                  cos_alpha=dot(VN2,rt-wr2)/d1;
                  d2=sqrt(dot(wr2-r,wr2-r));
                  cos_beta=dot(VN2,r-wr2)/d2;
                  cos_tetha1=dot(nr,wr2-r)/d2;
                  
                 
                 tau2=(d1+d2)/C;
                 index=find(round(tau2/delta_t)==t_vector);
                     if abs(acosd(cos_tetha1))<=FOV
                     h_vector(index)=h_vector(index)+(m+1)*p*Aef*dA*(cos_phi^m)*(cos_alpha)*(cos_beta)*(cos_tetha1)/(2*pi^2*d1^2*d2^2);
                     end
                 count=count+1;
                 end
            end
            %Reflexion de la tercera pared
            count=1;
            for k=1:length(x)
                 for l=1:length(z)

                  wr3=[x(k), largo, z(l)];
                  d1=sqrt(dot(rt-wr3,rt-wr3));
                  cos_phi=dot(nt,wr3-rt)/d1;
                  cos_alpha=dot(VN3,rt-wr3)/d1;
                  d2=sqrt(dot(wr3-r,wr3-r));
                  cos_beta=dot(VN3,r-wr3)/d2;
                  cos_tetha1=dot(nr,wr3-r)/d2;
                  
                 
                 tau3=(d1+d2)/C;
                 index=find(round(tau3/delta_t)==t_vector);
                     if abs(acosd(cos_tetha1))<=FOV
                     h_vector(index)=h_vector(index)+(m+1)*Aef*p*dA*cos_phi^m*cos_alpha*cos_beta*cos_tetha1/(2*pi^2*d1^2*d2^2);
                     end
                 count=count+1;
                 end
            end
            %Reflection de la cuarta pared
             count=1;
             for k=1:length(x)
                 for l=1:length(z)
                  wr4=[x(k), 0, z(l)];
                  d1=sqrt(dot(rt-wr4,rt-wr4));
                  cos_phi=dot(nt,wr4-rt)/d1;
                  cos_alpha=dot(VN4,rt-wr4)/d1;
                  d2=sqrt(dot(wr4-r,wr4-r));
                  cos_beta=dot(VN4,r-wr4)/d2;
                  cos_tetha1=dot(nr,wr4-r)/d2;
            
                 tau4=(d1+d2)/C;
                 index=find(round(tau4/delta_t)==t_vector);
                     if abs(acosd(cos_tetha1))<=FOV
                     h_vector(index)=h_vector(index)+(m+1)*p*Aef*dA*(cos_phi^m)*(cos_alpha)*(cos_beta)*(cos_tetha1)/(2*pi^2*d1^2*d2^2);
                     end
                     
                end
             end
             t_vector=t_vector*delta_t;
             tipo=h_vector;

             

%% Canal de comunicaciones multitrayecto.

    x_rayos=find(tipo);
    num_rayos = 10;
    T=1;
    alpha=tipo(x_rayos(1:10));
    tiempos = t_vector(x_rayos(1:10));
    
    T_sym = log2(M)*T_bit;           % periodo de simbolo en ns

    for num=2:num_rayos
        tao(:,num-1)=(tiempos(num)-tiempos(1))/T_sym;
        alpha_percent(:,num-1)=(alpha(num)/alpha(1));
    end
    
    x0 = val_pot;
    
    retardo1 = zeros(1, round(Rate*T*tao(1)));
    retardo2 = zeros(1, round(Rate*T*tao(2)));
    retardo3 = zeros(1, round(Rate*T*tao(3)));
    retardo4 = zeros(1, round(Rate*T*tao(4)));
    retardo5 = zeros(1, round(Rate*T*tao(5)));
    retardo6 = zeros(1, round(Rate*T*tao(6)));
    retardo7 = zeros(1, round(Rate*T*tao(7)));
    retardo8 = zeros(1, round(Rate*T*tao(8)));
    retardo9 = zeros(1, round(Rate*T*tao(9)));

    x1 = alpha_percent(1)*[retardo1 x0(1:length(x0)-length(retardo1))];
    x2 = alpha_percent(2)*[retardo2 x0(1:length(x0)-length(retardo2))];
    x3 = alpha_percent(3)*[retardo3 x0(1:length(x0)-length(retardo3))];
    x4 = alpha_percent(4)*[retardo4 x0(1:length(x0)-length(retardo4))];
    x5 = alpha_percent(5)*[retardo5 x0(1:length(x0)-length(retardo5))];
    x6 = alpha_percent(6)*[retardo6 x0(1:length(x0)-length(retardo6))];
    x7 = alpha_percent(7)*[retardo7 x0(1:length(x0)-length(retardo7))];
    x8 = alpha_percent(8)*[retardo8 x0(1:length(x0)-length(retardo8))];
    x9 = alpha_percent(9)*[retardo9 x0(1:length(x0)-length(retardo9))];

    xt=x0+x1+x2+x3+x4+x5+x6;%+x7+x8+x9;

    %responsividad del fotodetector
    
    R = 1;
    I_rx = R*(xt);  
    pot_rx = 1*I_rx;

    
    % Ruido AWGN: se agrega en el dominio electrico en Rx, debido a que un canal optico presenta muy bajo ruido
    %aux=0;

    %%
    for aux1=1:1:3
    aux=aux+1;
    pos=1;
    aux3=aux3+1;
    for snrdb=0:pasos_ebno:ebno_max

    % Calculo las curvas de SNR vs BER.
    
        if M == 4
            Es = 0.2345/256;
        elseif M==16
            Es = 1.1724/256;
        elseif M == 64
            Es = 4.9240/256;
        end
        

        ebno_v=10^(snrdb/10);
        
        sigma=sqrt(Es/(2*log2(M)*ebno_v));
    
        ruido=sigma*randn(1,length(Smod_tx_DC));% + i*sigma*randn(1,length(Smod_tx_canal));
     
        rx_ruido = pot_rx+ruido;
        

    % Se filtra la señal con el pulso de raiz de coseno alzado co ruido
%% Recepto OFDM 
 
    signal_rx_real = real(rx_ruido);
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
   
% Encuentro el inicio de la convolucion de la señal y sumo la duracion del
% simbolo OFDM Li-Fi para su posterior procesamiento


%% Conversión serial a paralelo y remuevo el CP.

sim_par=reshape(signal_rx_discreta, tam_ofdm(1),tam_ofdm(2));  % Conversion serial a paralelo

%sim_sin_cp=sim_par(num_cp+1:end,:);             % Elimino el prefijo ciclico

sim_sin_cp=sim_par;

tam_sin_cp = size(sim_sin_cp);

%% Convierto la señal al dominio de la frecuencia para su posterior demapeo.

sim_rx=fft(sim_sin_cp);

%% Estimación del canal (Ecualización)
%scatterplot(sim_rx);
% Respuesta en frecuencia al impulso

sim_equ=sim_rx;

%% Extraigo portadoras impares de información

    tam_sim_rx = size(sim_equ);
    sim_hermi_rx=[];
    
    for columnas=1:tam_sim_rx(2)
        for filas=1:(tam_sim_rx(1))/2
            sim_hermi_rx(filas,columnas) = sim_equ(2*filas,columnas); 
        end
    end

% extraigo los simbolos conjugados para obtener los simbolos de información
    
sim_sin_hermi=sim_hermi_rx(1:N_inf,:);

%% Conversión paralelo a serie, demodulación de simbolos, descisor

sim_rx_par=reshape(sim_sin_hermi,[tam_sim(1)*tam_sim(2),1]);

sim_rx_par=sim_rx_par;

        sim_rx_par_des=[];
        
        %sim_rx_par=round(sim_rx_par,1);
        
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
    
% %% Demodulación M-QAM

%scatterplot(sim_rx_par_des);

sim_rx_par_dem=reshape(sim_rx_par_des, tam_y(1), tam_y(2));

sim_dem = qamdemod(sim_rx_par_dem, M, 'OutputType', 'bit');      % Demodulación.

tam_sim_dem=size(sim_dem);

sim_dem = reshape(sim_dem, 1, tam_sim_dem(1)*tam_sim_dem(2));

sim_dem= sim_dem(1:length(sim_dem)-length(ceros_fluj));    % Quito el relleno de cerosm utilizadas para conformar bloques exactos


%% Conformar los pixeles en formato decimal uint 8

tam_sim_rx=size(sim_dem);            % variable para calcular el tamaño serial de la informacion rx

sim_rx_serie = reshape(sim_dem, 1, tam_sim_rx(1)*tam_sim_rx(2));          % Convierto la información rx de paralelo a serial

bits_errados=(sim_rx_serie==data);
bits_errados=sum(~bits_errados);
Ber(aux, pos)=bits_errados/length(data);
Ber_total(aux3, pos)=bits_errados/length(data);
pos=pos+1;
    end
    end
    Ber_sum(aux2,:)=mean(Ber);
    aux2=aux2+1;
    end
   
%%

%xlswrite('ACO_pared_2gbps', Ber_sum, 'Hoja1','A1');
%xlswrite('ACO_pared_2gbps', Ber_total, 'Hoja2','A1');

pos=1;
pasos_ebno=0.1;
ebno_max=10;

Ber_teo=[];
for ebno=0:pasos_ebno:ebno_max

M=4;
    
ebno_v=10^(ebno/10);  
pr_e=4*(1/log2(M))*(1-1/sqrt(M))*qfunc(sqrt((3*log2(M)*ebno_v)/(M-1)));

Ber_teo(pos)=pr_e;

pos=pos+1;

end

%Ber_prom=mean(Ber);
snrdb=0:pasos_ebno:ebno_max;
b_i=(ebno_max/pasos_ebno);
ber4=  Ber_sum(1,1:end);
ber16=  Ber_sum(2,1:end);
ber64= Ber_sum(3,1:end);
     
%semilogy(snrdb, Ber_teo,'r-*', snrdb, ber4,'k-*')

ebno_graf=0:pasos_ebno:ebno_max;

semilogy(ebno_graf, Ber_teo, 'r-*' ,snrdb, ber4, 'k-*', snrdb, ber16,'b-*', snrdb, ber64,'m-*')

title('ACO-OFDM')
grid on 
xlabel('EbNo(dB)')
ylabel('BER')
legend('4-QAM','16-QAM','64-QAM')
%legend('4-QAM teorica','4-QAM ACO-OFDM')

toc
