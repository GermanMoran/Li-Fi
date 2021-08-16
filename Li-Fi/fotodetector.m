%Funcion para la potencia optica en el Fotodetector

function[xt pot_rx]=fotodetector(t_vector,h_vector,val_pot,Rate,M,tasa)

    if tasa==1
        rb=0.1;
    elseif tasa==2
        rb=0.5;
    elseif tasa==3
        rb=1;
    elseif tasa==4
        rb=2;
    end

    tipo= h_vector;
    x_rayos=find(tipo);
    num_rayos = 10;
    
    alpha=tipo(x_rayos(1:10));
    tiempos = t_vector(x_rayos(1:10));
    
    T=1;
    Tb=1/rb;
    %T_sym = 2;          % Periodo se simbolo ns
    T_sym = log2(M)*Tb;           % periodo de simbolo en ns

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

    xt=x0+x1+x2+x3+x4+x5+x6; %+x7+x8+x9;

    %responsividad del fotodetector
    
    R = 1;
    I_rx = R*(xt);  
    pot_rx = 1*I_rx;
    
end