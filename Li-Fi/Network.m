
%Codigo Red Li-Fi con 4 Fuentes emitiendo simultaneamente
clear all
close all
clc
%Dimenciones de la Habitacion 
largo=12.21; ancho=12.96; alto=3.08; 

%Parametros del Transmisor
%rt=[2.5, 2.5, 3];                   %Posicion del transmisor

phi_g =60;                         %angulo de media potencia en grados
m = -log10(2)/(log10(cosd(phi_g)));      %Grado  del lobulo lambertiano
Pled=480;                           %Potencia optica media maxima(mw)
nLED=1;
P_total=nLED*nLED*Pled;
nt=[0,0,-1];                        %Vector normal del Transmisor

%Parametros Receptor 
hpr=0.85;           %Altura del plano receptor donde se encuentran las mesas(m)
h=alto-hpr;         %Altura entre el plano transmisor y receptor
Aef=7.02E-4;         %Area efectiva el fotodetector
Ts=1;               %Ganancia del Filtro Optico
index=1.5;          %Indice de refraccion del lente del fotodetector
FOV=60;             %Angulo de máxima incidencia(grados)
nr=[0,0,1];         %Vector normal del receptor
G_Con=(index.^2)/(sind(FOV)).^2;             %Ganacia del concentrador Optico sin imagenes

%Superficies
Nx=ancho*5; Ny=largo*5; Nz=round(alto*5);   %numero de grillas en c/d superficie
dA=h*largo/(Ny*Nz);                            %Calculo diferencial de area
% x=0:ancho/Nx:ancho;
% y=0:largo/Ny:largo;
% z=hpr:alto/Nz:alto;
x=linspace(0,ancho,Nx);
y=linspace(0,largo,Ny);
z=linspace(hpr,alto,Nz);


[XT,YT,ZT]=meshgrid([ancho/4 (3.5*ancho)/4],[largo/4 (3.5*largo)/4],alto);   %Posición de los transmisores LED
[XR,YR,ZR]=meshgrid(x,y,hpr);                                     %Plano recepción

rt1=[2.5 2.2 3.08]
rt2=[2.5 10.01 3.08];
rt3=[10.46 2.2 3.08];
rt4=[10.46 10.01 3.08];
% rt1=[XT(1,1,1) YT(1,1,1) ZT(1,1,1)];                            %Posicion TX1
% rt2=[XT(1,2,1) YT(1,2,1) ZT(1,2,1)];                            %Posicion TX1
% rt3=[XT(2,1,1) YT(2,1,1) ZT(2,1,1)];                            %Posicion TX1
% rt4=[XT(2,2,1) YT(2,2,1) ZT(2,2,1)];                            %Posicion TX1

%Coeficiente de reflexión de c/d  Pared
p=0.6;                                                             %Coeficiente de reflexion de todas las paredes
p2=0.06;                                                           %Coeficiente de reflexion de todas las paredes
                                                        

VN1=[1,0,0];
VN2=[-1,0,0];
VN3=[0,1,0];                                   %Vectores Normales de cada pared
VN4=[0,-1,0];


dA1=(largo*h)/((Ny)*(Nz));                   %Diferencial de area de cada pared

%Primera Bombilla
for i=1:length(x)
    for j=1:length(y)
        %Calculo distribucion de potencia cuando existe LOS
        r=[x(i), y(j), hpr];                     %Vector posicion del receptor
        d=sqrt(dot(rt1-r,rt1-r));                  %Distancia entre el trasmisor y receptor
        cos_phi=dot(nt,(r-rt1))/d;                %Angulo de irradianza
        cos_tetha=dot(nr,(rt1-r))/d;              %Angulo de incidencia
        
        if  abs(acosd(cos_tetha)) >= 0 && abs(acosd(cos_tetha)) <= FOV
            HLos(i,j)= Aef*(m+1).*(cos_phi.^m).*(cos_tetha)./(2*pi.*d.^2);
        else
            HLos(i,j)=0;
        end
        
        
        %Calculo de distribucion de potencia NLOS
         H1(i,j)=0;
         H2(i,j)=0;
         H3(i,j)=0;                             %Acumuladores ganancia DC c/d pared
         H4(i,j)=0;
     
        %Pared 1
        for k=1:length(y)
            for l=1:length(z)
                wr=[0,y(k),z(l)] ;             %Posicion del diferencial de area pared 1
                d1=sqrt(dot(rt1-wr ,rt1-wr));    %Distancia entre LED y diferencial de area
                cos_phi=dot(nt,wr-rt1)/d1;      %Angulo entre vector normal del LED y el dA
                cos_alpha=dot(VN1,rt1-wr)/d1;   %Angulo entre vector normal del dA y el LED
                d2=sqrt(dot(wr-r,wr-r));       %Distancia entre el diferencial de area y el fotodetector
                cos_beta=dot(VN1,r-wr)/d2;     %Angulo entre vector normal del dA y el fotodetector
                cos_tetha=dot(nr,wr-r)/d2;     %Angulo normla entre vector normal FD Y y el dA
               if abs(acosd(cos_tetha)) <=FOV  
                    H1(i,j)=H1(i,j)+((m+1)*p*Aef*dA1*(cos_phi^m)*(cos_alpha)*(cos_beta)*(cos_tetha))/(2*pi^2*d1^2*d2^2);
               elseif abs(acosd(cos_tetha)) > FOV 
                      H1(i,j)=H1(i,j);
               end  
            end
         end
        
         
       % Pared 2
        
          for k=1:length(y)
            for l=1:length(z)
                  wr2=[ancho, y(k),z(l)];
                  d1=sqrt(dot(rt1-wr2,rt1-wr2));
                  cos_phi=dot(nt,wr2-rt1)/d1;
                  cos_alpha=dot(VN2,rt1-wr2)/d1;
                  d2=sqrt(dot(wr2-r,wr2-r));
                  cos_beta=dot(VN2,r-wr2)/d2;
                  cos_tetha=dot(nr,wr2-r)/d2;
                  
                if abs(acosd(cos_tetha)) <= FOV 
                   H2(i,j)=H2(i,j)+((m+1)*p2*Aef*dA1*(cos_phi^m)*(cos_alpha)*(cos_beta)*(cos_tetha))/(2*pi^2*d1^2*d2^2);
                 elseif abs(acosd(cos_tetha)) > FOV 
                      H2(i,j)=H2(i,j);
                end
            end
           end
        
 %        Pared 3
        
        for k=1:length(x)
            for l=1:length(z)
                  wr3=[x(k),0, z(l)];
                  d1=sqrt(dot(rt1-wr3,rt1-wr3));
                  cos_phi=dot(nt,wr3-rt1)/d1;
                  cos_alpha=dot(VN3,rt1-wr3)/d1;
                  d2=sqrt(dot(wr3-r,wr3-r));
                  cos_beta=dot(VN3,r-wr3)/d2;
                  cos_tetha=dot(nr,wr3-r)/d2;
                 
                 if abs(acosd(cos_tetha)) <= FOV 
                    H3(i,j)=H3(i,j)+((m+1)*p*Aef*dA1*(cos_phi^m)*(cos_alpha)*(cos_beta)*(cos_tetha))/(2*pi^2*d1^2*d2^2);
                 elseif abs(acosd(cos_tetha)) > FOV 
                      H3(i,j)=H3(i,j);
                  end
            end
        end
        
        %Pared 4
        
        for k=1:length(x)
            for l=1:length(z)
                  wr4=[x(k),largo, z(l)];
                  d1=sqrt(dot(rt1-wr4,rt1-wr4));
                  cos_phi=dot(nt,wr4-rt1)/d1;
                  cos_alpha=dot(VN4,rt1-wr4)/d1;
                  d2=sqrt(dot(wr4-r,wr4-r));
                  cos_beta=dot(VN4,r-wr4)/d2;
                  cos_tetha=dot(nr,wr4-r)/d2;
                 
                  if abs(acosd(cos_tetha)) <= FOV 
                    H4(i,j)=H4(i,j)+((m+1)*p*Aef*dA1*(cos_phi^m)*(cos_alpha)*(cos_beta)*(cos_tetha))/(2*pi^2*d1^2*d2^2);
                  
                  elseif abs(acosd(cos_tetha)) > FOV 
                      H4(i,j)=H4(i,j);
                   end
                  
               

            end
        end
        
    end
end

H=H1+H2+H3;
P_rec_A1=(HLos+H)*P_total.*Ts.*G_Con;


%Segunda Bombilla
for i=1:length(x)
    for j=1:length(y)
        %Calculo distribucion de potencia cuando existe LOS
        r=[x(i), y(j), hpr];                     %Vector posicion del receptor
        d=sqrt(dot(rt2-r,rt2-r));                  %Distancia entre el trasmisor y receptor
        cos_phi=dot(nt,(r-rt2))/d;                %Angulo de irradianza
        cos_tetha=dot(nr,(rt2-r))/d;              %Angulo de incidencia
        
        if  abs(acosd(cos_tetha)) >= 0 && abs(acosd(cos_tetha)) <= FOV
            HLos(i,j)= Aef*(m+1).*(cos_phi.^m).*(cos_tetha)./(2*pi.*d.^2);
        else
            HLos(i,j)=0;
        end
        
        
        %Calculo de distribucion de potencia NLOS
         H1(i,j)=0;
         H2(i,j)=0;
         H3(i,j)=0;                             %Acumuladores ganancia DC c/d pared
         H4(i,j)=0;
     
        %Pared 1
        for k=1:length(y)
            for l=1:length(z)
                wr=[0,y(k),z(l)] ;             %Posicion del diferencial de area pared 1
                d1=sqrt(dot(rt2-wr ,rt2-wr));    %Distancia entre LED y diferencial de area
                cos_phi=dot(nt,wr-rt2)/d1;      %Angulo entre vector normal del LED y el dA
                cos_alpha=dot(VN1,rt2-wr)/d1;   %Angulo entre vector normal del dA y el LED
                d2=sqrt(dot(wr-r,wr-r));       %Distancia entre el diferencial de area y el fotodetector
                cos_beta=dot(VN1,r-wr)/d2;     %Angulo entre vector normal del dA y el fotodetector
                cos_tetha=dot(nr,wr-r)/d2;     %Angulo normla entre vector normal FD Y y el dA
               if abs(acosd(cos_tetha)) <=FOV  
                    H1(i,j)=H1(i,j)+((m+1)*p*Aef*dA1*(cos_phi^m)*(cos_alpha)*(cos_beta)*(cos_tetha))/(2*pi^2*d1^2*d2^2);
               elseif abs(acosd(cos_tetha)) > FOV 
                      H1(i,j)=H1(i,j);
               end  
            end
         end
        
         
       % Pared 2
        
          for k=1:length(y)
            for l=1:length(z)
                  wr2=[ancho, y(k),z(l)];
                  d1=sqrt(dot(rt2-wr2,rt2-wr2));
                  cos_phi=dot(nt,wr2-rt2)/d1;
                  cos_alpha=dot(VN2,rt2-wr2)/d1;
                  d2=sqrt(dot(wr2-r,wr2-r));
                  cos_beta=dot(VN2,r-wr2)/d2;
                  cos_tetha=dot(nr,wr2-r)/d2;
                  
                if abs(acosd(cos_tetha)) <= FOV 
                   H2(i,j)=H2(i,j)+((m+1)*p2*Aef*dA1*(cos_phi^m)*(cos_alpha)*(cos_beta)*(cos_tetha))/(2*pi^2*d1^2*d2^2);
                 elseif abs(acosd(cos_tetha)) > FOV 
                      H2(i,j)=H2(i,j);
                end
            end
           end
        
 %        Pared 3
        
        for k=1:length(x)
            for l=1:length(z)
                  wr3=[x(k),0, z(l)];
                  d1=sqrt(dot(rt2-wr3,rt2-wr3));
                  cos_phi=dot(nt,wr3-rt2)/d1;
                  cos_alpha=dot(VN3,rt2-wr3)/d1;
                  d2=sqrt(dot(wr3-r,wr3-r));
                  cos_beta=dot(VN3,r-wr3)/d2;
                  cos_tetha=dot(nr,wr3-r)/d2;
                 
                 if abs(acosd(cos_tetha)) <= FOV 
                    H3(i,j)=H3(i,j)+((m+1)*p*Aef*dA1*(cos_phi^m)*(cos_alpha)*(cos_beta)*(cos_tetha))/(2*pi^2*d1^2*d2^2);
                 elseif abs(acosd(cos_tetha)) > FOV 
                      H3(i,j)=H3(i,j);
                  end
            end
        end
        
        %Pared 4
        
        for k=1:length(x)
            for l=1:length(z)
                  wr4=[x(k),largo, z(l)];
                  d1=sqrt(dot(rt2-wr4,rt2-wr4));
                  cos_phi=dot(nt,wr4-rt2)/d1;
                  cos_alpha=dot(VN4,rt2-wr4)/d1;
                  d2=sqrt(dot(wr4-r,wr4-r));
                  cos_beta=dot(VN4,r-wr4)/d2;
                  cos_tetha=dot(nr,wr4-r)/d2;
                 
                  if abs(acosd(cos_tetha)) <= FOV 
                    H4(i,j)=H4(i,j)+((m+1)*p*Aef*dA1*(cos_phi^m)*(cos_alpha)*(cos_beta)*(cos_tetha))/(2*pi^2*d1^2*d2^2);
                  
                  elseif abs(acosd(cos_tetha)) > FOV 
                      H4(i,j)=H4(i,j);
                   end
                  
               

            end
        end
        
    end
end

H=H1+H2+H3;
P_rec_A2=(HLos+H)*P_total.*Ts.*G_Con;


%Tercera Bombilla
for i=1:length(x)
    for j=1:length(y)
        %Calculo distribucion de potencia cuando existe LOS
        r=[x(i), y(j), hpr];                     %Vector posicion del receptor
        d=sqrt(dot(rt3-r,rt3-r));                  %Distancia entre el trasmisor y receptor
        cos_phi=dot(nt,(r-rt3))/d;                %Angulo de irradianza
        cos_tetha=dot(nr,(rt3-r))/d;              %Angulo de incidencia
        
        if  abs(acosd(cos_tetha)) >= 0 && abs(acosd(cos_tetha)) <= FOV
            HLos(i,j)= Aef*(m+1).*(cos_phi.^m).*(cos_tetha)./(2*pi.*d.^2);
        else
            HLos(i,j)=0;
        end
        
        
        %Calculo de distribucion de potencia NLOS
         H1(i,j)=0;
         H2(i,j)=0;
         H3(i,j)=0;                             %Acumuladores ganancia DC c/d pared
         H4(i,j)=0;
     
        %Pared 1
        for k=1:length(y)
            for l=1:length(z)
                wr=[0,y(k),z(l)] ;             %Posicion del diferencial de area pared 1
                d1=sqrt(dot(rt3-wr ,rt3-wr));    %Distancia entre LED y diferencial de area
                cos_phi=dot(nt,wr-rt3)/d1;      %Angulo entre vector normal del LED y el dA
                cos_alpha=dot(VN1,rt3-wr)/d1;   %Angulo entre vector normal del dA y el LED
                d2=sqrt(dot(wr-r,wr-r));       %Distancia entre el diferencial de area y el fotodetector
                cos_beta=dot(VN1,r-wr)/d2;     %Angulo entre vector normal del dA y el fotodetector
                cos_tetha=dot(nr,wr-r)/d2;     %Angulo normla entre vector normal FD Y y el dA
               if abs(acosd(cos_tetha)) <=FOV  
                    H1(i,j)=H1(i,j)+((m+1)*p*Aef*dA1*(cos_phi^m)*(cos_alpha)*(cos_beta)*(cos_tetha))/(2*pi^2*d1^2*d2^2);
               elseif abs(acosd(cos_tetha)) > FOV 
                      H1(i,j)=H1(i,j);
               end  
            end
         end
        
         
       % Pared 2
        
          for k=1:length(y)
            for l=1:length(z)
                  wr2=[ancho, y(k),z(l)];
                  d1=sqrt(dot(rt3-wr2,rt3-wr2));
                  cos_phi=dot(nt,wr2-rt3)/d1;
                  cos_alpha=dot(VN2,rt3-wr2)/d1;
                  d2=sqrt(dot(wr2-r,wr2-r));
                  cos_beta=dot(VN2,r-wr2)/d2;
                  cos_tetha=dot(nr,wr2-r)/d2;
                  
                if abs(acosd(cos_tetha)) <= FOV 
                   H2(i,j)=H2(i,j)+((m+1)*p2*Aef*dA1*(cos_phi^m)*(cos_alpha)*(cos_beta)*(cos_tetha))/(2*pi^2*d1^2*d2^2);
                 elseif abs(acosd(cos_tetha)) > FOV 
                      H2(i,j)=H2(i,j);
                end
            end
           end
        
 %        Pared 3
        
        for k=1:length(x)
            for l=1:length(z)
                  wr3=[x(k),0, z(l)];
                  d1=sqrt(dot(rt3-wr3,rt3-wr3));
                  cos_phi=dot(nt,wr3-rt3)/d1;
                  cos_alpha=dot(VN3,rt3-wr3)/d1;
                  d2=sqrt(dot(wr3-r,wr3-r));
                  cos_beta=dot(VN3,r-wr3)/d2;
                  cos_tetha=dot(nr,wr3-r)/d2;
                 
                 if abs(acosd(cos_tetha)) <= FOV 
                    H3(i,j)=H3(i,j)+((m+1)*p*Aef*dA1*(cos_phi^m)*(cos_alpha)*(cos_beta)*(cos_tetha))/(2*pi^2*d1^2*d2^2);
                 elseif abs(acosd(cos_tetha)) > FOV 
                      H3(i,j)=H3(i,j);
                  end
            end
        end
        
        %Pared 4
        
        for k=1:length(x)
            for l=1:length(z)
                  wr4=[x(k),largo, z(l)];
                  d1=sqrt(dot(rt3-wr4,rt3-wr4));
                  cos_phi=dot(nt,wr4-rt3)/d1;
                  cos_alpha=dot(VN4,rt3-wr4)/d1;
                  d2=sqrt(dot(wr4-r,wr4-r));
                  cos_beta=dot(VN4,r-wr4)/d2;
                  cos_tetha=dot(nr,wr4-r)/d2;
                 
                  if abs(acosd(cos_tetha)) <= FOV 
                    H4(i,j)=H4(i,j)+((m+1)*p*Aef*dA1*(cos_phi^m)*(cos_alpha)*(cos_beta)*(cos_tetha))/(2*pi^2*d1^2*d2^2);
                  
                  elseif abs(acosd(cos_tetha)) > FOV 
                      H4(i,j)=H4(i,j);
                   end
                  
               

            end
        end
        
    end
end

H=H1+H2+H3;
P_rec_A3=(HLos+H)*P_total.*Ts.*G_Con;



%Cuarta Bombilla
for i=1:length(x)
    for j=1:length(y)
        %Calculo distribucion de potencia cuando existe LOS
        r=[x(i), y(j), hpr];                     %Vector posicion del receptor
        d=sqrt(dot(rt4-r,rt4-r));                  %Distancia entre el trasmisor y receptor
        cos_phi=dot(nt,(r-rt4))/d;                %Angulo de irradianza
        cos_tetha=dot(nr,(rt4-r))/d;              %Angulo de incidencia
        
        if  abs(acosd(cos_tetha)) >= 0 && abs(acosd(cos_tetha)) <= FOV
            HLos(i,j)= Aef*(m+1).*(cos_phi.^m).*(cos_tetha)./(2*pi.*d.^2);
        else
            HLos(i,j)=0;
        end
        
        
        %Calculo de distribucion de potencia NLOS
         H1(i,j)=0;
         H2(i,j)=0;
         H3(i,j)=0;                             %Acumuladores ganancia DC c/d pared
         H4(i,j)=0;
     
        %Pared 1
        for k=1:length(y)
            for l=1:length(z)
                wr=[0,y(k),z(l)] ;             %Posicion del diferencial de area pared 1
                d1=sqrt(dot(rt4-wr ,rt4-wr));    %Distancia entre LED y diferencial de area
                cos_phi=dot(nt,wr-rt4)/d1;      %Angulo entre vector normal del LED y el dA
                cos_alpha=dot(VN1,rt4-wr)/d1;   %Angulo entre vector normal del dA y el LED
                d2=sqrt(dot(wr-r,wr-r));       %Distancia entre el diferencial de area y el fotodetector
                cos_beta=dot(VN1,r-wr)/d2;     %Angulo entre vector normal del dA y el fotodetector
                cos_tetha=dot(nr,wr-r)/d2;     %Angulo normla entre vector normal FD Y y el dA
               if abs(acosd(cos_tetha)) <=FOV  
                    H1(i,j)=H1(i,j)+((m+1)*p*Aef*dA1*(cos_phi^m)*(cos_alpha)*(cos_beta)*(cos_tetha))/(2*pi^2*d1^2*d2^2);
               elseif abs(acosd(cos_tetha)) > FOV 
                      H1(i,j)=H1(i,j);
               end  
            end
         end
        
         
       % Pared 2
        
          for k=1:length(y)
            for l=1:length(z)
                  wr2=[ancho, y(k),z(l)];
                  d1=sqrt(dot(rt4-wr2,rt4-wr2));
                  cos_phi=dot(nt,wr2-rt4)/d1;
                  cos_alpha=dot(VN2,rt4-wr2)/d1;
                  d2=sqrt(dot(wr2-r,wr2-r));
                  cos_beta=dot(VN2,r-wr2)/d2;
                  cos_tetha=dot(nr,wr2-r)/d2;
                  
                if abs(acosd(cos_tetha)) <= FOV 
                   H2(i,j)=H2(i,j)+((m+1)*p2*Aef*dA1*(cos_phi^m)*(cos_alpha)*(cos_beta)*(cos_tetha))/(2*pi^2*d1^2*d2^2);
                 elseif abs(acosd(cos_tetha)) > FOV 
                      H2(i,j)=H2(i,j);
                end
            end
           end
        
 %        Pared 3
        
        for k=1:length(x)
            for l=1:length(z)
                  wr3=[x(k),0, z(l)];
                  d1=sqrt(dot(rt4-wr3,rt4-wr3));
                  cos_phi=dot(nt,wr3-rt4)/d1;
                  cos_alpha=dot(VN3,rt4-wr3)/d1;
                  d2=sqrt(dot(wr3-r,wr3-r));
                  cos_beta=dot(VN3,r-wr3)/d2;
                  cos_tetha=dot(nr,wr3-r)/d2;
                 
                 if abs(acosd(cos_tetha)) <= FOV 
                    H3(i,j)=H3(i,j)+((m+1)*p*Aef*dA1*(cos_phi^m)*(cos_alpha)*(cos_beta)*(cos_tetha))/(2*pi^2*d1^2*d2^2);
                 elseif abs(acosd(cos_tetha)) > FOV 
                      H3(i,j)=H3(i,j);
                  end
            end
        end
        
        %Pared 4
        
        for k=1:length(x)
            for l=1:length(z)
                  wr4=[x(k),largo, z(l)];
                  d1=sqrt(dot(rt4-wr4,rt4-wr4));
                  cos_phi=dot(nt,wr4-rt4)/d1;
                  cos_alpha=dot(VN4,rt4-wr4)/d1;
                  d2=sqrt(dot(wr4-r,wr4-r));
                  cos_beta=dot(VN4,r-wr4)/d2;
                  cos_tetha=dot(nr,wr4-r)/d2;
                 
                  if abs(acosd(cos_tetha)) <= FOV 
                    H4(i,j)=H4(i,j)+((m+1)*p*Aef*dA1*(cos_phi^m)*(cos_alpha)*(cos_beta)*(cos_tetha))/(2*pi^2*d1^2*d2^2);
                  
                  elseif abs(acosd(cos_tetha)) > FOV 
                      H4(i,j)=H4(i,j);
                   end
                  
               

            end
        end
        
    end
end

H=H1+H2+H3;
P_rec_A4=(HLos+H)*P_total.*Ts.*G_Con;


%%
P_rec_total_1ref=P_rec_A1+P_rec_A2 +P_rec_A3 +P_rec_A4;
P_rec_1ref_dBm=10*log10(P_rec_total_1ref);

figure(1)
grid on
surf(x,y,P_rec_1ref_dBm');
xlabel('Ancho [m]')
ylabel('Largo [m]')
title('Red Li-Fi','FontSize',12)
box on
colorbar

