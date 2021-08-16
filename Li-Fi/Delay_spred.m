
function [Drms,Rb,x,y]=Delay_spred(largo,ancho,alto,hpr,p,PxTx,PyTx)

%Parametros del transmisor
theta=60;
m=-log10(2)/log10(cosd(theta));
rt=[PxTx PyTx 3];
nt=[0 0 -1];

%Parametos del receptor
Aef=7.02e-3;
Ts=1;
FOV=60;
h_tx_rx=alto-hpr;
nr=[0 0 1];


%Numero de grillas
Nx=ancho*10; Ny=largo*10; Nz=round(h_tx_rx*10);
dA=h_tx_rx*largo/(Ny*Nz);

x=linspace(0,ancho,Nx);
y=linspace(0,largo,Ny);
z=linspace(hpr,alto,Nz);

% x=0:ancho/Nx:ancho;
% y=0:largo/Ny:largo;
% z=linspace(hpr,alto,Nz);

[X,Y,Z]=meshgrid(x,y,z);


%Vectores normales de c/d pared
VN1=[1 0 0];
VN2=[-1 0 0];
VN3=[0 -1 0];
VN4=[0 1 0];

%Nuevos paametros 
delta_t=1/2;        %Resolución en ns
C=3e8*1e-9;       %Velocidad de la luz



for ii=1:length(x)
    for jj=1:length(y)
     r=[x(ii) y(jj) hpr];                               %Posición del recpetor
     t_vector=0:100/(delta_t);                           %Vector de tiempos(ns)
     h_vector=zeros(1,length(t_vector));                %vector respuesta impulso del sistema
     d=sqrt(dot(rt-r,rt-r));                            %Distancia entre LED y FD
     cos_phi=dot(nt,(r-rt))/d;                          %Angulo de irradianza e incidencia
     tau0=d/C;                                          %Tiempo del rayo de luz en llegar al Rx
     index=find(round(tau0/delta_t)==t_vector);         %Indice de tiempo
         if abs(acosd(cos_phi))<=FOV
             h_vector(index)=Aef*(m+1).*(cos_phi.^m).*(cos_phi)./(2*pi.*d.^2);    %Respuesta al impulso
         end
         
        
         
        %Reflexión pared 1
        count=1;

            for k=1:length(y)
                for l=1:length(z)
                wr=[0,y(k),z(l)] ;                      %Posición del diferencial de area pared1
                d1=sqrt(dot(rt-wr, rt-wr));             %Distancia entre LED y diferencial de area
                cos_phi=dot(nt,wr-rt)/d1;               %Angulo entre vector normal del LED y el dA
                cos_alpha=dot(VN1,rt-wr)/d1;            %Angulo entre vector normal del dA y el LED
                d2=sqrt(dot(wr-r,wr-r));                %Distancia entre el diferencial de area y el fotodetector
                cos_beta=dot(VN1,r-wr)/d2;              %Angulo entre vector normal del dA y el fotodetector
                cos_tetha1=dot(nr,wr-r)/d2;             %Angulo entre el vector normal del FD  y el dA
                
                tau1=(d1+d2)/C;                         %Tiempo del rayo LED-PARED-FOTODECTOR
                index=find(round(tau1/delta_t)==t_vector);
                     if abs(acosd(cos_tetha1))<=FOV
                     h_vector(index)=h_vector(index)+(m+1)*p*Aef*dA*(cos_phi^m)*(cos_alpha)*(cos_beta)*(cos_tetha1)/(2*pi^2*d1^2*d2^2);
                     end
                count=count+1;
                 end
            end
            
            
            %Reflexión pared 2
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
            
  
            %Reflexión Pared 3
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
             t_vector=t_vector*delta_t;                                                              %Escala de tiempo
             mean_delay(ii,jj)=sum((h_vector).^2.*t_vector)./sum(h_vector.^2);                       %Media del retardo
             Drms(ii,jj)=sqrt(sum((t_vector-mean_delay(ii,jj)).^2.*h_vector.^2)./sum(h_vector.^2));  %Drms
             Rb(ii,jj)=(1/(10*Drms(ii,jj)));
             
           
                    
          
             
    end
end
end

%Graficas del DRMS
% 
% figure(4)
% surf(x,y,Drms);
% colorbar
% axis([0 largo 0 ancho min(min(Drms)) max(max(Drms))]);
% figure(6)
% surf(x,y,Rb)
% colorbar
