
function [Drms,x,y]= prueba(largo,ancho,alto,hpr,p,PxTx,PyTx)

%Funcion prueba 
%Tiempo en ns (Depuracion del programa
C=3e8*1e-9;
%Parametros iniciales (angulo media potencia, indice de refleccion paredes)
theta=60;
m=-log10(2)/log10(cosd(theta));
Aef=7.02e-3;
p1=p;
p2=p;
p3=p;
p4=p;

Ts=1;
%index=1.5;
FOV=60;


h_tx_rx=alto-hpr;

Nx=ancho*10; Ny=largo*10; Nz=round(3*alto);
%dA=17.4e-3;
dA=h_tx_rx*largo/((Ny-1)*(Nz-1));
%Construccion de la grilla

x=linspace(0,ancho,Nx);
y=linspace(0,largo,Ny);
z=linspace(hpr,alto,Nz);

% x=0:ancho/Nx:ancho;
% y=0:largo/Ny:largo;
% z=hpr:alto/Nz:round(alto);

[X,Y,Z]=meshgrid(x,y,z);


%Posición del transmisor
rt=[PxTx PyTx alto];
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
delta_t=1/2;


for ii=1:length(x)
    for jj=1:length(y)
     r=[x(ii) y(jj) hpr];
     t_vector=0:100/(delta_t); % time vector in ns
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
                     h_vector(index)=h_vector(index)+(m+1)*p1*Aef*dA*(cos_phi^m)*(cos_alpha)*(cos_beta)*(cos_tetha1)/(2*pi^2*d1^2*d2^2);
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
                     h_vector(index)=h_vector(index)+(m+1)*p2*Aef*dA*(cos_phi^m)*(cos_alpha)*(cos_beta)*(cos_tetha1)/(2*pi^2*d1^2*d2^2);
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
                     h_vector(index)=h_vector(index)+(m+1)*Aef*p3*dA*cos_phi^m*cos_alpha*cos_beta*cos_tetha1/(2*pi^2*d1^2*d2^2);
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
                     h_vector(index)=h_vector(index)+(m+1)*p4*Aef*dA*(cos_phi^m)*(cos_alpha)*(cos_beta)*(cos_tetha1)/(2*pi^2*d1^2*d2^2);
                     end
                     
                end
             end
             t_vector=t_vector*delta_t;
             mean_delay(ii,jj)=sum((h_vector).^2.*t_vector)./sum(h_vector.^2);
             Drms(ii,jj)=sqrt(sum((t_vector-mean_delay(ii,jj)).^2.*h_vector.^2)./sum(h_vector.^2));
             Rb(ii,jj)=(1/(10*Drms(ii,jj)));
    end
end

end