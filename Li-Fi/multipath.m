
function [Prx_dBLos,Prx_dBHnlos,Prx_dBHtot,x,y]=multipath(largo,ancho,alto,p,hpr,PxTx,PyTx)
    %Dimenciones de la Habitacion 
    %largo=5 ; ancho=5; alto=3; 

    %Parametros del Transmisor
        rt=[PxTx,PyTx,alto];                %Posicion del transmisor
        phi_g = 60;                         %angulo de media potencia en grados
        phi_r = phi_g*pi/180;               %angulo de media potencia radianes
        m = -log(2)/(log(cos(phi_r)));      %Grado  del lobulo lambertiano
        Pled=480;                           %Potencia optica media maxima(mw)
        nt=[0,0,-1];                        %Vector normal del Transmisor

        %Parametros Receptor 
        %hpr=0.85;            %Altura del plano receptor donde se encuentran las mesas(m)
        h=alto-hpr;         %Altura entre el plano transmisor y receptor
        Aef=7.02E-3;           %Area efectiva el fotodetector
        Ts=1;               %Ganancia del Filtro Optico
        index=1.5;          %Indice de refraccion del lente del fotodetector
        FOV=60;             %Angulo de máxima incidencia(grados)
        nr=[0,0,1];         %Vector normal del receptor
        G_Con=(index.^2)/(sind(FOV)).^2;             %Ganacia del concentrador Optico sin imagenes

        %Superficies
        Nx=ancho*10; Ny=largo*10; Nz=round(alto*10);   %numero de grillas en c/d superficie
        dA=h*largo/(Ny*Nz);                         %Calculo diferencial de area
        x=0:ancho/Nx:ancho;
        y=0:largo/Ny:largo;
        z=hpr:alto/Nz:alto;

        %Plano Recepecion
        [X,Y,Z]=meshgrid(x,y,hpr);
        %p=0.8;                                          %Coeficiente de reflexion de todas las paredes


        VN1=[1,0,0];
        VN2=[-1,0,0];
        VN3=[0,-1,0];                                   %Vectores Normales de cada pared
        VN4=[0,1,0];


        dA1=(largo*h)/((Ny-1)*(Nz-1));                   %Diferencial de area de cada pared

        for i=1:length(x)
            for j=1:length(y)
                %Calculo distribucion de potencia cuando existe LOS
                r=[x(i), y(j), hpr];                     %Vector posicion del receptor
                d=sqrt(dot(rt-r,rt-r));                  %Distancia entre el trasmisor y receptor
                cos_phi=dot(nt,(r-rt))/d;                %Angulo de irradianza
                cos_tetha=dot(nr,(rt-r))/d;              %Angulo de incidencia

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
                 A(i,j)=0;

                %Pared 1
                for k=1:length(y)
                    for l=1:length(z)
                        wr=[0,y(k),z(l)] ;             %Posicion del diferencial de area pared 1
                        d1=sqrt(dot(rt-wr ,rt-wr));    %Distancia entre LED y diferencial de area
                        cos_phi=dot(nt,wr-rt)/d1;      %Angulo entre vector normal del LED y el dA
                        cos_alpha=dot(VN1,rt-wr)/d1;   %Angulo entre vector normal del dA y el LED
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


                %Pared 2

                  for k=1:length(y)
                    for l=1:length(z)
                          wr2=[ancho, y(k),z(l)];
                          d1=sqrt(dot(rt-wr2,rt-wr2));
                          cos_phi=dot(nt,wr2-rt)/d1;
                          cos_alpha=dot(VN2,rt-wr2)/d1;
                          d2=sqrt(dot(wr2-r,wr2-r));
                          cos_beta=dot(VN2,r-wr2)/d2;
                          cos_tetha=dot(nr,wr2-r)/d2;

                        if abs(acosd(cos_tetha)) <= FOV 
                           H2(i,j)=H2(i,j)+((m+1)*p*Aef*dA1*(cos_phi^m)*(cos_alpha)*(cos_beta)*(cos_tetha))/(2*pi^2*d1^2*d2^2);
                         elseif abs(acosd(cos_tetha)) > FOV 
                              H2(i,j)=H2(i,j);
                        end
                    end
                   end

         %        Pared 3

                for k=1:length(x)
                    for l=1:length(z)
                          wr3=[x(k), largo, z(l)];
                          d1=sqrt(dot(rt-wr3,rt-wr3));
                          cos_phi=dot(nt,wr3-rt)/d1;
                          cos_alpha=dot(VN3,rt-wr3)/d1;
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
                          wr4=[x(k), 0, z(l)];
                          d1=sqrt(dot(rt-wr4,rt-wr4));
                          cos_phi=dot(nt,wr4-rt)/d1;
                          cos_alpha=dot(VN4,rt-wr4)/d1;
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

%          H1(1,:)=H1(2,:);
%          H2(length(H2),:)=H2(length(H2)-1,:);        %codigo correguir inderminacion en los bordes de las paredes
%          H3(:,length(H3))=H3(:,length(H3)-1);   
%          H4(:,1)=H4(:,2);

        Hnlos=H1+H2+H3+H4;                             %Gnancia total NLOS
        Htot=HLos+Hnlos;                               %Ganancia total: LOS+NLOS
    
        Prx_dBLos=10*log10((HLos)*Pled.*Ts.*G_Con);
        Prx_dBHnlos=10*log10((Hnlos)*Pled.*Ts.*G_Con);
        Prx_dBHtot=10*log10((Htot)*Pled.*Ts.*G_Con); 
    
    
end




