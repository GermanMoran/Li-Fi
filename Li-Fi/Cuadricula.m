%Codigo para calcular la distribución de Potencia.

function [Prx_dBLos,x,y]=Cuadricula(largo,ancho,alto,hpr,PxTx,PyTx)
    %Dimenciones de la Habitacion 
    %largo=5 ; ancho=5; alto=3; 

    %Parametros del Transmisor
        rt=[PxTx,PyTx,alto];                   %Posicion del transmisor
        phi_g = 60;                         %angulo de media potencia en grados
        phi_r = phi_g*pi/180;               %angulo de media potencia radianes
        m = -log(2)/(log(cos(phi_r)));      %Grado  del lobulo lambertiano
        Pled=236;                           %Potencia optica media maxima(mw)
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
        %dA=h*largo/(Ny*Nz);                         %Calculo diferencial de area
        x=0:ancho/Nx:ancho;
        y=0:largo/Ny:largo;
        z=hpr:alto/Nz:alto;

        %Plano Recepecion
        [X,Y,Z]=meshgrid(x,y,hpr);
        %p=0.5;                                          %Coeficiente de reflexion de todas las paredes


       
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


            end
        end

        Prx_dBLos=10*log10((HLos)*Pled.*Ts.*G_Con);

    
end




