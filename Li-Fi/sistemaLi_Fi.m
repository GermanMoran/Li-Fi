function varargout = sistemaLi_Fi(varargin)
%SISTEMALI_FI MATLAB code file for sistemaLi_Fi.fig
%      SISTEMALI_FI, by itself, creates a new SISTEMALI_FI or raises the existing
%      singleton*.
%
%      H = SISTEMALI_FI returns the handle to a new SISTEMALI_FI or the handle to
%      the existing singleton*.
%
%      SISTEMALI_FI('Property','Value',...) creates a new SISTEMALI_FI using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to sistemaLi_Fi_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      SISTEMALI_FI('CALLBACK') and SISTEMALI_FI('CALLBACK',hObject,...) call the
%      local function named CALLBACK in SISTEMALI_FI.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help sistemaLi_Fi

% Last Modified by GUIDE v2.5 03-Feb-2020 23:59:51

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @sistemaLi_Fi_OpeningFcn, ...
                   'gui_OutputFcn',  @sistemaLi_Fi_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
   gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before sistemaLi_Fi is made visible.
function sistemaLi_Fi_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for sistemaLi_Fi
handles.output = hObject;


set(handles.axes2,'xcolor','w') % Cambias el color de blanco a negro
set(handles.axes2,'ycolor','w') % Cambias el color de blanco a negro
set(handles.PS1,'enable','off');
set(handles.Hermitica,'enable','off');
set(handles.PS2,'enable','off');
set(handles.PS3,'enable','off');
set(handles.PS4,'enable','off');


set(handles.axes2,'XColor','W');
set(handles.axes2,'YColor','w')
set(handles.axes2,'GridColor','k')
grid on


guidata(hObject, handles);






% --- Outputs from this function are returned to the command line.
function varargout = sistemaLi_Fi_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in inicio.
function inicio_Callback(hObject, eventdata, handles)
% hObject    handle to inicio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global boton1
boton1 = get(handles.inicio,'value');
if boton1==1
     Inicio;
    close(sistemaLi_Fi);
    
end

% --- Executes on button press in aco_ofdm.
function aco_ofdm_Callback(hObject, eventdata, handles)
% hObject    handle to aco_ofdm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in dco_ofdm.
function dco_ofdm_Callback(hObject, eventdata, handles)
% hObject    handle to dco_ofdm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global boton1
boton1 = get(handles.dco_ofdm,'value');
if boton1==1
     DCO_OFDM;
    close(sistemaLi_Fi);
    
end

% --- Executes on button press in guardar.
function guardar_Callback(hObject, eventdata, handles)
% hObject    handle to guardar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Parametro del esenario de simulación
ParametrosEsenario=zeros(1,5);
ParametrosEsenario(1,1)=str2double(get(handles.large, 'String'));
ParametrosEsenario(1,2)=str2double(get(handles.ancho, 'String'));
ParametrosEsenario(1,3)=str2double(get(handles.alto, 'String'));
ParametrosEsenario(1,4)=str2double(get(handles.cr, 'String'));
ParametrosEsenario(1,5)=str2double(get(handles.hpr, 'String'));

 if isnan(ParametrosEsenario(1,1)) || isnan(ParametrosEsenario(1,2)) || isnan(ParametrosEsenario(1,3)) || isnan(ParametrosEsenario(1,4)) || isnan(ParametrosEsenario(1,5))
    errordlg('¡ERROR! DATOS INVALIDOS O CAMPOS VACÍOS','Error');
 else
     
    %set(handles.AtenuacionGases,'Enable','on')
    save ParametrosEsenario.mat ParametrosEsenario
 end
 

 
 parametrostx=zeros(1,8);
 parametrostx(1,1) = str2double(get(handles.N_sub, 'String')); 
 parametrostx(1,2) = str2double(get(handles.M, 'String')); 
 parametrostx(1,3) = get(handles.TM,'Value');
 parametrostx(1,4) = get(handles.TL,'Value');
 parametrostx(1,5) = get(handles.RB,'Value');
 parametrostx(1,6) = parametrostx(1,1)*16;
 parametrostx(1,7) = str2double(get(handles.xtx, 'String'));
 parametrostx(1,8) = str2double(get(handles.ytx, 'String')); 
 
 
 
 if isnan(parametrostx(1,1)) || isnan(parametrostx(1,2)) || isnan(parametrostx(1,7)) || isnan(parametrostx(1,8))
    errordlg('¡ERROR! DATOS INVALIDOS O CAMPOS VACÍOS','Error');
 else
    save parametrostx.mat parametrostx
 end
 

 ParametrosCanal=zeros(1,3);
 ParametrosCanal(1,1)=str2double(get(handles.snr, 'String')); 
 
  if isnan( ParametrosCanal(1,3)) 
    errordlg('¡ERROR! DATOS INVALIDOS O CAMPOS VACÍOS','Error');
 else
    save ParametrosCanal.mat ParametrosCanal
  end
  
 
 
 
  
 ParametrosReceptor=zeros(1,3);
 ParametrosReceptor(1,1)=get(handles.FD,'Value');
 save ParametrosReceptor.mat ParametrosReceptor
  
% --- Executes on selection change in FD.
function FD_Callback(hObject, eventdata, handles)
% hObject    handle to FD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns FD contents as cell array
%        contents{get(hObject,'Value')} returns selected item from FD


% --- Executes during object creation, after setting all properties.
function FD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in awgn.
function awgn_Callback(hObject, eventdata, handles)
% hObject    handle to awgn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of awgn


% --- Executes on button press in multitrayectoria.
function multitrayectoria_Callback(hObject, eventdata, handles)
% hObject    handle to multitrayectoria (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of multitrayectoria



function N_sub_Callback(hObject, eventdata, handles)
% hObject    handle to N_sub (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of N_sub as text
%        str2double(get(hObject,'String')) returns contents of N_sub as a double


% --- Executes during object creation, after setting all properties.
function N_sub_CreateFcn(hObject, eventdata, handles)
% hObject    handle to N_sub (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function M_Callback(hObject, eventdata, handles)
% hObject    handle to M (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of M as text
%        str2double(get(hObject,'String')) returns contents of M as a double


% --- Executes during object creation, after setting all properties.
function M_CreateFcn(hObject, eventdata, handles)
% hObject    handle to M (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in TM.
function TM_Callback(hObject, eventdata, handles)
% hObject    handle to TM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns TM contents as cell array
%        contents{get(hObject,'Value')} returns selected item from TM


% --- Executes during object creation, after setting all properties.
function TM_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in TL.
function TL_Callback(hObject, eventdata, handles)
% hObject    handle to TL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns TL contents as cell array
%        contents{get(hObject,'Value')} returns selected item from TL


% --- Executes during object creation, after setting all properties.
function TL_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in potencia.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to potencia (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in modulacion.
function modulacion_Callback(hObject, eventdata, handles)
% hObject    handle to modulacion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

datosCargadosTx=load('parametrostx.mat');
[y,data,ceros_fluj]=modulacionACO(datosCargadosTx.parametrostx(1,1),datosCargadosTx.parametrostx(1,2),datosCargadosTx.parametrostx(1,3));
tam_y=size(y);  
scatterplot(reshape(y,1,tam_y(1)*tam_y(2)));

save datosmodulados.mat y data ceros_fluj

% --- Executes on button press in PS1.
function PS1_Callback(hObject, eventdata, handles)
% hObject    handle to PS1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --- Executes on button press in Hermitica.
function Hermitica_Callback(hObject, eventdata, handles)
% hObject    handle to Hermitica (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes on button press in IFFT.
function IFFT_Callback(hObject, eventdata, handles)
% hObject    handle to IFFT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

limpiar()
datosmodulados=load('datosmodulados.mat');
ofdm_hermi=Bloqueifft(datosmodulados.y);
size(ofdm_hermi);
hermi_real = real(ofdm_hermi(:,1));
hermi_imag = imag (ofdm_hermi(:,1));

save ModulacionACO.mat  ofdm_hermi



% set(handles.axes2,'XColor','W');
% set(handles.axes2,'YColor','w');
axes(handles.axes2)
stem(hermi_real,'b','filled')
title('Parte real e imaginaria ACO-OFDM','color','w')
xlabel('Subportadoras','color','w')
ylabel('Amplitud','color','w')
set(handles.axes2,'XColor','W');
set(handles.axes2,'YColor','w')
set(handles.axes2,'GridColor','k')
grid on
hold on
stem(hermi_imag,'r','filled')
legend('Parte real','Parte Imaginaria')



% --- Executes on button press in PS2.
function PS2_Callback(hObject, eventdata, handles)
% hObject    handle to PS2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in DAC.
function DAC_Callback(hObject, eventdata, handles)
% hObject    handle to DAC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

datosCargadosTx=load('parametrostx.mat');
ModulacionACO=load('ModulacionACO.mat');



sim_serie=recortecero(ModulacionACO.ofdm_hermi);

[Smod_tx_DC,val_pot,Es]=ConversionDAC(datosCargadosTx.parametrostx(1,4),sim_serie,datosCargadosTx.parametrostx(1,1),datosCargadosTx.parametrostx(1,6));

save ConversionDAC.mat Smod_tx_DC val_pot

axes(handles.axes2);
plot(handles.axes2,Smod_tx_DC(1,1:datosCargadosTx.parametrostx(1,1)*4*datosCargadosTx.parametrostx(1,6)),'b')
grid on;
title(handles.axes2,'Señal de voltaje inicial','color','w')
xlabel(handles.axes2,'Tiempo','color','w')
ylabel(handles.axes2,'Voltaje (v)','color','w')
set(handles.axes2,'XColor','W');
set(handles.axes2,'YColor','w')
set(handles.axes2,'GridColor','k')




% --- Executes on button press in recorte_cero.
function recorte_cero_Callback(hObject, eventdata, handles)
% hObject    handle to recorte_cero (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Patametrostx=load('parametrostx.mat');
ModulacionACO=load('ModulacionACO.mat');

sim_serie=recortecero(ModulacionACO.ofdm_hermi);
sim_serie_real = real(sim_serie(1,1:Patametrostx.parametrostx(1,1)*4));
sim_serie_imag= imag (sim_serie(1,1:Patametrostx.parametrostx(1,1)*4));


size(sim_serie);


limpiar();
axes(handles.axes2);
stem(sim_serie_real(1,1:Patametrostx.parametrostx(1,1)*4),'b','filled');
title('Parte real e imaginaria ACO-OFDM','color','w')
xlabel('Subportadoras','color','w')
ylabel('Amplitud','color','w')
set(handles.axes2,'XColor','W');
set(handles.axes2,'YColor','w')
set(handles.axes2,'GridColor','k')
grid on
hold on
stem(sim_serie_imag,'r','filled');
legend('Parte real','Parte Imaginaria')
hold off

% --- Executes on button press in LED.
function LED_Callback(hObject, eventdata, handles)
% hObject    handle to LED (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

datosCargadosTx=load('parametrostx.mat');
ModulacionACO=load('ModulacionACO.mat');
ConversionDAC=load('ConversionDAC.mat');

limpiar();
axes(handles.axes2)
plot(handles.axes2,ConversionDAC.val_pot(1,1:datosCargadosTx.parametrostx(1,1)*4*datosCargadosTx.parametrostx(1,6)),'g')
grid on;
title(handles.axes2,'Señal de potencia optica','color','w')
xlabel(handles.axes2,'Tiempo','color','w')
ylabel(handles.axes2,'potencia optica (lm/W)','color','w')
set(handles.axes2,'XColor','W');
set(handles.axes2,'YColor','w')
set(handles.axes2,'GridColor','k')




% --- Executes on button press in pushbutton16.
function pushbutton16_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton16 (see GCBO)

limpiar();
parametrostx=load('parametrostx.mat');
datosEsenario=load('ParametrosEsenario.mat');


[Prx_dBLos,x,y]= Cuadricula(datosEsenario.ParametrosEsenario(1,1),datosEsenario.ParametrosEsenario(1,2),datosEsenario.ParametrosEsenario(1,3),datosEsenario.ParametrosEsenario(1,5),parametrostx.parametrostx(1,7),parametrostx.parametrostx(1,8));


figure(1)
contourf(x,y,Prx_dBLos')
xlabel('X (m)')
ylabel('Y (m)')
zlabel(texlabel('Potencia Recibida (dBm)'))
box on
%colorbar

%Seleccionar el punto para calcular el perfil de retardos/potencias
punto_rtaimpulso()
[l,m] = getpts(figure(1));
%Puntos redondeados
lr=round(l(1),1);
mr=round(m(1),1);
posicion=[lr mr datosEsenario.ParametrosEsenario(1,5)]


close(figure(1))

axes(handles.axes2);
[Drms,x,y,t_vector,h_vector]=drms(datosEsenario.ParametrosEsenario(1,1),datosEsenario.ParametrosEsenario(1,2),datosEsenario.ParametrosEsenario(1,3),datosEsenario.ParametrosEsenario(1,4),datosEsenario.ParametrosEsenario(1,5),posicion,handles.axes2,parametrostx.parametrostx(1,7),parametrostx.parametrostx(1,8));

set(handles.axes2,'XColor','W');
set(handles.axes2,'YColor','w')
set(handles.axes2,'GridColor','k')
grid on 
title('Perfil de potencias y retardos','color','w')
xlabel(handles.axes2,'Tiempo','color','w')
ylabel(handles.axes2,'Amplitud','color','w')
save rtaimpulso.mat  t_vector h_vector


% --- Executes on button press in pushbutton17.
function pushbutton17_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


ParametrosTx=load('parametrostx.mat');
rtaimpulso=load('rtaimpulso.mat');
ConversionDAC=load('ConversionDAC.mat');

[xt pot_rx]=fotodetector(rtaimpulso.t_vector,rtaimpulso.h_vector,ConversionDAC.val_pot,ParametrosTx.parametrostx(1,1),ParametrosTx.parametrostx(1,2),ParametrosTx.parametrostx(1,5));

save VariablesFD.mat xt pot_rx
axes(handles.axes2)
plot(xt,'b')
grid on
xlabel('Tiempo(ns)','color','w')
ylabel('Señal de Voltaje en Rx (V)','color','w')
set(handles.axes2,'XColor','W');
set(handles.axes2,'YColor','w')
set(handles.axes2,'GridColor','k')



function pushbutton20_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

limpiar();
ParametrosTx=load('parametrostx.mat');
bloquefft=load('bloquefft.mat');
datosmodulados=load('datosmodulados.mat');
 
 sim_rx_par_des=decisor(bloquefft.sim_rx_par, ParametrosTx.parametrostx(1,1));
 sim_rx_serie= demodulacion(sim_rx_par_des,datosmodulados.y,ParametrosTx.parametrostx(1,2),datosmodulados.ceros_fluj) ;
 bits_errados=(sim_rx_serie==datosmodulados.data);
 bits_errados=sum(~bits_errados);
 num_datos = length(datosmodulados.data);
 ber=bits_errados/num_datos;

set(handles.ber,'String',strcat('BER  = ',num2str(ber)));



%Obtengo la Fotocorriente en el receptor
% 
% datosparametrosCanal=load('parametrosCanal.mat');
% datosparametrosLed=load('ParametrosLed.mat');
%pot_rx=fotodetector(datosparametrosCana.parametrosCanal(1,1),datosparametrosCanal.parametrosCanal(1,3),M,datosparametrosLed.ParametrosLed(1,2))



% --- Executes on button press in pushbutton21.
function pushbutton21_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%Conversión analoga/Digital
ParametrosTx=load('parametrostx.mat');
rx_ruido=load('rx_ruido.mat');

signal_rx_discreta=ADC(rx_ruido.rx_ruido,ParametrosTx.parametrostx(1,6))

save SxADC.mat signal_rx_discreta


axes(handles.axes2)
stem(signal_rx_discreta(1,1:ParametrosTx.parametrostx(1,1)*4),'b')
grid on
xlabel('Tiempo','color','w')
ylabel('Señal de Voltaje discreta en Rx','color','w')
set(handles.axes2,'XColor','W');
set(handles.axes2,'YColor','w')
set(handles.axes2,'GridColor','k') 



% --- Executes on button press in PS3.
function PS3_Callback(hObject, eventdata, handles)
% hObject    handle to PS3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in FFT.
function FFT_Callback(hObject, eventdata, handles)
% hObject    handle to FFT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


ParametrosTx=load('parametrostx.mat');
SxADC=load('SxADC.mat');
ModulacionACO=load('ModulacionACO.mat');
datosmodulados=load('datosmodulados.mat');

sim_rx_par= Bloquefft(SxADC.signal_rx_discreta,ParametrosTx.parametrostx(1,1),ModulacionACO.ofdm_hermi,datosmodulados.y);

scatterplot(sim_rx_par)

save bloquefft.mat sim_rx_par

% --- Executes on button press in pushbutton24.
function pushbutton24_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

bloquefft=load('bloquefft.mat');
% --- Executes on button press in PS4.
function PS4_Callback(hObject, eventdata, handles)
% hObject    handle to PS4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function large_Callback(hObject, eventdata, handles)
% hObject    handle to large (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of large as text
%        str2double(get(hObject,'String')) returns contents of large as a double


% --- Executes during object creation, after setting all properties.
function large_CreateFcn(hObject, eventdata, handles)
% hObject    handle to large (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ancho_Callback(hObject, eventdata, handles)
% hObject    handle to ancho (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ancho as text
%        str2double(get(hObject,'String')) returns contents of ancho as a double


% --- Executes during object creation, after setting all properties.
function ancho_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ancho (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function alto_Callback(hObject, eventdata, handles)
% hObject    handle to alto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of alto as text
%        str2double(get(hObject,'String')) returns contents of alto as a double


% --- Executes during object creation, after setting all properties.
function alto_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cr_Callback(hObject, eventdata, handles)
% hObject    handle to cr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cr as text
%        str2double(get(hObject,'String')) returns contents of cr as a double


% --- Executes during object creation, after setting all properties.
function cr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function hpr_Callback(hObject, eventdata, handles)
% hObject    handle to hpr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hpr as text
%        str2double(get(hObject,'String')) returns contents of hpr as a double


% --- Executes during object creation, after setting all properties.
function hpr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hpr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in RB.
function RB_Callback(hObject, eventdata, handles)
% hObject    handle to RB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns RB contents as cell array
%        contents{get(hObject,'Value')} returns selected item from RB


% --- Executes during object creation, after setting all properties.
function RB_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ruido.
function ruido_Callback(hObject, eventdata, handles)
% hObject    handle to ruido (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%snr=30;
ParametrosCanal=load('ParametrosCanal.mat');
ParametrosTx=load('parametrostx.mat');
VariablesFD=load('VariablesFD.mat');
ConversionDAC=load('ConversionDAC.mat');

rx_ruido=ruidoawgn(ParametrosCanal.ParametrosCanal(1,1),VariablesFD.pot_rx,ConversionDAC.Smod_tx_DC,ParametrosTx.parametrostx(1,2));

save rx_ruido.mat rx_ruido

axes(handles.axes2)
plot(rx_ruido,'b')
grid on
xlabel('Tiempo (ns)','color','w')
ylabel('Señal de Voltaje en Rx (V)','color','w')
set(handles.axes2,'XColor','W');
set(handles.axes2,'YColor','w')
set(handles.axes2,'GridColor','k')

function snr_Callback(hObject, eventdata, handles)
% hObject    handle to snr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of snr as text
%        str2double(get(hObject,'String')) returns contents of snr as a double


% --- Executes during object creation, after setting all properties.
function snr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to snr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Limpiar.
function Limpiar_Callback(hObject, eventdata, handles)
% hObject    handle to Limpiar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --------------------------------------------------------------------
function Inicio_Callback(hObject, eventdata, handles)
% hObject    handle to Inicio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Leds_Callback(hObject, eventdata, handles)
% hObject    handle to Leds (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function resultados_Callback(hObject, eventdata, handles)
% hObject    handle to resultados (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function ayuda_Callback(hObject, eventdata, handles)
% hObject    handle to ayuda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function uipushtool7_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Inicio;
close(sistemaLi_Fi);



function ber_Callback(hObject, eventdata, handles)
% hObject    handle to ber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ber as text
%        str2double(get(hObject,'String')) returns contents of ber as a double


% --- Executes during object creation, after setting all properties.
function ber_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit15_Callback(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit15 as text
%        str2double(get(hObject,'String')) returns contents of edit15 as a double


% --- Executes during object creation, after setting all properties.
function edit15_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function uipushtool5_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.large, 'String','5');
set(handles.ancho, 'String','5');
set(handles.alto, 'String','3');
set(handles.cr, 'String','0.8');
set(handles.hpr, 'String','1');
set(handles.N_sub, 'String','64');
set(handles.M, 'String','4');
set(handles.snr, 'String','25');
set(handles.xtx, 'String','2.5');
set(handles.ytx, 'String','2.5');


% --------------------------------------------------------------------
function uipushtool9_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.large, 'String','');
set(handles.ancho, 'String','');
set(handles.alto, 'String','');
set(handles.cr, 'String','');
set(handles.hpr, 'String','');
set(handles.N_sub, 'String','');
set(handles.M, 'String','');
set(handles.snr, 'String','');
set(handles.xtx, 'String','');
set(handles.ytx, 'String','');

% --------------------------------------------------------------------
function leds_Callback(hObject, eventdata, handles)
% hObject    handle to leds (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function fotodetector_Callback(hObject, eventdata, handles)
% hObject    handle to fotodetector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function esenario_Callback(hObject, eventdata, handles)
% hObject    handle to esenario (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
open('Esenario de simulacion.pdf');

% --------------------------------------------------------------------
function aco_Callback(hObject, eventdata, handles)
% hObject    handle to aco (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function dco_Callback(hObject, eventdata, handles)
% hObject    handle to dco (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
DCO_OFDM
close(sistemaLi_Fi);

       
    
% --------------------------------------------------------------------
function BPW_21_Callback(hObject, eventdata, handles)
% hObject    handle to BPW_21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
open('BPW.pdf');

% --------------------------------------------------------------------
function golden_Callback(hObject, eventdata, handles)
% hObject    handle to golden (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
open('Goldendragon.pdf');

% --------------------------------------------------------------------
function nichia_Callback(hObject, eventdata, handles)
% hObject    handle to nichia (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
open('Nichia.pdf');


% --- Executes on button press in pot.
function pot_Callback(hObject, eventdata, handles)
% hObject    handle to pot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global boton1
boton1 = get(handles.pot,'value');
if boton1==1
    Potencia;
    close(sistemaLi_Fi);
    
end



function xtx_Callback(hObject, eventdata, handles)
% hObject    handle to xtx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xtx as text
%        str2double(get(hObject,'String')) returns contents of xtx as a double


% --- Executes during object creation, after setting all properties.
function xtx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xtx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ytx_Callback(hObject, eventdata, handles)
% hObject    handle to ytx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ytx as text
%        str2double(get(hObject,'String')) returns contents of ytx as a double


% --- Executes during object creation, after setting all properties.
function ytx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ytx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function potencia_Callback(hObject, eventdata, handles)
% hObject    handle to potencia (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Potencia;
close(sistemaLi_Fi);


% --------------------------------------------------------------------
function uipushtool6_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Inicio;
close(sistemaLi_Fi);


% --------------------------------------------------------------------
function uipushtool8_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
open('Anexo B.pdf')
