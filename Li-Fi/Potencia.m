function varargout = Potencia(varargin)
%POTENCIA MATLAB code file for Potencia.fig
%      POTENCIA, by itself, creates a new POTENCIA or raises the existing
%      singleton*.
%
%      H = POTENCIA returns the handle to a new POTENCIA or the handle to
%      the existing singleton*.
%
%      POTENCIA('Property','Value',...) creates a new POTENCIA using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to Potencia_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      POTENCIA('CALLBACK') and POTENCIA('CALLBACK',hObject,...) call the
%      local function named CALLBACK in POTENCIA.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Potencia

% Last Modified by GUIDE v2.5 04-Feb-2020 00:18:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Potencia_OpeningFcn, ...
                   'gui_OutputFcn',  @Potencia_OutputFcn, ...
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


% --- Executes just before Potencia is made visible.
function Potencia_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for Potencia
handles.output = hObject;


set(handles.axes2,'XColor','W');
set(handles.axes2,'YColor','w')
set(handles.axes2,'GridColor','k')
grid on
guidata(hObject, handles);






% --- Outputs from this function are returned to the command line.
function varargout = Potencia_OutputFcn(hObject, eventdata, handles)
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
    close(Potencia);
    
end
% --- Executes on button press in aco_ofdm.
 
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
 parametrostx(1,8) = str2double(get(handles.xty, 'String')); 
 
 
 
 if isnan(parametrostx(1,1)) || isnan(parametrostx(1,2))|| isnan(parametrostx(1,7)) || isnan(parametrostx(1,8))
    errordlg('¡ERROR! DATOS INVALIDOS O CAMPOS VACÍOS','Error');
 else
    save parametrostx.mat parametrostx
 end
 


 ParametrosCanal=zeros(1,3);
 ParametrosCanal(1,1)=get(handles.awgn,'value');
 ParametrosCanal(1,2)=get(handles.multitrayectoria,'value');
 ParametrosCanal(1,3)=str2double(get(handles.snr, 'String')); 
 
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
[y,data,ceros_fluj]=modulacionACO(datosCargadosTx.parametrostx(1,1),datosCargadosTx.parametrostx(1,2),datosCargadosTx.parametrostx(1,3),datosCargadosTx.parametrostx(1,5));
datosCargadosTx.parametrostx(1,5)
tam_y=size(y);  
scatterplot(reshape(y,1,tam_y(1)*tam_y(2)));

save datosmodulados.mat y data ceros_fluj

% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
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



axes(handles.axes2)
stem(hermi_real,'b','filled')
title('Parte real e imaginaria ACO-OFDM')
xlabel('Portadoras')
ylabel('Amplitud')
grid on
hold on
stem(hermi_imag,'r','filled')
legend('Parte real','Parte Imaginaria')



% --- Executes on button press in P.
function P_Callback(hObject, eventdata, handles)
% hObject    handle to P (see GCBO)
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


[Smod_tx_DC,val_pot]=ConversionDAC(datosCargadosTx.parametrostx(1,4),sim_serie,datosCargadosTx.parametrostx(1,1),datosCargadosTx.parametrostx(1,6));


axes(handles.axes2);
plot(handles.axes2,Smod_tx_DC(1,1:datosCargadosTx.parametrostx(1,1)*4*datosCargadosTx.parametrostx(1,6)),'b')
grid on;
title(handles.axes2,'Señal de voltaje inicial')
xlabel(handles.axes2,'Tiempo (ns)')
ylabel(handles.axes2,'Voltaje (v)')




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
%figure(3)
limpiar();
axes(handles.axes2);
stem(sim_serie_real(1,1:256),'b','filled');
title('Parte real e imaginaria ACO-OFDM')
xlabel('Tiempo de simbolo')
ylabel('Amplitud')
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


sim_serie=recortecero(ModulacionACO.ofdm_hermi);

[Smod_tx_DC,val_pot]=ConversionDAC(datosCargadosTx.parametrostx(1,4),sim_serie,datosCargadosTx.parametrostx(1,1),datosCargadosTx.parametrostx(1,6));

limpiar();
axes(handles.axes2)
plot(handles.axes2,val_pot(1,1:datosCargadosTx.parametrostx(1,1)*4*datosCargadosTx.parametrostx(1,6)),'g')
grid on;
title(handles.axes2,'Señal de potencia optica')
xlabel(handles.axes2,'Tiempo 4(ns)')
ylabel(handles.axes2,'potencia optica (lm/W)')

save LED.mat val_pot



% --- Executes on button press in pushbutton16.
function pushbutton16_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton16 (see GCBO)

limpiar();
datosCargadosTx=load('parametrostx.mat');
datosEsenario=load('ParametrosEsenario.mat');


[Prx_dBLos,Prx_dBHnlos,Prx_dBHtot,x,y]= multipath(datosEsenario.ParametrosEsenario(1,1),datosEsenario.ParametrosEsenario(1,2),datosEsenario.ParametrosEsenario(1,3),datosEsenario.ParametrosEsenario(1,4),datosEsenario.ParametrosEsenario(1,5));


figure(1)
contourf(x,y,Prx_dBLos)
xlabel('X (m)')
ylabel('Y (m)')
zlabel(texlabel('Potencia Recibida (dBm)'))
box on
colorbar

%Seleccionar el punto para calcular el perfil de retardos/potencias
punto_rtaimpulso()
[l,m] = getpts(figure(1));
%Puntos redondeados
lr=round(l(1),1);
mr=round(m(1),1);
posicion=[lr mr datosEsenario.ParametrosEsenario(1,5)]


close(figure(1))

axes(handles.axes2);
[Drms,x,y,t_vector,h_vector]=drms(datosEsenario.ParametrosEsenario(1,1),datosEsenario.ParametrosEsenario(1,2),datosEsenario.ParametrosEsenario(1,3),datosEsenario.ParametrosEsenario(1,4),datosEsenario.ParametrosEsenario(1,5),posicion,handles.axes2);


save rtaimpulso.mat  t_vector h_vector


% --- Executes on button press in pushbutton17.
function pushbutton17_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


ParametrosTx=load('parametrostx.mat');
rtaimpulso=load('rtaimpulso.mat');
LED=load('LED.mat');

[xt pot_rx]=fotodetector(rtaimpulso.t_vector,rtaimpulso.h_vector,LED.val_pot,ParametrosTx.parametrostx(1,1));


figure(1)
plot(xt,'b')
grid on
xlabel('Tiempo (ns)')
ylabel('Señal de Voltaje en Rx (V)')



function pushbutton20_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


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

signal_rx_discreta=ADC(rx_ruido.rx_ruido,ParametrosTx.parametrostx(1,6));

save SxADC.mat signal_rx_discreta


figure(1)
stem(signal_rx_discreta(1,1:ParametrosTx.parametrostx(1,1)*4),'b')
grid on
xlabel('Tiempo')
ylabel('Señal de Voltaje discreta en Rx')


% --- Executes on button press in pushbutton22.
function pushbutton22_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton22 (see GCBO)
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
% --- Executes on button press in pushbutton25.
function pushbutton25_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton25 (see GCBO)
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
ParametrosCanal=('ParametrosCanal.mat');
ParametrosTx=load('parametrostx.mat');
rtaimpulso=load('rtaimpulso.mat');
LED=load('LED.mat');
[xt pot_rx]=fotodetector(rtaimpulso.t_vector,rtaimpulso.h_vector,LED.val_pot,ParametrosTx.parametrostx(1,6));

rx_ruido=ruidoawgn(ParametrosCanal.ParametrosCanal(1,3),pot_rx);

save rx_ruido.mat rx_ruido
figure(1) 
plot(rx_ruido,'b')
grid on
xlabel('Tiempo (ns)')
ylabel('Señal de Voltaje en Rx (V)')

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
close(Potencia);



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
set(handles.xty, 'String','2.5');


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
set(handles.xty, 'String','');

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
function Untitled_8_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function potencia_Callback(hObject, eventdata, handles)
% hObject    handle to potencia (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


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


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1


% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2


% --- Executes on button press in awgn.


% --- Executes on button press in nlos.
function nlos_Callback(hObject, eventdata, handles)
% hObject    handle to nlos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of nlos


% --- Executes on button press in graficar.
function graficar_Callback(hObject, eventdata, handles)
% hObject    handle to graficar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



los=get(handles.awgn,'value');
nlos=get(handles.multitrayectoria,'value');
if  los==0 && nlos==0
    errordlg('¡ERROR! DATOS INVALIDOS O CAMPOS VACÍOS','Error');
else
    ParametrosTransmisor=load('parametrostx.mat');
    ParametrosEsenario=load('ParametrosEsenario.mat');
    [Prx_dBLos,Prx_dBHnlos,Prx_dBHtot,x,y]=multipath(ParametrosEsenario.ParametrosEsenario(1,1),ParametrosEsenario.ParametrosEsenario(1,2),ParametrosEsenario.ParametrosEsenario(1,3),ParametrosEsenario.ParametrosEsenario(1,4),ParametrosEsenario.ParametrosEsenario(1,5),ParametrosTransmisor.parametrostx(1,7),ParametrosTransmisor.parametrostx(1,8));
end

% ParametrosEsenario=load('ParametrosEsenario.mat');
%  [Prx_dBLos,Prx_dBHnlos,Prx_dBHtot,x,y]=multipath(ParametrosEsenario.ParametrosEsenario(1,1),ParametrosEsenario.ParametrosEsenario(1,2),ParametrosEsenario.ParametrosEsenario(1,3),ParametrosEsenario.ParametrosEsenario(1,4),ParametrosEsenario.ParametrosEsenario(1,5));

% %----Obtengo awgn valores de awgn botones






if los==1
    axes(handles.axes2)
    surfc(x,y,Prx_dBLos');
    xlabel('X (m)','color','w');
    ylabel('Y (m)','color','w');
    zlabel('Potencia Recibida (dBm)','color','w');
    title('Distribución de potencia LOS','color','w')
    set(handles.axes2,'XColor','W');
    set(handles.axes2,'YColor','w')
    set(handles.axes2,'GridColor','k')

    box on 
    colorbar('color','w')
    rotate3d on
    
end


if nlos==1
    axes(handles.axes2)
    surfc(x,y,Prx_dBHnlos');
    xlabel('X (m)','color','w');
    ylabel('Y (m)','color','w');
    zlabel('Potencia Recibida (dBm)','color','w')
    title('Distribución de potencia NLOS','color','w')
    set(handles.axes2,'XColor','W');
    set(handles.axes2,'YColor','w')
    set(handles.axes2,'GridColor','k');
    box on 
    colorbar('color','w')
    rotate3d on
end


if los==1 && nlos==1
    axes(handles.axes2)
    surfc(x,y,Prx_dBHtot');
    xlabel('X (m)','color','w');
    ylabel('Y (m)','color','w');
    zlabel('Potencia Recibida (dBm)','color','w')
    title('Distribución de potencia Total','color','w')
    set(handles.axes2,'XColor','W');
    set(handles.axes2,'YColor','w')
    set(handles.axes2,'GridColor','k');
    box on 
    colorbar('color','w')
    rotate3d on
end


% --- Executes on button press in checkbox6.
function checkbox6_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox6


% --- Executes on button press in checkbox7.
function checkbox7_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox7


% --- Executes on button press in aco_ofdm.
function aco_ofdm_Callback(hObject, eventdata, handles)
% hObject    handle to aco_ofdm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global boton1
boton1 = get(handles.aco_ofdm,'value');
if boton1==1
    sistemaLi_Fi;
    close(Potencia);
    
end


% --- Executes on button press in drms.
function drms_Callback(hObject, eventdata, handles)
% hObject    handle to drms (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ParametrosTransmisor=load('parametrostx.mat');
ParametrosEsenario=load('ParametrosEsenario.mat');

[Drms,x,y]=prueba(ParametrosEsenario.ParametrosEsenario(1,1),ParametrosEsenario.ParametrosEsenario(1,2),ParametrosEsenario.ParametrosEsenario(1,3),ParametrosEsenario.ParametrosEsenario(1,4),ParametrosEsenario.ParametrosEsenario(1,5),ParametrosTransmisor.parametrostx(1,7),ParametrosTransmisor.parametrostx(1,8));
% 
    axes(handles.axes2)
    surfc(y,x,Drms);
    xlabel('X (m)','color','w');
    ylabel('Y (m)','color','w');
    zlabel('Drms (ns)','color','w')
    title('DRMS','color','w')
    set(handles.axes2,'XColor','W');
    set(handles.axes2,'YColor','w')
    set(handles.axes2,'GridColor','k');
    box on 
    colorbar('color','w')
    rotate3d on
    axis([0 ParametrosEsenario.ParametrosEsenario(1,1) 0 ParametrosEsenario.ParametrosEsenario(1,2) min(min(Drms)) max(max(Drms))]);
    
    
% --- Executes on button press in pushbutton37.
function pushbutton37_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton37 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ParametrosTransmisor=load('parametrostx.mat');
ParametrosEsenario=load('ParametrosEsenario.mat');
[Drms,Rb,x,y]=Delay_spred(ParametrosEsenario.ParametrosEsenario(1,1),ParametrosEsenario.ParametrosEsenario(1,2),ParametrosEsenario.ParametrosEsenario(1,3),ParametrosEsenario.ParametrosEsenario(1,4),ParametrosEsenario.ParametrosEsenario(1,5),ParametrosTransmisor.parametrostx(1,7),ParametrosTransmisor.parametrostx(1,8));
% 
    axes(handles.axes2)
    surfc(x,y,Rb);
    xlabel('X (m)','color','w');
    ylabel('Y (m)','color','w');
    zlabel('Maxima Tasa Teorica (Gbps)','color','w')
    title('Rate Teorica','color','w')
    set(handles.axes2,'XColor','W');
    set(handles.axes2,'YColor','w')
    set(handles.axes2,'GridColor','k');
    box on 
    colorbar('color','w')
    rotate3d on
    axis([0 ParametrosEsenario.ParametrosEsenario(1,1) 0 ParametrosEsenario.ParametrosEsenario(1,2) min(min(Drms)) max(max(Drms))]);



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



function xty_Callback(hObject, eventdata, handles)
% hObject    handle to xty (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xty as text
%        str2double(get(hObject,'String')) returns contents of xty as a double


% --- Executes during object creation, after setting all properties.
function xty_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xty (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function uipushtool6_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
DCO_OFDM;
close(Potencia);


% --------------------------------------------------------------------
function uipushtool8_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
open('Anexo B.pdf');
