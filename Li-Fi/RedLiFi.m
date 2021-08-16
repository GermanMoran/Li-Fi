function varargout = RedLiFi(varargin)
%REDLIFI MATLAB code file for RedLiFi.fig
%      REDLIFI, by itself, creates a new REDLIFI or raises the existing
%      singleton*.
%
%      H = REDLIFI returns the handle to a new REDLIFI or the handle to
%      the existing singleton*.
%
%      REDLIFI('Property','Value',...) creates a new REDLIFI using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to RedLiFi_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      REDLIFI('CALLBACK') and REDLIFI('CALLBACK',hObject,...) call the
%      local function named CALLBACK in REDLIFI.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help RedLiFi

% Last Modified by GUIDE v2.5 27-Jan-2020 15:14:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @RedLiFi_OpeningFcn, ...
                   'gui_OutputFcn',  @RedLiFi_OutputFcn, ...
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


% --- Executes just before RedLiFi is made visible.
function RedLiFi_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for RedLiFi
handles.output = hObject;


set(handles.axes2,'XColor','W');
set(handles.axes2,'YColor','w')
set(handles.axes2,'GridColor','k')
grid on
guidata(hObject, handles);






% --- Outputs from this function are returned to the command line.
function varargout = RedLiFi_OutputFcn(hObject, eventdata, handles)
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
    close(RedLiFi);
    
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

%Parametros RED Li-Fi
ParametrosRed=zeros(1,16);
ParametrosRed(1,1)=str2double(get(handles.large, 'String'));
ParametrosRed(1,2)=str2double(get(handles.ancho, 'String'));
ParametrosRed(1,3)=str2double(get(handles.alto, 'String'));
ParametrosRed(1,4)=str2double(get(handles.hpr, 'String'));
ParametrosRed(1,5)=str2double(get(handles.cr, 'String'));
ParametrosRed(1,6)=str2double(get(handles.cr2, 'String'));
ParametrosRed(1,7)=str2double(get(handles.cr3, 'String'));
ParametrosRed(1,8)=str2double(get(handles.cr4, 'String'));
ParametrosRed(1,9) = get(handles.TL,'Value');
ParametrosRed(1,10) = get(handles.FD,'Value');
ParametrosRed(1,11) = get(handles.tx1,'Value');
ParametrosRed(1,12) = get(handles.tx2,'Value');
ParametrosRed(1,13) = get(handles.tx3,'Value');
ParametrosRed(1,14) = get(handles.tx4,'Value');
ParametrosRed(1,15) = get(handles.LineaVista,'Value');
ParametrosRed(1,16) = get(handles.Multipath,'Value');

 if isnan(ParametrosRed(1,1)) || isnan(ParametrosRed(1,2)) || isnan(ParametrosRed(1,3)) || isnan(ParametrosRed(1,4)) || isnan(ParametrosRed(1,5)) || isnan(ParametrosRed(1,6)) || isnan(ParametrosRed(1,7)) || isnan(ParametrosRed(1,8))
    errordlg('¡ERROR! DATOS INVALIDOS O CAMPOS VACÍOS','Error');
 else
     
    %set(handles.AtenuacionGases,'Enable','on')
    save ParametrosRed.mat ParametrosRed
 end
 
ParametrosRed
 

  
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


% --- Executes on button press in redlifi.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to redlifi (see GCBO)
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

set(handles.large, 'String','');
set(handles.ancho, 'String','');
set(handles.alto, 'String','');
set(handles.hpr, 'String','');
set(handles.cr, 'String','');
set(handles.cr2, 'String','');
set(handles.cr3, 'String','');
set(handles.cr4, 'String','');


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
close(RedLiFi);



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

set(handles.large, 'String','12.21');
set(handles.ancho, 'String','12.96');
set(handles.alto, 'String','3.08');
set(handles.hpr, 'String','0.85');
set(handles.cr, 'String','0.6');
set(handles.cr2, 'String','0.06');
set(handles.cr3, 'String','0.6');
set(handles.cr4, 'String','0');



% --------------------------------------------------------------------
function uipushtool9_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.large, 'String','');
set(handles.ancho, 'String','');
set(handles.alto, 'String','');
set(handles.hpr, 'String','');
set(handles.cr, 'String','');
set(handles.cr2, 'String','');
set(handles.cr3, 'String','');
set(handles.cr4, 'String','');


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
% hObject    handle to redlifi (see GCBO)
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



los=get(handles.awgn,'value')
nlos=get(handles.multitrayectoria,'value')
if  los==0 && nlos==0
    errordlg('¡ERROR! DATOS INVALIDOS O CAMPOS VACÍOS','Error');
else
    ParametrosEsenario=load('ParametrosEsenario.mat');
    [Prx_dBLos,Prx_dBHnlos,Prx_dBHtot,x,y]=multipath(ParametrosEsenario.ParametrosEsenario(1,1),ParametrosEsenario.ParametrosEsenario(1,2),ParametrosEsenario.ParametrosEsenario(1,3),ParametrosEsenario.ParametrosEsenario(1,4),ParametrosEsenario.ParametrosEsenario(1,5));
end

% ParametrosEsenario=load('ParametrosEsenario.mat');
%  [Prx_dBLos,Prx_dBHnlos,Prx_dBHtot,x,y]=multipath(ParametrosEsenario.ParametrosEsenario(1,1),ParametrosEsenario.ParametrosEsenario(1,2),ParametrosEsenario.ParametrosEsenario(1,3),ParametrosEsenario.ParametrosEsenario(1,4),ParametrosEsenario.ParametrosEsenario(1,5));

% %----Obtengo awgn valores de awgn botones






if los==1
    axes(handles.axes2)
    surfc(x,y,Prx_dBLos);
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
    surfc(x,y,Prx_dBHnlos);
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
    surfc(x,y,Prx_dBHtot);
    xlabel('X (m)','color','w');
    ylabel('Y (m)','color','w');
    zlabel('Potencia Recibida (dBm)','color','w')
    title('Distribución de potencia Total','color','w')
    set(handles.axes2,'XColor','W');
    set(handles.axes2,'YColor','w')
    set(handles.axes2,'GridColor','k');
    box on 
    colorbar('color','r')
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

ParametrosEsenario=load('ParametrosEsenario.mat');

[Drms,Rb,x,y]=Delay_spred(ParametrosEsenario.ParametrosEsenario(1,1),ParametrosEsenario.ParametrosEsenario(1,2),ParametrosEsenario.ParametrosEsenario(1,3),ParametrosEsenario.ParametrosEsenario(1,4),ParametrosEsenario.ParametrosEsenario(1,5));
% 
    axes(handles.axes2)
    surfc(x,y,Drms);
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

ParametrosEsenario=load('ParametrosEsenario.mat');
[Drms,Rb,x,y]=Delay_spred(ParametrosEsenario.ParametrosEsenario(1,1),ParametrosEsenario.ParametrosEsenario(1,2),ParametrosEsenario.ParametrosEsenario(1,3),ParametrosEsenario.ParametrosEsenario(1,4),ParametrosEsenario.ParametrosEsenario(1,5));
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


% --- Executes on button press in tx1.
function tx1_Callback(hObject, eventdata, handles)
% hObject    handle to tx1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tx1


% --- Executes on button press in tx2.
function tx2_Callback(hObject, eventdata, handles)
% hObject    handle to tx2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tx2


% --- Executes on button press in tx3.
function tx3_Callback(hObject, eventdata, handles)
% hObject    handle to tx3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tx3


% --- Executes on button press in tx4.
function tx4_Callback(hObject, eventdata, handles)
% hObject    handle to tx4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tx4



function cr2_Callback(hObject, eventdata, handles)
% hObject    handle to cr2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cr2 as text
%        str2double(get(hObject,'String')) returns contents of cr2 as a double


% --- Executes during object creation, after setting all properties.
function cr2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cr2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cr3_Callback(hObject, eventdata, handles)
% hObject    handle to cr3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cr3 as text
%        str2double(get(hObject,'String')) returns contents of cr3 as a double


% --- Executes during object creation, after setting all properties.
function cr3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cr3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cr4_Callback(hObject, eventdata, handles)
% hObject    handle to cr4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cr4 as text
%        str2double(get(hObject,'String')) returns contents of cr4 as a double


% --- Executes during object creation, after setting all properties.
function cr4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cr4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in LineaVista.
function LineaVista_Callback(hObject, eventdata, handles)
% hObject    handle to LineaVista (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of LineaVista


% --- Executes on button press in Multipath.
function Multipath_Callback(hObject, eventdata, handles)
% hObject    handle to Multipath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Multipath


% --- Executes on button press in pushbutton38.
function pushbutton38_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton38 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

limpiar();
ParametrosRed=load('ParametrosRed.mat');
ParametrosRedTx=load('ParametrosRedTx.mat');
xtx1=ParametrosRedTx.ParametrosRedTx(1,1);
ytx1=ParametrosRedTx.ParametrosRedTx(1,2);
xtx2=ParametrosRedTx.ParametrosRedTx(1,3);
ytx2=ParametrosRedTx.ParametrosRedTx(1,4);
xtx3=ParametrosRedTx.ParametrosRedTx(1,5);
ytx3=ParametrosRedTx.ParametrosRedTx(1,6);
xtx4=ParametrosRedTx.ParametrosRedTx(1,7);
ytx4=ParametrosRedTx.ParametrosRedTx(1,8);



 if ParametrosRed.ParametrosRed(1,15)==0 && ParametrosRed.ParametrosRed(1,16)==0
     errordlg('¡ERROR! DATOS INVALIDOS O CAMPOS VACÍOS','Error');
 
elseif ParametrosRed.ParametrosRed(1,11)==0 && ParametrosRed.ParametrosRed(1,12)==0 && ParametrosRed.ParametrosRed(1,13)==0 && ParametrosRed.ParametrosRed(1,14)==0
    errordlg('¡ERROR! DATOS INVALIDOS O CAMPOS VACÍOS','Error');

 else
    [P_rec_A1,P_rec_A2,P_rec_A3,P_rec_A4,x,y]=NetworkLiFi(ParametrosRed.ParametrosRed(1,1),ParametrosRed.ParametrosRed(1,2),ParametrosRed.ParametrosRed(1,3),ParametrosRed.ParametrosRed(1,4),ParametrosRed.ParametrosRed(1,5),ParametrosRed.ParametrosRed(1,6),ParametrosRed.ParametrosRed(1,7),ParametrosRed.ParametrosRed(1,8),ParametrosRed.ParametrosRed(1,15),ParametrosRed.ParametrosRed(1,16),xtx1,ytx1,xtx2,ytx2,xtx3,ytx3,xtx4,ytx4);
 end


if ParametrosRed.ParametrosRed(1,11)==1 && ParametrosRed.ParametrosRed(1,12)==0 && ParametrosRed.ParametrosRed(1,13)==0 && ParametrosRed.ParametrosRed(1,14)==0
    P_rec_total_1ref=P_rec_A1;
elseif ParametrosRed.ParametrosRed(1,11)==0 && ParametrosRed.ParametrosRed(1,12)==1 && ParametrosRed.ParametrosRed(1,13)==0 && ParametrosRed.ParametrosRed(1,14)==0
    P_rec_total_1ref=P_rec_A2;
elseif ParametrosRed.ParametrosRed(1,11)==0 && ParametrosRed.ParametrosRed(1,12)==0 && ParametrosRed.ParametrosRed(1,13)==1 && ParametrosRed.ParametrosRed(1,14)==0
    P_rec_total_1ref=P_rec_A3;
elseif ParametrosRed.ParametrosRed(1,11)==0 && ParametrosRed.ParametrosRed(1,12)==0 && ParametrosRed.ParametrosRed(1,13)==0 && ParametrosRed.ParametrosRed(1,14)==1
    P_rec_total_1ref=P_rec_A4;
elseif ParametrosRed.ParametrosRed(1,11)==1 && ParametrosRed.ParametrosRed(1,12)==1 && ParametrosRed.ParametrosRed(1,13)==0 && ParametrosRed.ParametrosRed(1,14)==0
    P_rec_total_1ref=P_rec_A1+P_rec_A2;
elseif ParametrosRed.ParametrosRed(1,11)==1 && ParametrosRed.ParametrosRed(1,12)==0 && ParametrosRed.ParametrosRed(1,13)==1 && ParametrosRed.ParametrosRed(1,14)==0
    P_rec_total_1ref=P_rec_A1+P_rec_A3;
elseif ParametrosRed.ParametrosRed(1,11)==1 && ParametrosRed.ParametrosRed(1,12)==0 && ParametrosRed.ParametrosRed(1,13)==0 && ParametrosRed.ParametrosRed(1,14)==1
    P_rec_total_1ref=P_rec_A1+P_rec_A4;
elseif ParametrosRed.ParametrosRed(1,11)==0 && ParametrosRed.ParametrosRed(1,12)==1 && ParametrosRed.ParametrosRed(1,13)==1 && ParametrosRed.ParametrosRed(1,14)==0
    P_rec_total_1ref=P_rec_A2+P_rec_A3;
elseif ParametrosRed.ParametrosRed(1,11)==0 && ParametrosRed.ParametrosRed(1,12)==1 && ParametrosRed.ParametrosRed(1,13)==0 && ParametrosRed.ParametrosRed(1,14)==1
    P_rec_total_1ref=P_rec_A2+P_rec_A4;
elseif ParametrosRed.ParametrosRed(1,11)==0 && ParametrosRed.ParametrosRed(1,12)==0 && ParametrosRed.ParametrosRed(1,13)==1 && ParametrosRed.ParametrosRed(1,14)==1
    P_rec_total_1ref=P_rec_A3+P_rec_A4;
elseif ParametrosRed.ParametrosRed(1,11)==1 && ParametrosRed.ParametrosRed(1,12)==1 && ParametrosRed.ParametrosRed(1,13)==1 && ParametrosRed.ParametrosRed(1,14)==0
    P_rec_total_1ref=P_rec_A1+P_rec_A2+P_rec_A3;
elseif ParametrosRed.ParametrosRed(1,11)==1 && ParametrosRed.ParametrosRed(1,12)==1 && ParametrosRed.ParametrosRed(1,13)==0 && ParametrosRed.ParametrosRed(1,14)==1
    P_rec_total_1ref=P_rec_A1+P_rec_A2+P_rec_A4;
elseif ParametrosRed.ParametrosRed(1,11)==0 && ParametrosRed.ParametrosRed(1,12)==1 && ParametrosRed.ParametrosRed(1,13)==1 && ParametrosRed.ParametrosRed(1,14)==1
    P_rec_total_1ref=P_rec_A2+P_rec_A3+P_rec_A4;
elseif ParametrosRed.ParametrosRed(1,11)==1 && ParametrosRed.ParametrosRed(1,12)==0 && ParametrosRed.ParametrosRed(1,13)==1 && ParametrosRed.ParametrosRed(1,14)==1
    P_rec_total_1ref=P_rec_A1+P_rec_A3+P_rec_A4;
else
    P_rec_total_1ref=P_rec_A1+P_rec_A2 +P_rec_A3 +P_rec_A4;
end


P_rec_1ref_dBm=10*log10(P_rec_total_1ref);
axes(handles.axes2)
surf(x,y,P_rec_1ref_dBm')
grid on
title(handles.axes2,'Red Li-Fi','color','w')
xlabel(handles.axes2,'Ancho [m]','color','w')
ylabel(handles.axes2,'Largo [m]','color','w')
zlabel(handles.axes2,'Potencia [dBm]','color','w')
set(handles.axes2,'XColor','W');
set(handles.axes2,'YColor','w')
set(handles.axes2,'GridColor','k')
box on
rotate3d on
colorbar('color','w')


% --- Executes on button press in PosicionTX.
function PosicionTX_Callback(hObject, eventdata, handles)
% hObject    handle to PosicionTX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global boton1
boton1 = get(handles.PosicionTX,'value');
if boton1==1
    RedLiFiPtx;
    close(RedLiFi);
    
end

% --- Executes on button press in pushbutton40.
function pushbutton40_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton40 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global boton1
boton1 = get(handles.dco_ofdm,'value');
if boton1==1
     RedLiFiPtx;
    close(RedLiFi);
    
end


% --------------------------------------------------------------------
function uipushtool6_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
RedLiFiPtx;
close(RedLiFi);
