function varargout = Inicio(varargin)
% INICIO MATLAB code for Inicio.fig
%      INICIO, by itself, creates a new INICIO or raises the existing
%      singleton*.
%
%      H = INICIO returns the handle to a new INICIO or the handle to
%      the existing singleton*.
%
%      INICIO('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in INICIO.M with the given input arguments.
%
%      INICIO('Property','Value',...) creates a new INICIO or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Inicio_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Inicio_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Inicio

% Last Modified by GUIDE v2.5 01-Dec-2019 21:01:14

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Inicio_OpeningFcn, ...
                   'gui_OutputFcn',  @Inicio_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
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


% --- Executes just before Inicio is made visible.
function Inicio_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Inicio (see VARARGIN)

% Choose default command line output for Inicio
handles.output = hObject;


%Cargar una imagen/Inicio del programa
axes(handles.axes1);
%Guardamos el parth de la primera imagen
img1 = ('fondo.jpg');
%Leer la imagen para convertirla en matriz
imag=imread(img1);
%Mostrar la imagen en el axes
imshow(imag)
%Desabilitamos los ejes verticales y horizontales
axis off;

% %insercion del imagen bombillo
% b=imread('bombilla.png');
% axes(handles.axes2);
% axis off
% imshow(b); 

guidata(hObject, handles);

% UIWAIT makes Inicio wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Inicio_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in transmisor.
function transmisor_Callback(hObject, eventdata, handles)
% hObject    handle to transmisor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global boton1
boton1 = get(handles.transmisor,'value');
if boton1==1
    sistemaLi_Fi;
    close(Inicio);
    
end


% --- Executes on button press in red.
function red_Callback(hObject, eventdata, handles)
% hObject    handle to red (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global boton1
boton1 = get(handles.red,'value');
if boton1==1
    RedLiFiPtx;
    close(Inicio);
    
end

% --- Executes on button press in receptor.
function receptor_Callback(hObject, eventdata, handles)
% hObject    handle to receptor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
