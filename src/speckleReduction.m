function varargout = speckleReduction(varargin)
% SPECKLEREDUCTION M-file for speckleReduction.fig
%      SPECKLEREDUCTION, by itself, creates a new SPECKLEREDUCTION or raises the existing
%      singleton*.
%
%      H = SPECKLEREDUCTION returns the handle to a new SPECKLEREDUCTION or the handle to
%      the existing singleton*.
%
%      SPECKLEREDUCTION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SPECKLEREDUCTION.M with the given input arguments.
%
%      SPECKLEREDUCTION('Property','Value',...) creates a new SPECKLEREDUCTION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before speckleReduction_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to speckleReduction_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help speckleReduction

% Last Modified by GUIDE v2.5 02-Dec-2010 21:20:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @speckleReduction_OpeningFcn, ...
                   'gui_OutputFcn',  @speckleReduction_OutputFcn, ...
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


% --- Executes just before speckleReduction is made visible.
function speckleReduction_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to speckleReduction (see VARARGIN)

% Choose default command line output for speckleReduction
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes speckleReduction wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = speckleReduction_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
FileName = get(handles.edit1,'String');
disp('Submit button pressed.');
if(exist(FileName)==2)
I = imread(FileName);
TypeList = get(handles.popupmenu1,'String');
IndexSelected = get(handles.popupmenu1,'Value');
TypeSelected = TypeList{IndexSelected}; 
switch IndexSelected
    case 2
        figure,imshow(I),title(FileName);
        O = srad(I);
        figure,imshow(O),title('Denoised image - SRAD filter');
    case 3
        figure,imshow(I),title(FileName);
        O = lee(I);
        figure,imshow(O),title('Denoised image - SRAD filter');
    case 4
        figure,imshow(I),title(FileName);
        O = frost(I);
        figure,imshow(O),title('Denoised image - SRAD filter');
    case 5
        figure,imshow(I),title(FileName);
        O = pde_noise_reduction(I);
        figure,imshow(O),title('Denoised image - SRAD filter');
    otherwise
        disp('Please select a filter technique.');
end
else
    disp('File not found');
end

% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


function compare_with_imnoise()
FileName = get(handles.edit1,'String');
disp('Submit button pressed.');
if(exist(FileName)==2)
I1 = imread(FileName);
disp(FileName);
figure,imshow(I1),title(FileName);
[x y z]=size(I1);
Q = 255;
j=1;
    for j=1:10
        i=j*0.002;
        disp(i);
        I2=imnoise(I1,'speckle',i);
        %figure,imshow(I),title(i);
        I1=pde_noise_reduction(I1);
        MSE= sum(sum((I2-I1) .^ 2)) / (x * y)
        %disp(MSE);
        rmse=sqrt(MSE);
        SNR=10*log10(sum(sum((I1).^2)))./(sum(sum((I2-I1) .^ 2)));
        PSNR = 10*log10(Q*Q./MSE)
        %disp(PSNR);
    end

end

% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


function compare_algos()
disp('Submit button pressed.');
FileName = 'gray-3.jpg'
if(exist(FileName)==2)
I2 = imread(FileName);
I1 = srad(I2);
ISNR_srad=10*log10(sum(sum((I1).^2)))./(sum(sum((I2-I1) .^ 2)))
I2 = lee(I1);
ISNR_lee=10*log10(sum(sum((I1).^2)))./(sum(sum((I2-I1) .^ 2)))
I2 = frost(I1);
ISNR_frost=10*log10(sum(sum((I1).^2)))./(sum(sum((I2-I1) .^ 2)))
I2 = pde_noise_reduction(I1);
ISNR_pde=10*log10(sum(sum((I1).^2)))./(sum(sum((I2-I1) .^ 2)))
else
    disp('file not found');
end


function S=srad(I)
tic
%I=imread('1.png');
T=3;
[x y]=size(I);
I=double(I);
Ic=double(I);
delta_t = 0.08;
t=1;
eps=0.00000000001;
for t=1:T     
    qt=exp(-t*.2);
    [Ix,Iy] = gradient(Ic);    
    di=sqrt(Ix.^2+Iy.^2);
    di2=del2(Ic);
    T1=0.5*((di./(Ic+eps)).^2);
    T2=0.0625*((di2./(Ic+eps)).^2);
    T3=(1+(0.25*(di2./(Ic+eps)))).^2;
    T=sqrt((T1-T2)./(T3+eps));
    dd=(T.^2-qt.^2)./((qt.^2*(1+qt.^2)+eps));
    cq=1./(1+dd);
    [D1,D2]=gradient(cq.*Ix);
    [D3,D4]=gradient(cq.*Iy);
    D=D1+D4;    
    Ic=real(Ic+delta_t .*D);   
end
toc
S=uint8(Ic);
%figure,imshow(S),title('SRAD - Image after Speckle Reduction');
imwrite(S,'srad.jpg');

function [ft]=frost(I)
tic
%I=imread('3.jpg');
[x y z]=size(I);
I=double(I);
K=1;
N=I;
for i=1:x
    for j=1:y                              
        if (i>1 & i<x & j>1 & j<y)
            mat(1)=I(i-1,j);
            mat(2)=I(i+1,j);
            mat(3)=I(i,j-1);
            mat(4)=I(i,j+1);
            d(1)=sqrt((i-(i-1))^2);
            d(2)=sqrt((i-(i+1))^2);
            d(3)=sqrt((j-(j-1))^2);
            d(4)=sqrt((j-(j+1))^2);
            mn=mean(mean(mat));
            c=mat-mn;
            c2=c.^2;
            c3=c/(c2+.0000001);
            Cs=0.25*sum(sum(c3));
            m(1)=exp(-K*Cs*d(1));
            m(2)=exp(-K*Cs*d(2));
            m(3)=exp(-K*Cs*d(3));
            m(4)=exp(-K*Cs*d(4));
            ms=sum(sum(m));
            mp=m/ms;
            N(i,j)=sum(sum(mp.*mat));                    
        end
     end
end
toc;
ft=uint8(N);
%figure,imshow(ft),title('Frost - Image after Speckle Reduction');
imwrite(ft,'frost2.jpg');

function [le]=lee(I)
%I=imread('3.jpg');
[x y z]=size(I);
I=double(I);
N=zeros(x,y,z);
for i=1:x
    
    for j=1:y
        % Checking first and last pixel of first row% 
        if (i==1 & j==1)
            mat(1)=0;
            mat(2)=0;
            mat(3)=0;
            mat(4)=0;
            mat(5)=I(i,j);
            mat(6)=I(i,j+1);
            mat(7)=0;
            mat(8)=I(i+1,j);
            mat(9)=I(i+1,j+1);
        end
        
        if (i==1 & j==y)
            mat(1)=0;
            mat(2)=0;
            mat(3)=0;
            mat(4)=I(i,j-1);
            mat(5)=I(i,j);
            mat(6)=0;
            mat(7)=I(i+1,j-1);
            mat(8)=I(i+1,j);
            mat(9)=0;
        end
        
        % Checking first and last pixel of last row% 
        if (i==x & j==1)
            mat(1)=0;
            mat(2)=I(i-1,j);
            mat(3)=I(i-1,j+1);
            mat(4)=0;
            mat(5)=I(i,j);
            mat(6)=I(i,j+1);
            mat(7)=0;
            mat(8)=0;
            mat(9)=0; 
        end
        
        if (i==x & j==y)
            mat(1)=I(i-1,j-1);
            mat(2)=I(i-1,j);
            mat(3)=0;
            mat(4)=I(i,j-1);
            mat(5)=I(i,j);
            mat(6)=0;
            mat(7)=0;
            mat(8)=0;
            mat(9)=0;
        end
        % Checking rest of the image%        
        if (i>1 & i<x & j>1 & j<y)
            mat(1)=I(i-1,j-1);
            mat(2)=I(i-1,j);
            mat(3)=I(i-1,j+1);
            mat(4)=I(i,j-1);
            mat(5)=I(i,j);
            mat(6)=I(i,j+1);
            mat(7)=I(i+1,j-1);
            mat(8)=I(i+1,j);
            mat(9)=I(i+1,j+1);
        end
        y1=I(i,j);
        ybar=mean(mean(mat));
        if ybar~=0
            ystad=std2(mat);
            ENL=(ybar/ystad)^2;
            sx2=((ENL*(ystad)^2)-(ybar)^2)/(ENL+1);
            xcap=ybar+(sx2*(y1-ybar)/(sx2+(ybar^2/ENL)));
            N(i,j)=xcap;
        else
            N(i,j)=y1;
        end
    end
end
le=uint8(N);  
%figure,imshow(le),title('Lee - Image after Speckle Reduction');
imwrite(le,'lee2.jpg');

function P=pde_noise_reduction(I)
%I=imread('1.jpg');
[x y]=size(I);
I=double(I);
I2=I;

T=7;
for t=1:T
    [Ix,Iy]=gradient(I);
    %c = 1./(1.+(Ix.^2+Iy.^2));
    c=exp((-1).*(sqrt(Ix.^2+Iy.^2))./7000);
    [div1,divt1]=gradient(c.*Ix);
    [divt2,div2]=gradient(c.*Iy);
    div=div1+div2;
    I1=I+(0.20).*div;
    I=I1;
end;
I=I2-I1;
P=uint8(I1);
imwrite(P,'pde.jpg');
%figure,imshow(P),title('PDE - Image after Speckle Reduction');


function wiener(f)
%f=imread('3.jpg');
f=im2double(f);
[r c]=size(f);
h=fspecial('gaussian',[r c],1.5);
g=imfilter(f,h,'circular');

a=0.1;
b=0.01;
n=a+b*randn(r,c);
%g=g+n;

F=fft2(f);
G=fft2(g);
H=psf2otf(h);
N=fft2(n);
H2=conj(H).*H;
Sn=conj(N).*N;
Sf=conj(F).*F;

snsf=Sn./Sf;
huv=H2./(H2+snsf);
tuv=(1./(H+eps)).*huv;
RA=tuv.*G;
ra=real(ifft2(RA));

smsn=sum(sum(Sn));
smsf=sum(sum(Sf));
sr=smsn/smsf;
huv=H2./(H2+sr);
tuv=(1./(H+eps)).*huv;
RC=tuv.*G;
rc=real(ifft2(RC));

figure,imshow(f),title('Original image');
figure,imshow(g),title('Noisy image');
figure,imshow(ra),title('Restored image - Using Autocorrelation');
figure,imshow(rc),title('Restored image - Using Constant Ratio');
imwrite(ra,'wiener3.jpg');
imwrite(rc,'wiener4.jpg');
