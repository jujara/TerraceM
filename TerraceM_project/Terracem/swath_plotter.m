function varargout = swath_plotter(varargin)%varargin
% SWATH_PLOTTER MATLAB code for swath_plotter.fig
%      SWATH_PLOTTER, by itself, creates a new SWATH_PLOTTER or raises the existing
%      singleton*.
%
%      H = SWATH_PLOTTER returns the handle to a new SWATH_PLOTTER or the handle to
%      the existing singleton*.
%
%      SWATH_PLOTTER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SWATH_PLOTTER.M with the given input arguments.
%
%      SWATH_PLOTTER('Property','Value',...) creates a new SWATH_PLOTTER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before swath_plotter_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to swath_plotter_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help swath_plotter

% Last Modified by GUIDE v2.5 12-Mar-2015 18:20:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @swath_plotter_OpeningFcn, ...
                   'gui_OutputFcn',  @swath_plotter_OutputFcn, ...
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


% --- Executes just before swath_plotter is made visible.
function swath_plotter_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to swath_plotter (see VARARGIN)

% Choose default command line output for swath_plotter
handles.output = hObject;

%TerraceM link 
T=load('terracem_path.mat');
handles.station=T.station;
handles.profsel=T.profsel;
%set directories
dirstations=T.stationsdir;
maindir=T.terracem_dir;

%detect selected profile and station
station=handles.station;
profnum=str2double(handles.profsel); 

cd(dirstations);
cd(station);

%load swaths file
g=load(sprintf('%s_swaths.mat',station)); %load swaths.mat database
%select inputs
swat=g.(sprintf('swat%u',profnum));
ZI=g.(sprintf('ZI%u',profnum));
Dprof=g.(sprintf('Dprof%u',profnum));
%KsY=g.(sprintf('KsY%u',profnum));


handles.swat=swat;  %.txt export output
handles.xsw=swat(:,1); %distance along profile
handles.misw=swat(:,2); %min values%
handles.mesw=swat(:,3); %mean values
handles.masw=swat(:,4); %max values

handles.ZI=ZI; %mean values
handles.Dprof=Dprof; %max values

%handles.Ks=Ks;
%handles.KsX=KsX;
%handles.KsY=KsY;



%Initial plot
hold on
plot(swat(:,1),swat(:,4),'k')
plot(swat(:,1),swat(:,3),'color',[0.3 0.3 0.3])
plot(swat(:,1),swat(:,2),'color',[0.6 0.6 0.6])
%hleg1 = legend('Max topography','Mean topography','Min topography','Location','SouthEast');
xlim([min(swat(:,1)) max(swat(:,1))]); 
xlabel('Distance along swath [m]');
ylabel('Elevation [m]'); 

tit=sprintf('Station %s, Profile %u',station,profnum);
title(tit);

%legend handles
%handles.hleg1=hleg1;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes swath_plotter wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = swath_plotter_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles) %flip profile
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1

flipi=get(hObject,'Value'); %flipi is refered to flip profile
handles.flipi=flipi; 

ax_plot=handles.axes1;
axes(ax_plot);
if handles.flipi==1 %if checkbox is turned on

set(gca,'XDir','reverse')

else
    
set(gca,'XDir','normal')
end

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.output = hObject;
ax_plot=handles.axes1;
axes(ax_plot);

station=handles.station;
profnum=str2double(handles.profsel); %selected profile

rect=[-1,-2.5,30,18]; %[xmin ymin width height] 
%set(gcf,'PaperType','A4');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf,'paperorientation','portrait')
set(gcf,'PaperUnits','centimeters'); 
set(gcf, 'PaperSize', [27.3 13]); %(width height)
set(gcf,'Paperposition',rect);
fout = sprintf('%s_%u_swath.pdf',station,profnum);
saveas(gcf,fout,'pdf')

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.output = hObject;


station=handles.station;
profnum=str2double(handles.profsel); %selected profile

swat=handles.swat;
%Export swath ASCII
fout = sprintf('%s_%u_swath.txt',station,profnum);
save(fout,'swat','-ASCII')

% --- Executes when selected object is changed in uipanel1.
function uipanel1_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel1 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

handles = guidata(hObject); %define handles
a1=handles.axes1;
axes(a1);

station=handles.station;
profnum=str2double(handles.profsel); %selected profile
%hleg1=handles.hleg1;

x=handles.xsw;
mi=handles.misw;
me=handles.mesw;
ma=handles.masw;
ZI=handles.ZI;
ddg=handles.Dprof;
%Ks=handles.Ks;
%KsX=handles.KsX;
%KsY=handles.KsY;

newButton=get(eventdata.NewValue,'tag');

switch newButton
    
case 'R1' %line swath plot
cla(a1) %clear plot

hold on
plot(x,ma,'k')
plot(x,me,'color',[0.3 0.3 0.3])
plot(x,mi,'color',[0.6 0.6 0.6])
hleg2 = legend('Max topography','Mean topography','Min topography','Location','SouthEast');
xlim([min(x) max(x)]);  
xlabel('Distance along swath [m]');
ylabel('Elevation [m]');  

tit=sprintf('Station %s, Profile %u',station,profnum);
title(tit);

case 'R2' %grey swath plot 
cla(a1) %clear plot
%cla(hleg2)
  
hold on
ma(isnan(ma)) = 0; %convert nan values to zero
mi(isnan(mi)) = 0; %convert nan values to zero

U7=area(x,ma);
U6=area(x,mi);
set(U7,'FaceColor',[0.6 0.6 0.6],'edgeColor',[0.6 0.6 0.6]);
set(U6,'FaceColor',[1 1 1],'edgeColor',[0.6 0.6 0.6]);
set(gca,'Layer','top')
plot(x,me,'-k','linewidth',1.4)
xlim([min(x) max(x)]); 
ylim([min(mi)-(min(mi)/8) max(ma)+(max(ma)/8)]); 

xlabel('Distance along swath [m]');
ylabel('Elevation [m]');    
    
tit=sprintf('Station %s, Profile %u',station,profnum);
title(tit);    
    
case 'R3' %density swath plot
cla(a1) %clear plot

ZI=ZI';
total_line_dist=max(x);
sma=max(ZI); 
sme=mean(ZI); 
smi=min(ZI);

sma(isnan(sma)) = 0; %convert nan values to zero
smi(isnan(smi)) = 0; %convert nan values to zero

sx=linspace(0,total_line_dist,numel(sma)); %distance along swath
%swat=[sx'./1e3 smi' sme' sma'];

% compute pdf
[~, nxZ]=size(ZI);
hx=linspace(min(ZI(:)),max(ZI(:)),200);
for i=1:nxZ
    Ks(:,i)=(ksdensity(ZI(:,i),hx,'function','pdf')); 
end


Ks(Ks>=1)=1;%remove outliers from grid

[KsX,KsY]=meshgrid(sx,hx);

%display in logarithmic color scale


hold on
pcolor(KsX,KsY,Ks)
colormap(flipud(hot))%hot
hcb=colorbar;

shading flat


ma(isnan(ma)) = 0; %convert nan values to zero
mi(isnan(mi)) = 0; %convert nan values to zero

plot(sx,ma,'color',[0.6 0.6 0.6])
plot(sx,mi,'color',[0.6 0.6 0.6])
plot(sx,me,'-k','linewidth',1.4)

xlim([min(sx) max(sx)]);  
xlabel('Distance along swath [m]');
ylabel('Elevation [m]');  

tit=sprintf('Station %s, Profile %u',station,profnum);
title(tit);  
end    


% --- Executes on button press in loadsh.
function loadsh_Callback(hObject, eventdata, handles)
% hObject    handle to loadsh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of loadsh

%path to stations
%dirstations=handles.stationsdir; 
%cd(dirstations);

handles = guidata(hObject); %define handles
a1=handles.axes1;

oldFolder = pwd; %store present folder location
axes(a1);

cd('..') %go up one level to stations folder
%load table
file=('shorelines.txt');
if exist(file,'file')==1 || exist(file,'file')==2
    
%load shorelines table
fid=fopen(file,'r');
i=1;
while ~feof(fid)
    A=textscan(fid,'%s %u %u %f %f %f %f %f %s %f/','delimiter',' ','HeaderLines',1);
    
    for j=1:10 %number of columns (default from shorelines.txt)
        data(i,j) = A(1,j);
    end
    i=i+1;
end
fclose(fid);
cd(oldFolder)

%extract data from table
nim=numel(data(:,10)); 
u=1;
for k=1:nim
md1(u,1)=data{k,4}; %x
md1(u,2)=data{k,5}; %Y
md1(u,3)=data{k,7}; %z
md1(u,4)=data{k,8}; %ze
md1(u,5)=data{k,3}; %shore
md1(u,6)=data{k,2}; %profnum
md1(u,7)=data{k,6}; %distance along profile
vars{u,1}=data{k,1}{1,1}; %station names
u=1+u;
end


name=(vars);

%shorelines clasification of names
name2=cellfun(@double,name,'uniformoutput',0);
code_shore=cell2mat(name2);
tab(:,1)=code_shore(:,1); %code_name  station
tab(:,2)=code_shore(:,2); %code_name  station
%tab(:,3)=code_shore(:,3); %code_name  station
%tab(:,4)=code_shore(:,4); %code_name  station
tab(:,5)=md1(:,3); %shoreline elevation
tab(:,6)=md1(:,4); %shoreline error
tab(:,7)=md1(:,5); %terrace shore
tab(:,8)=md1(:,6); %profile
tab(:,9)=md1(:,7); %distance along profile

statsel=cellstr(handles.station);
name3=(cellfun(@double,statsel,'uniformoutput',0));

name3=cell2mat(name3);
profilese=str2num(handles.profsel);


nwename(:,1)=code_shore(:,1);
nwename(:,2)=code_shore(:,2);
nwename(:,3)=code_shore(:,3);
nwename(:,4)=code_shore(:,4);
t=1;
xc=1;
for j=1:numel(nwename(:,1))
hold on   
if nwename(j,1:4)==name3(1,1:4) 

xc=xc+1;
if tab(j,8)==profilese
count(j,1)=1;
%error bar    
errorbar(tab(j,9),tab(j,5),tab(j,6),'.k') 

plot(tab(j,9),tab(j,5),'o','markerfacecolor',[1 0 1])
 
end

xtop=get(gca,'xlim');
ytop=get(gca,'ylim');
%color coding by disntace (not required in this version)
XB(:,1)=md1(:,1);
XB(:,2)=md1(:,2);
XA(:,1)=tab(j,1);
XA(:,2)=tab(j,2);


Xd = [XA;XB(j,1:2)];
d(t) = pdist(Xd,'euclidean');       
    
t=t+1;
end
end


profnum=str2double(handles.profsel);
%load line from structure
statsel=(handles.station);
det=sprintf('%s_fits.mat',statsel);
load(det);

%loop for shorenum
%shorenum=1;
try
count(count==0)=[];

for shorenum=1:numel(count)
if isstruct(Fits.(sprintf('%s_%u_%u_fit',statsel,profnum,shorenum)))==1    
    
linesplot=Fits.(sprintf('%s_%u_%u_fit',statsel,profnum,shorenum));

Xlinesplot1=linesplot.children.children(4,1).properties;
Xlinesplot2=linesplot.children.children(8,1).properties; %fields exploration for fits 6 8 9 10 

plot(Xlinesplot1.XData,Xlinesplot1.YData,'--r')%load from fig-structure
plot(Xlinesplot2.XData,Xlinesplot2.YData,'--r')%load from fig-structure

xlim(xtop); ylim(ytop)
      
end  
end

catch %detect error when loading shoreline angles

end
tab(:,10)=d';

end
