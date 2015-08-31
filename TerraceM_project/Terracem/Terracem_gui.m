function varargout = Terracem_gui(varargin)
% TERRACEM_GUI MATLAB code for Terracem_gui.fig
%
%%%%%%%%%%%%%%%%%%%%
%to run TERRACEM type this line in the command window:
%>> Terracem_gui
%%%%%%%%%%%%%%%%%%%%
%     
% Begin initialization code - DO NOT EDIT (GUIDE INITIALIZATION)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Terracem_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @Terracem_gui_OutputFcn, ...
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

%% INTERFACE FUNCTIONS %%

%% Opening function
function Terracem_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% Choose default command line output for Terracem_gui
handles.output = hObject;
set(handles.sureornot,'Value',0)
set(handles.t2,'String','35000');
set(handles.ti,'String','50');
set(handles.t1,'String','0.5');
set(handles.difu,'String','0.8');
set(handles.slope,'String','50');
set(handles.Search_ratio,'String','0.5');
set(handles.edit14,'String','70');
set(handles.edit15,'String','5');
set(handles.resoluti,'String','2');
set(handles.azimuti,'String','300');
set(handles.svpdf,'Value',0);
%set handles for image control in listbox2, all desactivated 
handles.demsel=0;
handles.shadsel=0; 
handles.imsel=0; 
% condition if to check existence of .path file
path = 'terracem_path.mat'; %path file
cdx=pwd; %local directory
%check if path file exist
if exist(fullfile(cdx, 'terracem_path.mat'), 'file')==2     
    load(path) %load path file       
    try
    if pwd==terracem_dir 
        handles.maindir=(terracem_dir);
        handles.stationsdir=stationsdir;
        handles.zone=zoneutm;     
        %set directories readed from terracem_path.mat
        set(handles.fullpath1,'String',handles.maindir) %print terracem directory path in fullpath1 bar
        set(handles.fullpath2,'String',handles.stationsdir) %print stations directory in fullpath2 bar
        set(handles.utmzonen,'String',handles.zone)    
        %prompt message
        Commands=('Welcome to TerraceM. The directories have been automatically loaded, push "Read Inputs" to proceed');
        T = evalc('Commands');
        Tcell = regexp(T, '\n', 'split');
        set(handles.text38,'String',Tcell);    
    else %If path file differs from directory achitecture    
        %prompt message    
        msgbox('Welcome. seems to be an error in the directories, please set manually the directories and UTM zone and press "Read Inputs"');         
    end     
    catch
    msgbox('the path file is corrupt, please introduce again TerraceM inputs folders')
    delete('terracem_path.mat')  
    end
else
    %command window prompt 
    Commands=('Welcome. Set the imput directories for the TerraceM scripts and your stations');
    T = evalc('Commands');
    Tcell = regexp(T, '\n', 'split');
    set(handles.text38,'String',Tcell);
end  

%load logos
lgt=handles.logotype;
axes (lgt);
[cdata,map] = imread('logo.tif'); 
imshow(cdata,map)

ax1=handles.map;
axes (ax1);
imshow(cdata,map)

handles.edit15=5;
handles.edit14=70;
% Update handles structure
guidata(hObject, handles);

%% outputs function DO NOT EDIT
function varargout = Terracem_gui_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

%% Small logo TerraceM
function logotype_CreateFcn(hObject, eventdata, handles)
handles = guidata(hObject);
handles.logotype=axis;

guidata(hObject, handles);

%% main plot area
function map_CreateFcn(hObject, eventdata, handles)
handles = guidata(hObject);
handles.map=axis;

guidata(hObject, handles);



%% About information display
function about_Callback(hObject, eventdata, handles)

    Commands=('About the authors. Please send us a feedback (jara@geo.uni-potsdam.de)');
    T = evalc('Commands');
    Tcell = regexp(T, '\n', 'split');
    set(handles.text38,'String',Tcell);
[cdata,map] = imread('About.tif'); 
hold on
figure
imshow(cdata,map)

%% clean project button
function cleanstat_Callback(hObject, eventdata, handles)

Commands=('Debbugging deletes all path.m files in directories. We suggest debugg your project before moving it into a different directory ');
    T = evalc('Commands');
    Tcell = regexp(T, '\n', 'split');
    set(handles.text38,'String',Tcell);

station=handles.statsel; %Selected station directory
dirstations=handles.stationsdir; %Stations directory
maindir=handles.maindir;

cd(dirstations)
cd(station)
cdx1=pwd;

if exist(fullfile(cdx1, 'terracem_path.mat'), 'file')==2 
    delete('terracem_path.mat')
end

cd(maindir)
cdx2=pwd;
if exist(fullfile(cdx2, 'terracem_path.mat'), 'file')==2 
    delete ('terracem_path.mat')
end

%% Help buttons 
function pushbutton23_Callback(hObject, eventdata, handles)
msgbox('Swath profile extraction. This is the second step during first-time use of TerraceM. Please select a Station and press the "SWATH PROFILE EXTRACTION" button to start.')

function pushbutton43_Callback(hObject, eventdata, handles)
msgbox('This function needs to be executed each time the profile polygons (Shapefiles) are edited.')

%%
function pushbutton20_Callback(hObject, eventdata, handles)
msgbox('This is the first step to start using TerraceM. Please define the directories.')
guidata(hObject, handles);

%%
function pushbutton21_Callback(hObject, eventdata, handles)
msgbox('Profile selection list. Used to select the profile for analysis and refresh the map viewer. Selected profile will appear in red')

%%
function pushbutton22_Callback(hObject, eventdata, handles)
msgbox('Select the Station to extract the swath profiles, to display the profile list and also to refresh the mapviewer.')

%%
function pushbutton24_Callback(hObject, eventdata, handles)
msgbox('Parameters box. Used to visualize and modify the parameters for the analysis functions. ')

%%
function pushbutton25_Callback(hObject, eventdata, handles)
msgbox('Staircase analysis. This routine extimates the shoreline angle by back-projecting the slope of the abrasion platform to the paleocliff.')

%%
function pushbutton26_Callback(hObject, eventdata, handles)
msgbox('This routine is created to analyze rough coasts, calculating the maximum elevation of stacks. Modify the peak search parameter according to the DEM resolution (peak search ratio ~ DEM pixel size).')

%%
function pushbutton27_Callback(hObject, eventdata, handles)
msgbox('This routine calculates the geomorphic age of a terrace (Kt) and the shoreline angle.')

%%
function pushbutton29_Callback(hObject, eventdata, handles)
msgbox('Select the output format for shoreline angles')

%%
function pushbutton41_Callback(hObject, eventdata, handles)
msgbox('This routine uses a defined cliff slope.')

%%
function pushbutton35_Callback(hObject, eventdata, handles)
msgbox('Shoreline-angles and Kt values can be reviewed and edited using the table buttons')

%%
function pushbutton34_Callback(hObject, eventdata, handles)
msgbox('Post processing functions')

%%
function pushbutton33_Callback(hObject, eventdata, handles)
msgbox('Export your shoreline-angle results to PDF')

%% Command window Information display (static text)
function text38_CreateFcn(hObject, eventdata, handles)
handles = guidata(hObject);    
handles.text38=str2double(get(hObject,'String')); %returns contents of edit14 as a double
guidata(hObject, handles);

%% Directories interface 

%% TerraceM directory
function maindir_Callback(hObject, eventdata, handles)

handles = guidata(hObject);
dir1=char(uigetdir(pwd,'Set the path to Terracem script files'))
handles.maindir=char(dir1); %define handles for the path to terracem directory
set(handles.fullpath1,'String',dir1) %print terracem directory path in fullpath1 bar
% Update handles structure
guidata(hObject, handles);

%% Set station directory button
function stationsdir_Callback(hObject, eventdata, handles)
handles = guidata(hObject);
dir2=char(uigetdir(pwd,'Set folder with stations files')) %navigate to stations folder
handles.stationsdir=char(dir2); %define handles for stations directory
set(handles.fullpath2,'String',dir2) %print stations directory in fullpath2 bar
% Save handles 
guidata(hObject, handles);

%% Set path 1
function path2dir1_Callback(hObject, eventdata, handles)

%% Set path 1
function path2dir1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Set path 2
function path2dir2_Callback(hObject, eventdata, handles)

%% Set path 2
function path2dir2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Main directory
function maindir_KeyPressFcn(hObject, eventdata, handles)

%% Display TerraceM directories in textbox
function fullpath1_Callback(hObject, eventdata, handles)
guidata(hObject, handles);

%% set TerraceM directories textbox properties
function fullpath1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Display station directory
function fullpath2_Callback(hObject, eventdata, handles)
guidata(hObject, handles);

%% Display stations directory
function fullpath2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% DEM display interface %%%%%%%%%%%%%%%%

%% DEM resolution textbox
function resoluti_Callback(hObject, eventdata, handles)
handles = guidata(hObject);
handles.resoluti=str2double(get(hObject,'String')); %returns contents of edit14 as a double
guidata(hObject, handles);

%% DEM resolution textbox properties
function resoluti_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% DEM Shaded textbox
function azimuti_Callback(hObject, eventdata, handles)
handles = guidata(hObject);
handles.azimuti=str2double(get(hObject,'String')); %returns contents of edit14 as a double
guidata(hObject, handles);

%% DEM shaded azimut textbox properties
function azimuti_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% plot options button group properties
function plot_options_CreateFcn(hObject, eventdata, handles)
handles = guidata(hObject);
handles.radio = hObject;
guidata(hObject, handles);

%% plot options button group interface
function plot_options_SelectionChangeFcn(hObject, eventdata, handles)
ax1=axis;
newButton=get(eventdata.NewValue,'tag');
maindir=handles.maindir;
stationsdir=handles.stationsdir;
station=handles.statsel;
UTM_zone=get(handles.utmzonen,'String');

if exist('handles.shpoint')==1
    D=handles.shpoint;
else
    D=[];
end

SF=handles.SF;
nim=numel(SF);
cd(char(maindir));
switch newButton
    
case 'radiobutton1'   
%set handles for image control in listbox2 
handles.demsel=0; %denote non dem activation
handles.shadsel=0; %GoogleMaps shaded desactivated
handles.imsel=0; %GoogleMaps image desactivated

%prompt imformation
Commands=('The display interface may be clean each time the display is changed (e.g. satellite, shaded or DEM)');
T = evalc('Commands');
Tcell = regexp(T, '\n', 'split');
set(handles.text38,'String',Tcell);
cla

%plot profiles again
for tk=1:nim   
    obj.h=plot(SF(tk,1).Lon, SF(tk,1).Lat,'-k'); %plot shapefiles in 2D latlong system
    hold on
end
axis (ax1)
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'Position', [4 8.7 111 44]);
refreshdata(obj.h)

%plot shorelines again
if exist('handles.shpoint')==1
obj.p=plot(D(:,1),D(:,2),'o','MarkerEdgeColor','k','MarkerFaceColor','b');
end
axis (ax1)
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'Position', [4 8.7 111 44]);
axis equal
drawnow

case 'radiobutton2' %Load shaded relief GoogleMaps
%set handles for image control in listbox2 
handles.demsel=0; %denote non dem activation
handles.shadsel=1; %GoogleMaps shaded activated
handles.imsel=0; %GoogleMaps image desactivated
%set googlemaps style
param='&style=feature:landscape.natural.terrain&style=element:geometry.fill|color:0xffffff|saturation:-100|lightness:50|visibility:off|gamma:0.14';
%call plot googlemaps function
h1=plot_google_map('MapType','terrain','style',param);
%promt information on display
Commands=('GoogleMaps shaded relief image loaded');
T = evalc('Commands');
Tcell = regexp(T, '\n', 'split');
set(handles.text38,'String',Tcell);       
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'Position', [4 8.7 111 44]);
axis equal

case 'radiobutton3' %Load shaded rsatellite image GoogleMaps
%set handles for image control in listbox2 
handles.demsel=0; %denote non dem activation
handles.shadsel=0; %GoogleMaps shaded desactivated
handles.imsel=1; %GoogleMaps image activated    
%call plot googlemaps function
plot_google_map('MapType','satellite','style',0);
%promt message
Commands=('GoogleMaps satellite image loaded');
T = evalc('Commands');
Tcell = regexp(T, '\n', 'split');
set(handles.text38,'String',Tcell);       
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'Position', [4 8.7 111 44]);
axis equal

case 'radiobutton4' %load dem
%set handles for image control in listbox2 
handles.demsel=1; %dem display activation reflected in handles
handles.shadsel=0; %GoogleMaps shaded desactivated
handles.imsel=0; %GoogleMaps image desactivated
%change directory to DEM_XXX.mat location
cd(char(stationsdir)); 
cd(station)
%call to axis of main plot
axis (ax1)
% prompt information
Commands=('The DEM displays pixelated the first time is loaded, change the resolution and clean plot to adjust the display');
T = evalc('Commands');
Tcell = regexp(T, '\n', 'split');
set(handles.text38,'String',Tcell);
%DEM loading procedure
%load &/or combine structures
DEMS=(sprintf('DEM_%s.mat',station));

if exist((DEMS),'file')==2    
load(DEMS)
X=DEMS.X;
Y=DEMS.Y;
Z=DEMS.Z;
xa=X'; 
xa(:,2)=Y(1,1);
ka=numel(xa(:,1));

    for f=1:ka; 
    xut(f,:)=(UTM_zone); %generate  UTM zone matrix
    end

ya=Y;
ya(:,2)=X(1,1); %generate sinthetic utm coordinates for LatLon transformation
ky=numel(ya(:,1));

    for k=1:ky; 
    yut(k,:)=(UTM_zone); %generate  UTM zone matrix
    end

% Convert bounding coordinates of raster area to LatLong
[xLa xLo]=utm2deg(xa(:,1),xa(:,2),xut);
[yLa yLo]=utm2deg(ya(:,2),ya(:,1),yut);

%plot using Beaducel (2014) DEM shading routine
Zf=flipud(DEMS.Z);
hold on
%manual inputs
N=handles.resoluti; %resampling 1/N
Az=handles.azimuti; %Azimut
[H I]=dem(xLo,yLa,Zf,'Azimuth',Az,'Contrast',0.5,'Lake',...
    'Decim',N,'Watermark',0,'LandColor',landcolor,'AxisEqual','off','NaNColor',[1 1 1],'ZCut',0);
%change directory to stationsdir
cd(char(stationsdir));

%plot profiles again
for tk=1:nim   
obj.h=plot(SF(tk,1).Lon, SF(tk,1).Lat,'-k'); %plot shapefiles in 2D latlong system
hold on
end

set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'Position', [4 8.7 111 44]);
refreshdata(obj.h)

%plot shorelines again
if exist('handles.shpoint')==1
obj.p=plot(D(:,1),D(:,2),'o','MarkerEdgeColor','k','MarkerFaceColor','b');%plot shorelines
refreshdata(obj.p)
end
refreshdata(H)
else    
msgbox('Error, please extract swath profiles again')
end    
end

guidata(hObject, handles);

%% Parameters for analysis functions

%% Difussion original slope textbox
function slope_Callback(hObject, eventdata, handles)

%% Difussion original slope textbox properties
function slope_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Diffusion parameters
function difu_Callback(hObject, eventdata, handles)

%% Diffusion parameters
function difu_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Diffussion start time
function t1_Callback(hObject, eventdata, handles)

%% Diffussion start time properties
function t1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Diffussion time interval
function ti_Callback(hObject, eventdata, handles)

%% Diffussion time interval properties
function ti_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Diffusion endtime
function t2_Callback(hObject, eventdata, handles)

%% Diffussion endtime properties
function t2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Stack analysis peak search ratio
function Search_ratio_Callback(hObject, eventdata, handles)

%% Stack analysis peak search ratio properties
function Search_ratio_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Free-face slope parameter
function edit14_Callback(hObject, eventdata, handles)
handles.edit14=str2double(get(hObject,'String')) %returns contents of edit14 as a double
guidata(hObject, handles);

%% Free-face slope properties
function edit14_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Free-face slope range 
function edit15_Callback(hObject, eventdata, handles)
handles.edit15=str2double(get(hObject,'String')); %returns contents of edit15 as a double
guidata(hObject, handles);

%% Free-face slope range properties
function edit15_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Directories management buttons %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Read directories button
function pushbutton5_Callback(hObject, eventdata, handles)
handles = guidata(hObject);
%get handles
utmzone=get(handles.utmzonen,'String');
structdata = struct('terracem_dir', handles.maindir, 'stationsdir', handles.stationsdir,'zoneutm',utmzone); %define handles to save
addpath(handles.maindir) %set the main directory
save('terracem_path.mat', '-struct', 'structdata') %save path file as a structure

Commands=('Select the station to display the profiles');
T = evalc('Commands');
Tcell = regexp(T, '\n', 'split');
set(handles.text38,'String',Tcell);

%only put folder names new version
d=dir((handles.stationsdir)); %set the station directory
isub = [d(:).isdir]; % returns logical vector
nameFolds = {d(isub).name}'; %obtain fold names
nameFolds(ismember(nameFolds,{'.','..'})) = []; 
set(handles.listbox_1,'string', nameFolds) % define listbox_1 elements (Stations)

%create shorelines output table if not exist
cd(handles.stationsdir) %before with handles.maindir
%load output table
file=('shorelines.txt');

if exist(file,'file')==0 
        fid = fopen(file, 'a');     
        fprintf(fid,'Station Profile_number shoreline_number east north distance_along_swath shoreline_elevation error analysis_type time');

fclose(fid);
type 'shorelines.txt' 
end

guidata(hObject, handles);

%% UTM zone textbox 
function utmzonen_Callback(hObject, eventdata, handles)
Commands=('Set the UTM zone of the dataset');
T = evalc('Commands');
Tcell = regexp(T, '\n', 'split');
set(handles.text38,'String',Tcell);

guidata(hObject, handles);

%% UTM zone textbox properties
function utmzonen_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Set utm zone sub-GUI interface
function utmbutton_Callback(hObject, eventdata, handles)
handles = guidata(hObject);
%call UTM GUI
zone=utmzone;

if numel(zone)==0    
    disp('load aborted') 
    Commands=('UTM not defined');
    T = evalc('Commands');
    Tcell = regexp(T, '\n', 'split');
    set(handles.text38,'String',Tcell);
else
%separate strings
if numel(zone)==3
    A=zone(1);
    B=zone(2);
    C=zone(3);
    utmzt=horzcat([A,B,' ',C]);
else
    A=zone(1);
    B=zone(2);  
    utmzt=horzcat(['0',A,' ',B]);
end
    set(handles.utmzonen,'String',utmzt)
end

guidata(hObject, handles);

%% select profile interface
%function profnum_Callback(hObject, eventdata, handles)

%% select profile properties
%function profnum_CreateFcn(hObject, eventdata, handles)
%if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%    set(hObject,'BackgroundColor','white');
%end

%% Select station listbox
function listbox_1_Callback(hObject, eventdata, handles)
Commands=('Select a station to display it contents. If profles are not displayed use the "Swath Profile Extraction"');
T = evalc('Commands');
Tcell = regexp(T, '\n', 'split');
set(handles.text38,'String',Tcell);
handles = guidata(hObject);
handles.output = hObject;
%store the selected station name from listbox as station (input for
%terracem scripts)
contents = cellstr(get(hObject,'String')); %returns stationslist contents as cell array
station_selection=char(contents{get(hObject,'Value')}) %returns selected item from stationslist
handles.statsel=station_selection;
UTM_zone=get(handles.utmzonen,'String');
%put profile names of profiles in LISTBOX2%%%%%%%%%%
station=handles.statsel; %selected station
dirstations=handles.stationsdir; %stations directory 
try
cd(dirstations);
cd(station);
catch

Commands=('To run TerraceM set first the directories ');
T = evalc('Commands');
Tcell = regexp(T, '\n', 'split');
set(handles.text38,'String',Tcell);
    
end

if exist(sprintf('%s_swaths.mat',station))~=0
g=load(sprintf('%s_swaths.mat',station)); % load swaths.mat dataset
as=g.nim;  %number of profiles per station
handles.nim=as; %add number of profiles to handles
af=(1:as);
a=num2cell(af); %transform number to cell
set(handles.listbox_2,'String',a); %set the string inputs for listbox 2
%load all boxes in black colour
file=sprintf('%s_clip.shp',station); %search for profiles in shapefile format
S3 = shaperead(file,'UseGeoCoords', true);
nim=numel(S3);

for i=1:nim %polygon entry
lon11=S3(i,1).Lon; %utm
lon11=lon11';
lat11=S3(i,1).Lat; %utm
lat11=lat11';
bb=S3(i,1).BoundingBox;%in utm
xy=[lon11 lat11];
xy(6,:)=[]; %remove redundant values of shp

for k=1:5; 
    utmzo2(k,:)=(UTM_zone); %zone
end
zone=utmzo2;
%utm2deg function
[LatF LonF]=utm2deg(xy(:,1),xy(:,2),zone);
%save as an structure
SF(i,1).Lon=LonF';
SF(i,1).Lat=LatF';
handles.SF=SF;
end

%plot rectangular profiles in 2D latlong system
for tk=1:nim   
obj.h=plot(SF(tk,1).Lon,SF(tk,1).Lat,'-k'); 
hold on
end

set(gca,'xtick',[])
set(gca,'ytick',[])
refreshdata(obj.h) %refresh data if plot changes
set(gca,'Position', [4 8.7 111 44]);
axis equal
%avoid listbox_2 disapearing when changing station
selection = get(handles.listbox_2, 'Value')
    if as < selection
    set(handles.listbox_2, 'Value', 1);
    end    
end
guidata(hObject, handles);

%% Select profile Listbox
function listbox_2_Callback(hObject, eventdata, handles)
handles = guidata(hObject);
handles.output = hObject;
%store the selected profile name from listbox as station 
contents = cellstr(get(hObject,'String')); %returns profile contents as cell array
profile_selection=char(contents{get(hObject,'Value')}) %returns selected item from stationslist
handles.profsel=profile_selection;   
%plot selected profile in red colour
station=handles.statsel;
dirstations=handles.stationsdir;
%change to station directory

try
cd(dirstations);
cd(station);
catch
Commands=('In order to run TerraceM set first the directories');
T = evalc('Commands');
Tcell = regexp(T, '\n', 'split');
set(handles.text38,'String',Tcell);
end

file=sprintf('%s_clip.shp',station); 
try
S3 = shaperead(file,'UseGeoCoords', true);
nim=handles.nim;
UTM_zone=get(handles.utmzonen,'String');

for i=1:nim %polygon entry
lon11=S3(i,1).Lon; %
lon11=lon11';
lat11=S3(i,1).Lat; %
lat11=lat11';
bb=S3(i,1).BoundingBox;%in utm
xy=[lon11 lat11];
xy(6,:)=[]; %remove redundant values of shp
%set UTM zone matrix
for k=1:5; 
    utmzo2(k,:)=(UTM_zone); %zone
end

%Convert everything to LatLong for plotting
zone=utmzo2;
[LatF LonF]=utm2deg(xy(:,1),xy(:,2),zone);
SF(i,1).Lon=LonF';
SF(i,1).Lat=LatF';
end

%manage profiles
for i=1:nim
    [SF(i,1).Id]=i; %assign an id    
if str2double(profile_selection) == SF(i,1).Id;
 handles.a=([SF(i,1).Lon]'); %selected swath x
 handles.b=([SF(i,1).Lat]'); %selected swath y
       
end
end

% Plot DEM if dem plot is selected
ax1=axis;
if (handles.demsel)==1 %check if dem was activated from menu
%call to axis of main plot
axis (ax1)
% prompt information
Commands=('DEM display');
T = evalc('Commands');
Tcell = regexp(T, '\n', 'split');
set(handles.text38,'String',Tcell);
%DEM loading procedure
%load &/or combine structures
DEMS=(sprintf('DEM_%s.mat',station));
if exist((DEMS),'file')==2  
load(DEMS)
X=DEMS.X;
Y=DEMS.Y;
Z=DEMS.Z;
xa=X'; 
xa(:,2)=Y(1,1);
ka=numel(xa(:,1));
    for f=1:ka; 
    xut(f,:)=(UTM_zone); %generate  UTM zone matrix
    end
ya=Y;
ya(:,2)=X(1,1); %generate sinthetic utm coordinates for LatLon transsformation
ky=numel(ya(:,1)); 
    for k=1:ky; 
    yut(k,:)=(UTM_zone); %generate  UTM zone matrix
    end
% Tranform UTM to LatLon
[xLa xLo]=utm2deg(xa(:,1),xa(:,2),xut);
[yLa yLo]=utm2deg(ya(:,2),ya(:,1),yut);
%plot using Beaducel (2014) DEM shading routine
Zf=flipud(DEMS.Z);
hold on
%user defined inputs
N=handles.resoluti; %resampling 1/N
Az=handles.azimuti; %Azimut
[H I]=dem(xLo,yLa,Zf,'Azimuth',Az,'Contrast',0.5,'Lake',...
    'Decim',N,'Watermark',0,'LandColor',landcolor,'AxisEqual','off','NaNColor',[1 1 1],'ZCut',0);
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'Position', [4 8.7 111 44]);
refreshdata(H)
%prompt message
Commands=('To switch between different displays Clean the display interface (Diplays: DEM, satellite, and shaded )');
T = evalc('Commands');
Tcell = regexp(T, '\n', 'split');
set(handles.text38,'String',Tcell);
else    
%prompt message
Commands=('Error!! check the existence of XXX_DEM.mat in your stations folder path, otherwise run "Extract Swath profiles" again');
T = evalc('Commands');
Tcell = regexp(T, '\n', 'split');
set(handles.text38,'String',Tcell);
end
end

% if plot shaded selected
if (handles.shadsel)==1 %check if shaded googlemaps was activated from menu
%set googlemaps style
param='&style=feature:landscape.natural.terrain&style=element:geometry.fill|color:0xffffff|saturation:-100|lightness:50|visibility:off|gamma:0.14';
%call plot googlemaps function
h1=plot_google_map('MapType','terrain','style',param);
%axis(ax1)         
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'Position', [4 8.7 111 44]);    
%prompt message
Commands=('To switch between different displays "Clean" the display interface (Diplays: DEM, satellite, and shaded )');
T = evalc('Commands');
Tcell = regexp(T, '\n', 'split');
set(handles.text38,'String',Tcell);    
end

% if satellite image is selected
if (handles.imsel)==1 %check if shaded googlemaps was activated from menu
%call plot googlemaps function
plot_google_map('MapType','satellite','style',0);        
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'Position', [4 8.7 111 44]);
axis equal      
%prompt message
Commands=('To switch between different displays Clean the display interface (Diplays: DEM, Satellite, and Shaded )');
T = evalc('Commands');
Tcell = regexp(T, '\n', 'split');
set(handles.text38,'String',Tcell);        
end

% plot profiles in all cases
for tk=1:nim   
obj.h=plot(SF(tk,1).Lon,SF(tk,1).Lat,'-k','LineWidth',1); %plot shapefiles in 2D latlong system
hold on
end

set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'Position', [4 8.7 111 44]);
%plot selected profile in red
obj.i=plot(handles.a,handles.b,'-r','LineWidth',1); %'LineWidth',2
handles.objh=(obj.h);
handles.obji=(obj.i);
set(gca,'xtick',[])
set(gca,'ytick',[])
%prompt message
Commands=('To visualize a selected porfile use the "Profile plotter" button. The analysis functions allows analyzing the selected profile and estimate the shoreline angle');
T = evalc('Commands');
Tcell = regexp(T, '\n', 'split');
set(handles.text38,'String',Tcell);

% load shorelines table
if true
cd(dirstations);
file=('shorelines.txt');
fid=fopen(file,'r');
i=1;
while ~feof(fid)
A=textscan(fid,'%s %u %u %f %f %f %f %f %s %f/','delimiter',' ','HeaderLines',1);    
   for j=1:10 %number of columns (default)
       data(i,j) = A(1,j);
   end
i=i+1;
end
fclose(fid);
clear j i 
shore_pts(:,1)=cell2mat(data(:,4));
shore_pts(:,2)=cell2mat(data(:,5));
%convert shorelines to latlong
nix=numel(shore_pts(:,1));

if nix>0
for k=1:nix; 
    utmzo(k,:)=(UTM_zone); %zone
end

zone=utmzo;
[LatFp LonFp]=utm2deg(shore_pts(:,1),shore_pts(:,2),zone);
SP(:,1)=LonFp;
SP(:,2)=LatFp;
handles.shpoint=SP;
obj.p=plot(LonFp(:,1),LatFp(:,1),'o','MarkerEdgeColor','k','MarkerFaceColor','b');
axis (ax1)
hold off
set(gca,'xtick',[])
set(gca,'ytick',[])
end
end
%save data
guidata(hObject, handles);
catch
    Commands=('To run TerraceM set first the directories ');
T = evalc('Commands');
Tcell = regexp(T, '\n', 'split');
set(handles.text38,'String',Tcell);
end

%% Select profile listbox properties
function listbox_2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Select stations listbox properties
function listbox_1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% PROCESSING AND ANALYSIS FUNCTIONS %%

%% Extract swaths 
% Confirm swath extraction checkbox
function sureornot_Callback(hObject, eventdata, handles)
handles = guidata(hObject);
    Commands=('Prepare to extract swath profiles. This function also to update modifications of swath profile shapefiles');
    T = evalc('Commands');
    Tcell = regexp(T, '\n', 'split');
    set(handles.text38,'String',Tcell);

sure=get(hObject,'Value');
handles.sureornot=sure;
%save data to handles
guidata(hObject, handles);

%% Swath profile extraction button
function swatbuton_Callback(hObject, eventdata, handles)
handles = guidata(hObject);

if handles.sureornot==1  
%read inputs
path = 'terracem_path.mat';

if exist(path,'file')
     load(path);
else msgbox('Problem with path, enter again the input paths')
end

% select station 
station=handles.statsel
dirstations=handles.stationsdir;
Commands=('Please wait, processing data');
T = evalc('Commands');
Tcell = regexp(T, '\n', 'split');
set(handles.text38,'String',Tcell);
cd(dirstations);
cd(station);
% load asciigrid or geotiff
asc_file=sprintf('%s_dem.asc',station); 
tiff_file=sprintf('%s_dem.tiff',station);
cdx=pwd;

if exist(fullfile(cdx,asc_file),'file') == 2; %if ascii exist load it
    [Z,R] = arcgridread(asc_file);
    Z2=flipud(Z);
    dd=R(2,1); [ny nx]=size(Z);
    x1=R(3,1); x2=R(3,1)+dd*(nx-1);
    y2=R(3,2); y1=R(3,2)-dd*(ny-1);
    [X,Y]=meshgrid(x1:dd:x2,y1:dd:y2);
    X2=X;
    Y2=Y;
    X=X(1,:);
    Y=Y(:,1);
else %search for geotiff if ascii dont exist
[Z,R] = geotiffread(tiff_file);
Z = double(Z);
Z=flipud(Z);
Z(Z == -9999) = NaN; %outliers
Z(Z == -32767) = NaN;
dd=R.DeltaX; 
[ny nx]=size(Z);
x1=R.XLimWorld(1,1); 
x2=R.XLimWorld(1,1)+dd*(nx-1);
y2=R.YLimWorld(1,2); 
y1=R.YLimWorld(1,2)-dd*(ny-1);        
[X,Y]=meshgrid(x1:dd:x2,y1:dd:y2);  
X2=X;
Y2=Y;
X=X(1,:);
Y=Y(:,1);
end    

%save formated dem
DEMS=(sprintf('DEM_%s.mat',station));

%load table &/or combine structures
if exist(fullfile(DEMS,pwd),'file')==2    
load(DEMS)
DEMS.X=X;
DEMS.Y=Y;
DEMS.Z=Z;
save((sprintf('DEM_%s.mat',station)),'DEMS'); 
%in case dont exist, create file
else
DEMS.X=X;
DEMS.Y=Y;
DEMS.Z=Z;
save((sprintf('DEM_%s.mat',station)),'DEMS');    
end

file=sprintf('%s_clip.shp',station); %load the .shp rectangular profile areas
p=shaperead(file); 
nim=numel(p);
% Clip topography data with polygons
% Extract corners for each rectangular profile and clip with topography

Q = waitbar(0.1,'Please wait... procesing the swath data','windowstyle', 'modal');
frames = java.awt.Frame.getFrames();
frames(end).setAlwaysOnTop(1); 
for i=1:nim    
    data=[];
    waitbar(i / nim);
   
    data=p(i,1);     
    x=data.X;       %extract x values from shapefile
    x=(x');         %transpose elements
    x(end)=[];
    y=data.Y;       %extract y values from shapefile
    y=(y');         %transpose elements
    y(end)=[]; 
    xy=[x y];       %conform the resulting matrix
    xy(end,:)=[];   %delete the redundant last row and fix a 4 vertices polygon    
    taga=num2str(i); %assign names and vertice coordinates for each polygon
    v = genvarname(['Prof' taga]);
    eval([v '=xy;']);   
    in = inpolygon(X2,Y2,xy(:,1),xy(:,2)); %clip only X Y data of topography        
      
    %create cropped xyz for rotation input
    xyz(:,1)=X2(in);
    xyz(:,2)=Y2(in);
    xyz(:,3)=Z2(in);    
    %call function to extract swaths
    [swat,pto,thetaB,Dprof,ZI]=Terracem_extract_swaths(xy,xyz,station,i);    
    clear xyz xy in
    taga=num2str(i);
    v = genvarname(['swat' taga]); %Asign names for each cliped area
    eval([v '=swat;']); %Asign xyz values for each cliped area
    v = genvarname(['thetaB' taga]); %angle longest side
    eval([v '=thetaB;']); %angle of longest side
    v = genvarname(['pto' taga]); %point
    eval([v '=pto;']); %Asign minor xy point values for each cliped area
    v = genvarname(['Dprof' taga]); %point
    eval([v '=Dprof;']); %Asign minor xy point values for each cliped area    
    v = genvarname(['ZI' taga]); % mesh swath
    eval([v '=ZI;']); %  
    fprintf('...Extracting profile %u / %u...\n',i,nim)
end
close(Q);
clear DEMS frames data X2 Y2 Xp Xt Yp Yt h lati loni xy x y taga a out ex i swat Q R X Y Z Zc dd file nx ny v x1 x2 y1 y2 Commands T Tcell asc_file
clear cdx dirstations eventdata path p stationsdir terracem_dir tiff_file zoneutm pto thetaB Dprof ZI
%clearvars -except swat thetaB pto Dprof

% Save workspace for shoreline fit routine
save(sprintf('%s_swaths.mat',station),'-regexp','^(?!(hObject|handles|station|profsel)$).');

end

if handles.svpdf==1; %check if save swath to PDF is required
SWW=load(sprintf ('%s_swaths.mat',station)); 
%nim=numel(p)
%nim2=numel(SWW)
ps=3;% plots per page 
pg=ceil(nim/ps);%number of pages

%load swaths
for i=1:nim    
    SW=SWW.(sprintf('swat%u',i));    
    Xll=SW(:,1);%x
    Yll=SW(:,4);%max topography
    Yl2=SW(:,2);%minimum
    Yl3=SW(:,3);%mean    
    %loop subplot
    xt(i)=rem(i-1,ps)+1;    
    %generate new figure page after 3 plots  
    
    if xt(i)==1
    figure;
    end    
    ax1(i)=subplot(3,1,xt(i)); %loop subplot size 
    hold on
    plot(Xll,Yll,'-k');
    plot(Xll,Yl2,'-k');
    text(min(Xll)+200,max(Yll)-((max(Yll)-min(Yll))/4),num2str(i-1),'FontSize',18);
    plot(Xll,Yl3,'-b','linewidth',0.7); 
    xlim([min(Xll) max(Xll)]);   
    xlabel('Distance along profile (m)');
    ylabel('Elevation (m)');
    
    %set page dimensions
    if xt(i)==ps || i==nim
    rect=[1 0 19 29];
    set(gcf,'paperunits','centimeters');
    set(gcf,'papertype','A4');    
    set(gcf,'paperposition',rect);    
    fout = sprintf('Swaths_%s_page_%u',station,ceil(i/ps));
    saveas(gcf,fout,'pdf');
    end     
end     
end    
guidata(hObject, handles);

%% Save swaths as pdf checkbox
function svpdf_Callback(hObject, eventdata, handles)
handles = guidata(hObject);
%promt information
    Commands=('Swath profiles will be saved in .pdf format');
    T = evalc('Commands');
    Tcell = regexp(T, '\n', 'split');
    set(handles.text38,'String',Tcell);
vaono=get(hObject,'Value');
handles.svpdf=vaono;
%save data to handles
guidata(hObject, handles);


%% CLASSIC TERRACE ANALYSIS
function pushbutton8_Callback(hObject, eventdata, handles)
handles = guidata(hObject);
%promt information
    Commands=('The user should indicate the terrace level to map');
    T = evalc('Commands');
    Tcell = regexp(T, '\n', 'split');
    set(handles.text38,'String',T);
    
    
    
station=handles.statsel; %selected station
dirstations=handles.stationsdir; %path to stations
cd(dirstations);
cd(station);


g=load(sprintf('%s_swaths.mat',station)); %load swaths.mat database
count=0; %number of shorelines 
profnum=str2double(handles.profsel); %selected profile
nim=handles.nim; %nim is the number of profiles for each station



%while nim>0	   
    batchmode=handles.batchmode;
    
    if batchmode==1
    %Graphic interface for shoreline level in batch mode analysis
    shorenum = str2double(inputdlg('Enter shoreline level for batch mode mapping. [0=exit]'));   
    
    for profnum=1:nim        
    swat=g.(sprintf('swat%u',profnum));
    thetaB=g.(sprintf('thetaB%u',profnum)); 
    pto=g.(sprintf('pto%u',profnum));  
    Dprof=g.(sprintf('Dprof%u',profnum)); 
    %promt information
    Commands=('Staricase analysis in batch mode allows performing a continuous analysis of swath profiles');
    T = evalc('Commands');
    Tcell = regexp(T, '\n', 'split');
    set(handles.text38,'String',T);
    %change directories
    cd(dirstations);
    cd(station);
    
    
    
    %call to function Sfitprof  
    shoreline=Terracem_sfitprof(swat,station,profnum,shorenum,pto,thetaB,Dprof,batchmode); 
    handles.t=now;
    %change directory    
    cd(dirstations); 
    %Read data file 
    file=('shorelines.txt');
    %write header if file not exist
    fid = fopen(file, 'a');    
    fprintf(fid,'\n%s %u %u %6.0f %7.0f %5.2f %4.2f %4.2f staircase %8.5f',station,profnum,shoreline,now);
    fclose(fid);
    type 'shorelines.txt' 
    %promt information   
    Commands=('The shoreline angle has been saved');
    T = evalc('Commands');
    Tcell = regexp(T, '\n', 'split');
    set(handles.text38,'String',Tcell);     
    count=count+1;
    profnum=profnum+1;   
    end    
    else    
    %Graphic interface for shoreline level
    shorenum = str2double(inputdlg(sprintf('CONTINUE..? [0=exit]. Enter shoreline level (profile %u has %u shorelines) ',profnum,count)));   
    if shorenum== 0
        return          
    else    
    % select inputs
    swat=g.(sprintf('swat%u',profnum));
    thetaB=g.(sprintf('thetaB%u',profnum)); 
    pto=g.(sprintf('pto%u',profnum));  
    Dprof=g.(sprintf('Dprof%u',profnum));         
    cd(dirstations);
    cd(station);
    
  %promt message
    Commands=('Staricase analysis: 4 points are used to determine the shoreline angle from right-to-left. (1-2) Along the palecliff zone 2 points should be marked with red crosshair; and (3-4) along the paleo platform two points with blue crosshair.');
    T = evalc('Commands');
    Tcell = regexp(T, '\n', 'split');
    set(handles.text38,'String',T);
    
    %call to the function Sfitprof %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    shoreline=Terracem_sfitprof(swat,station,profnum,shorenum,pto,thetaB,Dprof,batchmode); 
    handles.t=now;
    pause(2)    
    cd(dirstations); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Read data file 
    file=('shorelines.txt');
        fid = fopen(file, 'a');    
        fprintf(fid,'\n%s %u %u %6.0f %7.0f %5.2f %4.2f %4.2f staircase %8.5f',station,profnum,shoreline,now);    
    fclose(fid);
    type 'shorelines.txt' 
    disp('Shoreline Saved')
    %promt information
    Commands=('The shoreline has been saved');
    T = evalc('Commands');
    Tcell = regexp(T, '\n', 'split');
    set(handles.text38,'String',Tcell);
    count=count+1;
    end
    end
%end

%% Batchmode checkbox for staircase analysis
function batchmode_Callback(hObject, eventdata, handles)
%promt information    
    Commands=('Staricase analysis in batch mode allows performing a continuous analysis of swath profiles');
    T = evalc('Commands');
    Tcell = regexp(T, '\n', 'split');
    set(handles.text38,'String',Tcell);

batchmode=get(hObject,'Value')
handles.batchmode=batchmode;
%save data to handles
guidata(hObject, handles);

%% Cliff Free-face analysis call button
function fxcliff_Callback(hObject, eventdata, handles)
handles = guidata(hObject);
station=handles.statsel; %selected station
dirstations=handles.stationsdir; %path to stations
cd(dirstations);
cd(station);
g=load(sprintf('%s_swaths.mat',station)); %load swaths.mat database
count=0; %number of shorelines 
profnum=str2double(handles.profsel); %selected profile
%promt information
    Commands=('The user should indicate the terrace level to map');
    T = evalc('Commands');
    Tcell = regexp(T, '\n', 'split');
    set(handles.text38,'String',Tcell);    
nim=handles.nim; %nim is the number of profiles for each station

%while nim>0	    
    shorenum = str2double(inputdlg(sprintf('[0=exit]. Enter shoreline number (profile %u has %u shorelines) ',profnum,count)));    
    if shorenum== 0
        return          
    else          
    swat=g.(sprintf('swat%u',profnum));
    thetaB=g.(sprintf('thetaB%u',profnum)); 
    pto=g.(sprintf('pto%u',profnum));  
    Dprof=g.(sprintf('Dprof%u',profnum));  
    slo=handles.edit15;
    ra=handles.edit14;   
    cd(dirstations);
    cd(station);
    %promt information
    Commands=('Cliff Free Face: Extrapolations of a custom cliff-slope and paleo-platform are intersected by defining 3 points form right-to-left: (1) The coluvium-bedrock contact at the cliff with red crosshair (2) Two points defining the paleo platform with blue corsshair');
    T = evalc('Commands');
    Tcell = regexp(T, '\n', 'split');
    set(handles.text38,'String',T);     
    %call to function Sfitprof %%  
    shoreline=Terracem_fxcliff(swat,station,profnum,shorenum,pto,thetaB,Dprof,slo,ra); 
    handles.t=now;
    cd(dirstations); 
    %Read data file 
    file=('shorelines.txt');
    %write header if file not exist
    fid = fopen(file, 'a');    
    fprintf(fid,'\n%s %u %u %6.0f %7.0f %5.2f %4.2f %4.2f fxcliff %8.5f',station,profnum,shoreline,now);
    fclose(fid);
    type 'shorelines.txt' 
    %disp('Shoreline Saved')
    Commands=('The shoreline have been saved, congratulations ^_^');
    T = evalc('Commands');
    Tcell = regexp(T, '\n', 'split');
    set(handles.text38,'String',Tcell);    
    count=count+1;       
    end
%end

%% Stack analysis call button
function pushbutton9_Callback(hObject, eventdata, handles)
handles = guidata(hObject);
station=handles.statsel; %selected station
dirstations=handles.stationsdir; %path to stations
cd(dirstations);
cd(station);
g=load(sprintf('%s_swaths.mat',station)); %load swaths.mat database
count=0; %number of shorelines 
profnum=str2double(handles.profsel); %selected profile
nim=handles.nim; %nim is the number of profiles for each station

Commands=('Detection of sea-stacks and extrapolation towards the cliff. (1) The stacks zone is defined by two points from right-to-left with blue crosshair. (2) Levels or fringes from bottom-to-top with the green crosshair. (3) Slope inflection at the cliff (black corsshair)');
    T = evalc('Commands');
    Tcell = regexp(T, '\n', 'split');
    set(handles.text38,'String',Tcell);  

while profnum>0   
while nim>0   
% select inputs
    swat=g.(sprintf('swat%u',profnum));
    thetaB=g.(sprintf('thetaB%u',profnum)); 
    pto=g.(sprintf('pto%u',profnum));  
    Dprof=g.(sprintf('Dprof%u',profnum));          
    cd(dirstations);
    cd(station);  
% call to the function Stacks Analysis
    peak=str2num(get(handles.Search_ratio,'String')); %peak search ratio input    
    [stacks]=Terracem_stacks(swat,station,profnum,pto,thetaB,Dprof,peak);
handles.t=now;      
cd(dirstations); 


%Read shorelines data file 
    file=('shorelines.txt'); 
    
    if exist (file,'file')==0
    sfprint('Error, the shoreline file has been deleted by error, please load the directories again')    
    else        
    numb=numel(stacks(:,1)); 
    
    for i=1:numb;    
    fid = fopen(file, 'a');    
    fprintf(fid,'\n%s %u %u %6.0f %7.0f %5.2f %4.2f %4.2f stack %8.5f',station,profnum,stacks(i,:),now);
    end   
    fclose(fid);
    end
    
    type 'shorelines.txt' 
    disp('Shoreline Saved')
    %promt information
    Commands=('The shoreline have been saved :p');   
    T = evalc('Commands');
    Tcell = regexp(T, '\n', 'split');
    set(handles.text38,'String',Tcell);
    count=count+1;  
%save secondary outputs table for each station (check vertcat)  
     cd(station)    
     file2=sprintf ('%s_stackpeaks.mat',station);
     return     
end
end            
guidata(hObject, handles);

%% Difussion analysis call button
function difusionbuton_Callback(hObject, eventdata, handles)
handles = guidata(hObject);
station=handles.statsel;
dirstations=handles.stationsdir;
cd(dirstations)
cd(station)
g=load(sprintf('%s_swaths.mat',station));
count=0; %number of shorelines 
profnum=str2double(handles.profsel);
nim=handles.nim;
Commands=('The user should indicate the terrace level to map');
    T = evalc('Commands');
    Tcell = regexp(T, '\n', 'split');
    set(handles.text38,'String',Tcell);  
  

%while profnum>0    
%while nim>0	  
    %Graphic interface
    shorenum = str2double(inputdlg(sprintf('Enter shoreline level (profile %u has %u shorelines)',profnum,count)));
    
    
    
    
    if shorenum==0
          return         
    else   
    % select inputs
    swat=g.(sprintf('swat%u',profnum));
    thetaB=g.(sprintf('thetaB%u',profnum)); 
    pto=g.(sprintf('pto%u',profnum));  
    Dprof=g.(sprintf('Dprof%u',profnum));      
    %diffusion parameters
    thetadegrees=str2num(get(handles.slope,'String'));
    k=str2num(get(handles.difu,'String'));
    t1=str2num(get(handles.t1,'String'));
    ti=str2num(get(handles.ti,'String'));
    t2=str2num(get(handles.t2,'String'));
    %change directory
    cd(dirstations);
    cd(station);
    %promt information
    Commands=('The analysis area should be set adjusting the zoom. (1) One point is defined at the top of paleo-cliff, (2) The paleo-platform is defined by two points from right-to left, (3) finally one point defines the center of profile, near the dashed line.');
    T = evalc('Commands');
    Tcell = regexp(T, '\n', 'split');
    set(handles.text38,'String',T);    
    %call fit function
    [shoreline Difusion]=Terracem_terracediff(swat,station,profnum,shorenum,pto,thetaB,Dprof,thetadegrees,k,t1,ti,t2);
    handles.t=now;    
    cd(dirstations);
    %Read outputs table file 
    file=('shorelines.txt');
    fid = fopen(file, 'a');    
    fprintf(fid,'\n%s %u %u %6.0f %7.0f %5.2f %4.2f %4.2f diffusion %8.5f',station,profnum,shoreline,now);    
    fclose(fid);
    type 'shorelines.txt' 
    disp('Shoreline Saved')
    %promt information
    Commands=('The shoreline have been saved');
    T = evalc('Commands');
    Tcell = regexp(T, '\n', 'split');
    set(handles.text38,'String',Tcell);   
    count=count+1;
    
    %save secondary outputs table   
    cd(station)    
    file2=sprintf ('%s_diffusion.mat',station);
   
    if exist(file2,'file')==0 
    DIFFUSION=Difusion;     
    save(file2,'DIFFUSION')
    else
    load(file2);        
    DIFFUSION=vertcat(DIFFUSION,Difusion); 
    save(file2,'DIFFUSION');
    end      
   end
%end
%end
guidata(hObject, handles);
    



%% Post-processing/export buttons 
%% table edition interface
function table_edition_Callback(hObject, eventdata, handles)

dirstations=handles.stationsdir; %path to stations
cd(dirstations);

Commands=('This table contains the shoreline angle outputs. The input file of the table -shorelines.txt- is stored in the stations folder');
T = evalc('Commands');
Tcell = regexp(T, '\n', 'split');
set(handles.text38,'String',Tcell);

handles = guidata(hObject);
%call table editing function
cnames = {'Station', 'Profile_number', 'shoreline_number', 'east', 'north', 'distance_along_swath', 'shoreline_elevation', 'error', 'analysis_type', 'time'};%column names
myTable('shorelines.txt',cnames)

%% open diffusion results table
function tablekt_Callback(hObject, eventdata, handles)
% hObject    handle to tablekt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dirstations=handles.stationsdir; %path to stations
cd(dirstations);

Commands=('This table contains the shoreline angles and diffusion analysis outputs. The input file for this table -diffusion.txt- is stored in the station folder');
T = evalc('Commands');
Tcell = regexp(T, '\n', 'split');
set(handles.text38,'String',Tcell);

station=handles.statsel;
cd(station)
try
load(sprintf('%s_diffusion.mat',station))

% set variables 
profnum=DIFFUSION(:,1);
shorenum=DIFFUSION(:,2);
x=DIFFUSION(:,3);
y=DIFFUSION(:,4);
z=DIFFUSION(:,5);
ze=DIFFUSION(:,6);
kt=DIFFUSION(:,7);
rms=DIFFUSION(:,8);
a=DIFFUSION(:,9);
b=DIFFUSION(:,10);

%Read data file  
if exist('difussion.txt','file') ==2
%delete file difussion.txt
delete('difussion.txt')
    
file1=('difussion.txt');   
for i=1:numel(DIFFUSION(:,1))
    
    titles={'station', 'profile', 'shoreline', 'X','Y','Elevation (Z)','Ze','Kt','RMS','a'};
    
        fid = fopen(file1, 'a');    
        fprintf(fid,'\n%s %u %u %6.0f %7.0f %4.2f %4.2f %4.2f %4.2f %4.2f,%4.2f',station,DIFFUSION(i,:)');
end    
    fclose(fid);
    type 'difussion.txt' 
    disp('Difussion data Saved')
    
    Commands=('Difussion data Saved');
    T = evalc('Commands');
    Tcell = regexp(T, '\n', 'split');
    set(handles.text38,'String',Tcell);
    %save table
    myTable('difussion.txt',titles)
else

file1=('difussion.txt');   
for i=1:numel(DIFFUSION(:,1))
        fid = fopen(file1, 'a');    
        fprintf(fid,'\n%s %u %u %6.0f %7.0f %4.2f %4.2f %4.2f %4.2f %4.2f,%4.2f',station,DIFFUSION(i,:)');
end    
    fclose(fid);
    type 'difussion.txt' 
    disp('Difussion data Saved')
    
    Commands=('Difussion data Saved');
    T = evalc('Commands');
    Tcell = regexp(T, '\n', 'split');
    set(handles.text38,'String',Tcell);
    
    
myTable('difussion.txt')
end
catch

msgbox('No scarp difussion have been calculated yet')
end    
%% Cretae data repository  in PDF
function pushbutton32_Callback(hObject, eventdata, handles)

%set directories
station=handles.statsel;
dirstations=handles.stationsdir;
cd(dirstations)
cd(station)

Commands=('The output plots of shoreline angle fits may be exported to .pdf files');
T = evalc('Commands');
Tcell = regexp(T, '\n', 'split');
set(handles.text38,'String',Tcell);

%load .mat container of figures
load(sprintf ('%s_fits.mat',station)); 
fields = (fieldnames(Fits)); %names of figures
nim=numel(fields); %numel figures in directory

ps=8;% plots per page
pg=ceil(nim/ps);%number of pages

%decompose names using textscan & detect profnum and shorenum values
for i=1:nim
    str=fields{i,1};
    all(i,:)=textscan(str,'%s','delimiter','_');
    profnum=str2double(all{i,1}{2,1}); %identifier of profile
    shorenum=str2double(all{i,1}{3,1}); %identifier of shoreline
  
    %define structures for each figure based on shorenum and profnum
    myFigstruct=Fits.(sprintf('%s_%u_%u_fit',station,profnum,shorenum));
try
    Xll(i,:)=myFigstruct.children(1,1).properties.XLim; %xlim from structure
    Yll(i,:)=myFigstruct.children(1,1).properties.YLim; %ylim from structure
catch
msgbox('mark the last unloaded fit of the DR again')    
end
    %loop subplot
    xt(i)=rem(i-1,ps)+1;
    
    %generate new figure page after 8 plots  
    if xt(i)==1
    figure;
    end
        
    ax1(i)=subplot(4,2,xt(i)); %loop subplot size   
    xlim(myFigstruct.children(1,1).properties.XLim)
    ylim(myFigstruct.children(1,1).properties.YLim)
    xlabel('Distance (m)'); ylabel('Elevation (m)')
    
    %copy fig handles to multiplot
    h1(i) = struct2handle(myFigstruct,0);       
    copyobj(allchild(get(h1(i),'CurrentAxes')),ax1(i)); 
    
    xlabel('Distance along profile (m)')
    ylabel('Elevation (m)')
    close(h1(i))
    
    %set page dimensions
    if xt(i)==ps || i==nim
    rect=[1 0 19 29];
    set(gcf,'paperunits','centimeters');
    set(gcf,'papertype','A4');    
    set(gcf,'paperposition',rect);    
    fout = sprintf('Data-Repository_%s_%u',station,ceil(i/ps));
    saveas(gcf,fout,'pdf');
    end    
end
clear xt it ps np pg fout 

%% generate histograms (shoreline statistics)
function shore_statistic_Callback(hObject, eventdata, handles)
guidata(hObject, handles);

Commands=('The Statistics routine create histograms of shoreline-angle elevations for each station of the project, the figure is stored as a .pdf file in the station folder');
T = evalc('Commands');
Tcell = regexp(T, '\n', 'split');
set(handles.text38,'String',Tcell);

%path to stations
dirstations=handles.stationsdir; 
cd(dirstations);

%load table
file=('shorelines.txt');
if exist(file,'file')==1 || exist(file,'file')==2
    
%load shorelines table
fid=fopen(file,'r');
i=1;
while ~feof(fid)
    A=textscan(fid,'%s %u %u %f %f %f %f %f %s %f/','delimiter',' ','HeaderLines',1);
    
    for j=1:10 %number of columns (default of shorelines.txt)
        data(i,j) = A(1,j);
    end
    i=i+1;
end
fclose(fid);
%extract data from table
nim=numel(data(:,10)); u=1;

for k=1:nim
md1(u,1)=data{k,4}; %x
md1(u,2)=data{k,5}; %Y
md1(u,3)=data{k,7}; %z
md1(u,4)=data{k,8}; %ze
md1(u,5)=data{k,3}; %shore
%md1(u,6)=data{k,2}; %profnum
vars{u,1}=data{k,1}{1,1}; %station names

u=1+u;
end
end

name=(vars);
%shorelines clasification of names
name2=cellfun(@double,name,'uniformoutput',0);
code_shore=cell2mat(name2);
tab(:,1)=code_shore(:,1); %code_name  station
tab(:,2)=code_shore(:,2); %code_name  station
tab(:,3)=code_shore(:,3); %code_name  station
tab(:,4)=code_shore(:,4); %code_name  station
tab(:,5)=md1(:,3); %shoreline elevation
tab(:,6)=md1(:,4); %shoreline error
tab(:,7)=md1(:,5); %terrace shore
% loop to obtain all codes of avaliable stations
C = unique(tab(:,1:4),'rows');
numc=numel(C(:,1)); %number of name inputs & plots
for w=1:numc
%figure(w)  
r=1;
    for Q=1:nim %total elements in table   
    
    if tab(Q,(1:4)) == C(w,(1:4))
       vars2=(vars(Q,:));

       h_sh(r,:)=(tab(Q,5));       
       nx=numel(h_sh);      
       ppp(w).sh(r,:)=(tab(Q,5));
       ppp(w).st(r,:)=(vars2{1,1});
       r=r+1;    
    end
    end
end
cd(handles.maindir);
ps=6;% plots per page (par number)
np=numc; %total number of plots

for it=1:np
    xt=rem(it-1,ps)+1;
    if xt==1
    figure;
    end
    %subplots
    ax1=subplot(ps/2,2,xt);
    hold on
    hist(ppp(it).sh,sqrt(numel(ppp(it).sh)))
    hold on
    na=(ppp(it).st(1,:));
    %histogram statistics
    vline(mean(ppp(it).sh),'r')
    vline(mean(ppp(it).sh)-(std(ppp(it).sh)/2),'r--')
    vline(mean(ppp(it).sh)+(std(ppp(it).sh)/2),'r--')    
    texto1=num2str(mean(ppp(it).sh));
    text(mean(ppp(it).sh+0.5),30,texto1)
    xlabel('Elevation (m)')
    ylabel('Number of shorelines')
    title(sprintf('Shoreline statistics of station %s',na));  
    cd(dirstations);    
    %save to pdf routine
    if xt(it)==ps || it==np
    rect=[1 0 19 29];
    set(gcf,'paperunits','centimeters');
    set(gcf,'papertype','A4');    
    set(gcf,'paperposition',rect);     
    fout = sprintf('Station-statistics_%u',ceil(it/ps));
    saveas(gcf,fout,'pdf');
    end 
end

%% export in shapefile format
function shapexport_Callback(hObject, eventdata, handles)
%set directories
dirstations=handles.stationsdir; %path to stations
cd(dirstations);
Commands=('Shoreline angle outputs may be exported to ESRI shapefile format.');
T = evalc('Commands');
Tcell = regexp(T, '\n', 'split');
set(handles.text38,'String',Tcell);
%read table
file=('shorelines.txt');
fid=fopen(file,'r');
data=textscan(fid,'%4s %f %f %f %f %f %f %f %s %f','HeaderLines',1);
fclose(fid);
% Cell2Mat
id(:,1)=data{1,1}(:,1);
pnum(:,1)=data{1,2}(:,1);
snum(:,1)=data{1,3}(:,1);
j=1;

for i=4:8
    xydze(:,j)=data{1,i}(:,1); 
    j=j+1;
end

ana(:,1)=data{1,9}(:,1); 
%Export shape
[Tracks(1:length(xydze(:,1))).Geometry] = deal('Point');
%utmzonen='18 H';
utmzone=get(handles.utmzonen,'String');

for i = 1:length(xydze(:,1))
    [Tracks(i).station]= id{i,1};
    [Tracks(i).profile]= pnum(i,1);
    [Tracks(i).shoreline]= snum(i,1);
    [Tracks(i).X]= xydze(i,1);
    [Tracks(i).Y]= xydze(i,2);
    [Tracks(i).distance]= xydze(i,3);
    [Tracks(i).Z]= xydze(i,4);
    [Tracks(i).Ze]= xydze(i,5);        
    [Tracks(i).type] = ana{i,1};  
    [Lat,Lon] = utm2deg(xydze(i,1),xydze(i,2),utmzone);   
end

%outputs
fout=('shorelines_points.shp');
shapewrite(Tracks, fout);

%% export in google earth format
function exportkmz_Callback(hObject, eventdata, handles)
%set directories
dirstations=handles.stationsdir; %path to stations
cd(dirstations);
%promt information
Commands=('Shoreline angle outputs may be exported to GooGle Earth (.kml) format.');
T = evalc('Commands');
Tcell = regexp(T, '\n', 'split');
set(handles.text38,'String',Tcell);
%load table
file=('shorelines.txt');
fid=fopen(file,'r');
data=textscan(fid,'%4s %f %f %f %f %f %f %f %s %f','HeaderLines',1);
fclose(fid);
% Cell2Mat
id(:,1)=data{1,1}(:,1);
pnum(:,1)=data{1,2}(:,1);
snum(:,1)=data{1,3}(:,1);
j=1;

for i=4:8
    xydze(:,j)=data{1,i}(:,1); 
    j=j+1;
end

ana(:,1)=data{1,9}(:,1); %analysis method
time(:,1)=data{1,10}(:,1); %time of execution
%Export kml
[Tracks(1:length(xydze(:,1))).Geometry] = deal('Point');
utmzone=get(handles.utmzonen,'String');

for i = 1:length(xydze(:,1))
    [Lat,Lon] = utm2deg(xydze(i,1),xydze(i,2),utmzone);
    [Tracks(i).Lat]= Lat;
    [Tracks(i).Lon]= Lon;
    [Tracks(i).station]= id{i,1};
    [Tracks(i).profile]= pnum(i,1);
    [Tracks(i).shoreline]= snum(i,1);
    [Tracks(i).X]= xydze(i,1);
    [Tracks(i).Y]= xydze(i,2);
    [Tracks(i).distance]= xydze(i,3);
    [Tracks(i).Z]= xydze(i,4);
    [Tracks(i).Ze]= xydze(i,5);        
    [Tracks(i).type] = ana{i,1};  
end

fout=(sprintf('shorelines.kml'));
kmlwrite(fout, Tracks, 'Name', id);

%% Sub-GUI button calls %%%%%%%%%%%%%%

%% Swath profile display sub-GUI properties
function profile_Callback(hObject, eventdata, handles)

%% Swath profile display sub-GUI call button
function profile_button_Callback(hObject, eventdata, handles)
maindir=handles.maindir;
cd(maindir)
dirstations=handles.stationsdir;
station=handles.statsel;
profsel=handles.profsel; 
save('terracem_path.mat','station','profsel','-append')
Commands=('Swath previewer display swath profiles and shoreline angle fits. Swath profiles can be displayed as lineas, grey patch, or color-coded by stacked probability distribution functions. Shorelines-angles, if already calculated, can be displayed');
T = evalc('Commands');
Tcell = regexp(T, '\n', 'split');
set(handles.text38,'String',Tcell);
%call swath plotter sub-GUI
swath_plotter

%% SH filter sub-GUI call button
function pushbutton42_Callback(hObject, eventdata, handles)
    %promt information
    Commands=('This routine allows improving interpolations between shoreline angles by removing outliers. Shoreline angles should be projected before proceed. Interpolations and outliers can be saved.');
    T = evalc('Commands');
    Tcell = regexp(T, '\n', 'split');
    set(handles.text38,'String',Tcell);
cd(handles.stationsdir) 
file=('sh_projected.txt');

if exist(file,'file')==1 || exist(file,'file')==2;
fid=fopen(file,'r');
i=1;
while ~feof(fid)
    A=textscan(fid,'%u %u %u %u/','delimiter',' ','HeaderLines',1);
    
    for j=1:4 %number of columns (default from shorelines.txt)
        data(i,j) = A(1,j);
    end
    i=i+1;
end
%This function can run only if the number of extracted shoreline-angles is
%higher than 10
 if numel(data(:,1))>=10;
  SH_filter;  
 else
 %promt information
    Commands=('Insuficient shorelines have been mapped for the analysis, minimum 10 shoreline angles');
    T = evalc('Commands');
    Tcell = regexp(T, '\n', 'split');
    set(handles.text38,'String',Tcell);      
 end
end
 
%% project difussion values routine
function difuse_analysis_Callback(hObject, eventdata, handles)
cd(handles.stationsdir)
station=handles.statsel;
cd(station)

if exist(sprintf('%s_diffusion.mat',station))==0
disp('No diffusion calculated yet')   
%promt information
    Commands=('No difussion values have been analyzed yet');
    T = evalc('Commands');
    Tcell = regexp(T, '\n', 'split');
    set(handles.text38,'String',Tcell);
else    
%load difussion secondary table
load(sprintf('%s_diffusion.mat',station))
figure
ax1=subplot(2,1,1:2);
% set variables 
x=DIFFUSION(:,2);
y=DIFFUSION(:,3);
z=DIFFUSION(:,4);
ze=DIFFUSION(:,5);
dif=DIFFUSION(:,6);
rms=DIFFUSION(:,7);
a=DIFFUSION(:,8);
%plot values for projection
hold on
scatter(x,y,100,dif,'filled');
col=colorbar;
tite = get(col,'title');
set(tite,'string','KT (m^2)');
axis equal
box on
ylabel('UTM   N')
xlabel('UTM   E')
%graphic input profile
[xl,yl]=ginput(2);
l(:,1)=xl;
l(:,2)=yl;
plot(xl,yl,'-k') %plot of profile
plot(l(1,1),l(1,2),'or'); %point1
plot(l(2,1),l(2,2),'ob'); %point2
%angle of profile
dx=l(1,1)-l(2,1);
dy=l(1,2)-l(2,2);
m=(dy/dx);
m2=-1/m;
an=-atan(m)*180/pi;
anr=an*pi/180;
ddg=sqrt((l(1,1)-l(2,1))^2+(l(1,2)-l(2,2))^2);
dkm=(ddg);
b=l(1,2)-m*l(1,1);
%text in plot
str1=num2str(dkm/1000); %m to km
text(min(xl)+(dx/2),min(yl)+(dy/2),sprintf('%s km',str1),'FontSize',12)
%rotate and reproject data

for i=1:length(x(:,1))
    b2=y(i)-m2*x(i);
    Pnx=(b2-b)/(m-m2);
    Pny=m*Pnx+b;
    Pn=[Pnx Pny];
    dp(i,1)=Pn(1,1); %xp
    dp(i,2)=Pn(1,2); %yp
    dp(i,3)=((l(1,2)-dp(i,2))/sin(anr)); % distance in km
    dp(i,4)=z(i);
    dp(i,5)=ze(i);
    dp(i,6)=dif(i);
    dp(i,7)=rms(i);
    dp(i,8)=a(i);
end
%plot projected Kt values
figure(2)
clf
subplot(2,1,1)
hold on
box on
scatter(dp(:,3),dp(:,6),60,dif,'filled');
plot(dp(:,3),dp(:,6),'o')
%cartesian plane conditions for profiles
if dx>0 && dy>0   
xlim([-dkm/1000 0])
a=1; %a is a parameter for data export
elseif dx<0 && dy<0 
xlim([0 dkm/1000])    
a=3;
elseif dx>0 && dy<0  
xlim([-dkm/1000 0]) 
a=4;
else
a=2;
xlim([0 dkm/1000]) 
end
ylabel('Elevation (m)')
xlabel('Distance along profile (km)')
end

%% Point projection sub-GUI call button
function projection_button_Callback(hObject, eventdata, handles)
%promt information
    Commands=('Shoreline angles can be projected along custom lines or profiles. User defines the line by two points, or using a ESRI(R) shapefile containing the profile.');
    T = evalc('Commands');
    Tcell = regexp(T, '\n', 'split');
    set(handles.text38,'String',Tcell);
%set directories
cd(handles.stationsdir)
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
if numel(data(:,1))>=10
SH_projection
else 
    %promt information
    Commands=('Insuficient shorelines');
    T = evalc('Commands');
    Tcell = regexp(T, '\n', 'split');
    set(handles.text38,'String',Tcell); 
    
end
end






%% TERRACEM ANALYSIS AND PROCESSING FUNCTIONS %%%%%%%%%

%% Extract swaths function
function [swat,pto,thetaB,Dprof,ZI]=Terracem_extract_swaths(p,a,station,profnum)
%=======================================================================
% extract swath profiles from rectangular boxed drwan in arcgis
% Designed by by Melnick and Jara-Mu?oz
% inputs: cropped xyz topography and polygon
% UTM projection only
% Extracts multiple swath profiles from a shapefile with several polygons
% inputs
% p: Vectors with the corners of clip poligons
% a: cropped xyz file
% station: site name
% profnum: profile number
%========================================================================

format long;
X=p(:,1); Y=p(:,2);
%Define coordinates for each vertice
x1=p(1,1); x2=p(2,1); x3=p(3,1); x4=p(4,1);
y1=p(1,2); y2=p(2,2); y3=p(3,2); y4=p(4,2);
%minimum and maximum for the entire dataset
minx=min(X); miny=min(Y);
maxx=max(X); maxy=max(Y);
nim=numel(p);
% LENGHT
% Alternative 1, calculate all diagonals, the middle value should be the lenght
A(1,1)=sqrt((x2-x1)^2+(y2-y1)^2);
A(2,1)=sqrt((x3-x2)^2+(y3-y2)^2);
A(3,1)=sqrt((x3-x4)^2+(y3-y4)^2);
A(4,1)=sqrt((x4-x2)^2+(y4-y2)^2);
A(5,1)=sqrt((x3-x1)^2+(y3-y1)^2);
A(6,1)=sqrt((x4-x1)^2+(y4-y1)^2);
% respective coordinate of each measured diagonal is added to A, (coordinate dragging)
A(1,2)=x1'; A(1,3)=y1'; A(1,4)=x2'; A(1,5)=y2';
A(2,2)=x2'; A(2,3)=y2'; A(2,4)=x3'; A(2,5)=y3';
A(3,2)=x3'; A(3,3)=y3'; A(3,4)=x4'; A(3,5)=y4';
A(4,2)=x2'; A(4,3)=y2'; A(4,4)=x4'; A(4,5)=y4';
A(5,2)=x1'; A(5,3)=y1'; A(5,4)=x3'; A(5,5)=y3';
A(6,2)=x1'; A(6,3)=y1'; A(6,4)=x4'; A(6,5)=y4';
% preparing data to obtain rectangle diagonals and angles
A(:,1)=round(A(:,1)); %eliminate decimal places in the first column, this cause problems with near 180? theta value
length=A(:,1);
X1=A(:,2); Y1=A(:,3);
X2=A(:,4); Y2=A(:,5);
H=[length X1 Y1 X2 Y2]; %defines the matrix H which will be used only for length calculation
L=H;     
% Calculation of the lengh of the longest side
% select the min side length
minV=min(H(:,1)); % find min value in H
H(H==minV)=0; % set all min values to zero
maxV=max(H(:,1)); % find max value in H
H(H==maxV)=0; % set all max values to zero
ind=find(H(:,1)~=0, 1, 'first'); %finds the address of the first non-zero value
ddg=H(ind); %lenght of the rectangle "Requested answer"
Dprof=ddg; %for use in sfitprof
% Find coordinate pairs of the shortest and nearest to the origin
G=sortrows(L,1);
minimum(1,:)=G(1,:);
minimum(2,:)=G(2,:);
%extract vertices coordinates of the shortest side nearest to the origin
x1cord=minimum(1,2);%pair x1,y1
y1cord=minimum(1,3);%pair x1,y1
x2cord=minimum(1,4);%pair x2,y2
y2cord=minimum(1,5);%pair x2,y2
x3cord=minimum(2,2);%pair x3,y3
y3cord=minimum(2,3);%pair x3,y3
x4cord=minimum(2,4);%pair x4,y4
y4cord=minimum(2,5);%pair x4,y4
%distance
De1=(x1cord^2)+(y1cord^2);
De2=(x2cord^2)+(y2cord^2);
De3=(x3cord^2)+(y3cord^2);
De4=(x4cord^2)+(y4cord^2);
%detect shorter side
minimum(1,6)=De1;
minimum(1,7)=De2;
minimum(2,6)=De3;
minimum(2,7)=De4;
%second euclidian distance
minimum(1,8)=((minimum(1,6)^2)+(minimum(1,7)^2));
minimum(2,8)=((minimum(2,6)^2)+(minimum(2,7)^2));
% definition of the lower coordinates for the shortest side
Q=sortrows(minimum,8);
x1c=Q(1,2);
y1c=Q(1,3);
x2c=Q(1,4);
y2c=Q(1,5);
%determining x1 and x2, x1 must be the nearest point to the origin
QF(1,1)=x1c;
QF(1,2)=y1c;
QF(2,1)=x2c;
QF(2,2)=y2c;
%determination of which of the points of the shorter side is near to the
%origin
DF=sortrows(QF,2); % modified, now use Y as discriminant
x1cord=DF(1,1);
y1cord=DF(1,2);
x2cord=DF(2,1);
y2cord=DF(2,2);
% Rotation angle theta and the middle point of the shortest side 
if x1cord-x2cord ~= 0; %any case & 90?
    
    if ((y2cord>y1cord) && (x2cord>x1cord)); 
    g1=(atand((y2cord-y1cord)/(x2cord-x1cord))); %calculate angle of the shortest side positive slope
    theta=-g1+180; %requested angle theta for swath
    %point coordinate
    d1=(sqrt((y2cord-y1cord)^2+(x2cord-x1cord)^2))/2; %1/2 hypotenusa     
    xpto=x1cord+((cosd(g1))*d1); %check bug in g1
    ypto=y1cord+((sind(g1))*d1);%check bug in g1
    thetaB=(90+g1);
    pto=[xpto, ypto]; %requested angle (point location)

    else 
    g2=(atand((y1cord-y2cord)/(x2cord-x1cord))); %calculate angle of the shortest side with negative slope
    theta=g2; %theta requested answer for matrix rotation
    d2=(sqrt((y2cord-y1cord)^2+(x1cord-x2cord)^2))/2; %1/2 hypotenusa    
    xpto=x1cord+((cosd(g2+180))*d2); 
    ypto=y1cord-((sind(g2+180))*d2); 
    thetaB=(90-g2);    
    pto=[xpto, ypto]; %requested answer (point location)
        
    end    
end

if x1cord-x2cord == 0; %horizontal rectangle (180?)

    g1=(atand((y2cord-y1cord)/(x2cord-x1cord))); %calculate angle of the shortest side positive slope
    theta=-g1+180; %requested angle theta for swath
     %point coordinate
     d1=(sqrt((y2cord-y1cord)^2+(x2cord-x1cord)^2))/2; %1/2 hypotenusa     
     xpto=x1cord+((cosd(g1))*d1); 
     ypto=y1cord+((sind(g1))*d1);
     thetaB=(90+g1);
     pto=[xpto, ypto]; %requested angle (point location)          
end 
% rotation matrix
for i=1:(numel(a(:,1))-1); 
    aa = [a(i,1) ; a(i,2)]; %"aa" will be equal to i rows in the column 1, and the same for the column 2
    b = [cosd(theta)  -sind(theta) ; sind(theta)  cosd(theta)] * aa; %define b including the line angle
    ar=b';
    ro(i,1) = ar(1,1); 
    ro(i,2) = ar(1,2); 
    ro(i,3) = a(i,3); %keeps Z in the colum 3 intact
end
clear i

if exist('ro','var')==1;
%gridding rotated xyz to matrix
ddx=abs(ro(1,1)-ro(2,1)); %define ddx as the absolute value of the diference of the proyected x
ddy=abs(ro(1,2)-ro(2,2)); %define ddy as the absolute value of the diference of the proyected y
xlim_min=min(ro(:,1)); %define the boundaries of proyected x coordinate
xlim_max=max(ro(:,1));
ylim_min=min(ro(:,2)); %define the boundaries of proyected Y coordinate
ylim_max=max(ro(:,2));
x=(xlim_min:ddx:xlim_max); %define x and y vectors
y=(ylim_min:ddy:ylim_max);
[XI,YI]=meshgrid(x,y); %mesh with x and y
ZI=griddata(ro(:,1),ro(:,2),ro(:,3),XI,YI,'nearest'); %grid rotated data with nearest neighbors 
% Swath profile extraction
ma=max(ZI'); %maximum raster value
me=mean(ZI'); %mean raster value
mi=min(ZI'); %minimum raster value
% use projected latitude as x
[ny nx]=size(ZI);
% use km as x
xx2=0:ddg/(ny-1):ddg;%X for normal swath

% DEFINE VARIABLE swat 
%and check errors during run
try
swat(:,1)=xx2(1,:)';
swat(:,2)=mi(1,:)';
swat(:,3)=me';
swat(:,4)=ma';
swat(:,5)=abs(atand((gradient(ma))));%,mean(diff(xx2')))))); %slope
catch
    msgbox('Some profiles are outside the topography extent or are not rectangular')
end

h2=figure;
% plot swaths
subplot(3,1,1); 
hold on
U7=area(swat(:,1),swat(:,4));
U6=area(swat(:,1),swat(:,2));
set(U7,'FaceColor',[0.6 0.6 0.6],'edgeColor',[0.6 0.6 0.6]);
set(U6,'FaceColor',[1 1 1],'edgeColor',[0.6 0.6 0.6]);
set(gca,'Layer','top');
plot(swat(:,1),swat(:,3),'-k','linewidth',1.4);
axis([0 ddg min(swat(:,2)) max(swat(:,4))]);
xlabel('Distance along swath [m]');
ylabel('Elevation [m]');
tit=sprintf('Station %s, Profile %u',station,profnum);
title(tit);
grid on, box on
xlim([0 ddg]);

subplot(3,1,2); 
plot(swat(:,1),swat(:,5),'.r'); 
grid on, box on
xlim([0 ddg]);
ylabel('Slope (deg)'); xlabel('Distance along swath (m)');

subplot(3,1,3)
pcolor(YI,XI,ZI);
colormap (demcmap(ZI));
colorbar('location','east');
shading flat;

close(h2)
else   
msgbox('The topography not match with the extent of swath profiles, check input topography (..._dem.asc) or profiles (..._clip.shp), if the problem persist try deleting the file ...._swaths.mat from the station directory')    
return
end
    
%% Terracem_Free-face function
function shoreline=Terracem_fxcliff(swat,station,profnum,shorenum,pto,thetaB,Dprof,slo,ra) 
%==================================================================
% D. Melnick 2012, modified by Jara-Munoz, 2015
%
% Intersect cliff's free-face and interpolated abrasion platform to 
% find shoreline angle.
% 
% v. 1.8
%==================================================================
% Set free-face slope and error
fface(1,1)=ra; %Plio-Pleist. Konya limestone (70 5) (slo ra)
fface(1,2)=slo;
format long
%load swath profile values
x=swat(:,1); %distance along profile
mi=swat(:,2); %min values%
me=swat(:,3); %mean values
ma=swat(:,4); %max values
slp=swat(:,5); %slope values

% Check if flipped, high part must be at right
if ma(1)>ma(end)
   type_profile=2; %identifies if profile is flipped
    
   x=max(x)-x+1;   
else   type_profile=1; %identifies if profile is in correct position
 
    x=x;
end

%plot of topographic distributions from swaths
f1=figure;
hold on
ax1=subplot(3,1,1:2);
hold on
plot(x,ma,'k')
plot(x,me,'color',[0.3 0.3 0.3])
plot(x,mi,'color',[0.6 0.6 0.6])
hleg1 = legend('Max topography','Mean topography','Min topography','Location','SouthEast');
xlim([min(x) max(x)]);
grid on 
box on
yl1=ylabel('Elevation (m)'); 
xl1=xlabel('Distance along profile (m)');
tit=sprintf('Station %s, Profile %u, Shoreline %u',station,profnum,shorenum);
title(tit);
%plot of slope of maximum topography
ax2=subplot (3,1,3);
hold on
plot(x,slp,'.k') 
plot(x,moving(slp,50,'mean'),'r-');
grid on
box on
leyendita = legend('slope values','Location','NorthEast');
xlim([min(x) max(x)]); 
ylabel('Slope (deg)'); xlabel('Distance along profile (m)'); 
linkaxes([ax1 ax2], 'x'); 
% UNICONTROL RADIO BUTTON for select maximum or minimum topo for analysis
% Create the button group.
ht = uibuttongroup('visible','on','Position',[0 0 1 0.1],'Units','Normalized'); %radio buttons group
% Create two radio buttons inside the button group.
u1 = uicontrol('Style','radiobutton','String','Maximum swath topo','pos',[300 5 130 30],'parent',ht,'HandleVisibility','on');
u2 = uicontrol('Style','radiobutton','String','Minimum swath topo','pos',[100 5 130 30],'parent',ht,'HandleVisibility','on');
disp('Zoom in, select area, and press [Enter]')
%pause
%set(ht, 'Visible','off') %hide buttons when drawing paleocliff and paleo platform
u1p=get(u1,'Value'); %detect max topo radiobutton handles
u2p=get(u2,'Value'); %%detect min topo radiobutton handles

%check which button was pushed and associate with type of topography used
if u1p==1
answer=1; %max topography
else
answer=0; %minumum topography   
end

if exist(sprintf ('%s_Clicks.mat',station),'file')==2
load(sprintf ('%s_Clicks.mat',station));  %assumes that Clicks empty matrix was already saved  
else
Clicks=[]; %if clicks matrix not exist creates an empty  matrix for the next condition
end
%check structure existence
doesVarExist = true;
try
    Clicks.(sprintf('prof_%u',profnum)).(sprintf('shore_%u',shorenum)).fxcliff
catch
    doesVarExist = false;
end

if doesVarExist == 1
prompt1 = {'would you like to overwrite the data? Yes(0) / No (1)'};
dlg_title1 = 'rewrite or re-read the data?';
num_lines1 = 1;
def1 = {'1'}; %default values
answer1 = str2double(inputdlg(prompt1,dlg_title1,num_lines1,def1));

if answer1(1,1)==0; 
        existence=0; %parameter for data existence situation, existence=0 no previous data recognized.
        disp('enter outcrop position of paleocliff free face');
%input free-face
        [face(1,1),face(1,2)]=ginputc(1,'Color','r');       
        plot(face(1,1),face(1,2),'sk','markerfacecolor','k');
        disp('enter top and bottom of paleo abrasion platform');
        [plat(1,1),plat(1,2)]=ginputc(1,'Color','b');       
        [plat(2,1),plat(2,2)]=ginputc(1,'Color','b');       
else
        existence=1;  %parameter for data existence situation, existence=1 previous data recognized.  
        face=Clicks.(sprintf('prof_%u',profnum)).(sprintf('shore_%u',shorenum)).fxcliff.face;
        plat=Clicks.(sprintf('prof_%u',profnum)).(sprintf('shore_%u',shorenum)).fxcliff.plat;
end
else
% input clif
existence=0; 
disp('enter outcrop position of paleocliff free face');
%input free-face
[face(1,1),face(1,2)]=ginputc(1,'Color','r');
disp('enter top and bottom of paleo abrasion platform');
[plat(1,1),plat(1,2)]=ginputc(1,'Color','b'); 
[plat(2,1),plat(2,2)]=ginputc(1,'Color','b');
end

if answer == 1
da(:,1)=x; 
da(:,2)=ma; 
md=swat; 
ny=length(md(:,1)); 

else
da(:,1)=x; 
da(:,2)=mi;
da(isnan(da)) = 0; %convert nan values to zero
md=swat; 
ny=length(md(:,1));     
end

%correct inverted profile
if type_profile==2
da=flipud(da);  
end

% select points for fitting cliff
daclif=da(da(:,1)>=face(1,1),:);
% select points for fitting platform
daplat1=da(da(:,1)<=plat(1,1),:);
daplat=daplat1(daplat1(:,1)>=plat(2,1),:);


%linear regressions
% free-face and cliff
xc=(daclif(:,1));
yc=daclif(1,2)-((daclif(:,1)-daclif(1,1))./tand(90+fface(1,1)));
yc1=daclif(1,2)-((daclif(:,1)-daclif(1,1))./tand(90+fface(1,1)+fface(1,2)));
yc2=daclif(1,2)-((daclif(:,1)-daclif(1,1))./tand(90+fface(1,1)-fface(1,2))); 


%first order fit for cliff
[pclif,sclif]=polyfit(xc,yc,1);
[pclif1,sclif1]=polyfit(xc,yc1,1);
[pclif2,sclif2]=polyfit(xc,yc2,1);

% first order fit for platform
xp=daplat(:,1); yp=daplat(:,2);
[pplat,splat]=polyfit(xp,yp,1);

% extrapolate both
d=mean(diff(xc));
xx=min(xp):d:face(1,1); xx=xx';

[p_clif,d_clif]=polyval(pclif,xx,sclif);
[p_clif1,d_clif1]=polyval(pclif1,xx,sclif1);
[p_clif2,d_clif2]=polyval(pclif2,xx,sclif2);
[p_plat,d_plat]=polyval(pplat,xx,splat);
% locate intersection
mc=pclif(1,1); 
ic=pclif(1,2); 
mp=pplat(1,1); 
ip=pplat(1,2);
sh=[(-ic+ip)/(mc-mp),(ip*mc-ic*mp)/(mc-mp)];
shx=sh(1,1); 
shz=sh(1,2);

[shzp,shzep]=polyval(pplat,shx,splat);
% propagate 2s error by intersecting 2s regressions
% intersect max
mc=pclif1(1,1); ic=pclif1(1,2);
sh1=[(-ic+ip)/(mc-mp),(ip*mc-ic*mp)/(mc-mp)];
shx1=sh1(1,1); shz1=sh1(1,2);
[shzc1,shzec1]=polyval(pplat,shx1,splat);

% intersect min
mc=pclif2(1,1); ic=pclif2(1,2);
sh2=[(-ic+ip)/(mc-mp),(ip*mc-ic*mp)/(mc-mp)];
shx2=sh2(1,1); shz2=sh2(1,2);
[shzc2,shzec2]=polyval(pplat,shx2,splat);

% calculate error
shze1=shz1-shzec1;
shze2=shz2+shzec2;
shze=(shze2-shze1);
close(f1)
% plot intersection
f2=figure;
clf
hold on
plot(x,ma,'-b'), grid on, box on %ma
plot(x,mi,'-k'), grid on, box on %ma
plot(xp,yp,'ok','markerfacecolor','b','MarkerSize',4)
plot(xx,p_plat,'r-',xx,p_plat+2*d_plat,'r--',xx,p_plat-2*d_plat,'r--') %platform
plot(xx,p_clif,'r-')
plot(xx,p_clif1,'r--')
plot(xx,p_clif2,'r--')

plot(face(1,1),face(1,2),'sk','markerfacecolor','r')
errorbar(shx,shz,shze,'ok','markerfacecolor','k')
axis([min(da(:,1)) max(da(:,1)) min(da(:,2)) max(da(:,2))]); 
ylabel('Elevation (m)'); xlabel('Distance along profnumile (m)');

v = axis;
xpos1=v(1)+((v(2)-v(1))/6);
ypos1=v(4)-((v(4)-v(3))/4);

clear v
text(xpos1,ypos1,sprintf('Cliff slope:\n%3.1f +/- %2.1f deg.\n SA: %6.2f +/- %4.2f m',fface(1,1),fface(1,2),shz,shze))
%text(min(da(:,1))+(max(da(:,1)-min(da(:,1))/5)),max(da(:,2))-max(da(:,2))/4,sprintf('Free-face slope\n%3.1f +/- %2.1f deg.\n\nShoreline angle\n elevation\n%6.2f +/- %4.2f m',fface(1,1),fface(1,2),shz,shze))
box on
yl1=ylabel('Elevation (m)'); 
xl1=xlabel('Distance along swath(m)');
%convert to utm
point2utm=Terracem_point2utm(sh,pto,thetaB, Dprof, type_profile); %call to point2utm function
Rx=point2utm(:,1); %UTM E
Ry=point2utm(:,2); %UTM N


msgbox('SAVE FIT: Adjust zoom if necesary and push enter when ready to save')


pause  
if ishandle(f2) ==1
tit=sprintf('Station %s, Profile %u, Shoreline %u',station,profnum,shorenum);
titi=title(tit);
v = axis;
xpos=v(1)+(v(2)-v(1))/2;
ypos=v(4)-((v(4)-v(3))/14);
set(titi,'Position',[xpos ypos]);
legend('off')
%delete(tx1);
delete(yl1)%delete ylabel
delete(xl1)%elete xlabel
% The line below converts the current figure handle into a struct.
this_fig = handle2struct(gcf); 
else
msgbox('ATTENTION: Your shoreline was not saved, you must push enter to save your shoreline data')
end

%load &/or combine structures
if exist((sprintf ('%s_fits.mat',station)),'file')==2  
load(sprintf ('%s_fits.mat',station))
Fits.(sprintf('%s_%u_%u_fit',station,profnum,shorenum))=this_fig;
save((sprintf ('%s_fits.mat',station)),'Fits') 
else   
Fits.(sprintf('%s_%u_%u_fit',station,profnum,shorenum))=this_fig;
save((sprintf ('%s_fits.mat',station)),'Fits') 
end
close(findobj('type','figure','name','this_fig'))
close(findobj('type','figure','name','f2'))
close(findobj('type','figure','name','figure(2)'))
close(figure)
close(f2)
%Primary output
shoreline=[shorenum Rx Ry shx shz shze]; %number of shoreline,E, N, distance along profile,z ze
%Clicks output
%check if table exist
if exist((sprintf ('%s_Clicks.mat',station)),'file')==2 
%load table & combine structures
load(sprintf ('%s_Clicks.mat',station))
Clicks.(sprintf('prof_%u',profnum)).(sprintf('shore_%u',shorenum)).fix.clif=face;
Clicks.(sprintf('prof_%u',profnum)).(sprintf('shore_%u',shorenum)).fix.plat=plat;
%save the final combination
save((sprintf ('%s_Clicks.mat',station)),'Clicks') 
else
Clicks.(sprintf('prof_%u',profnum)).(sprintf('shore_%u',shorenum)).fix.clif=face;
Clicks.(sprintf('prof_%u',profnum)).(sprintf('shore_%u',shorenum)).fix.plat=plat;      
save((sprintf ('%s_Clicks.mat',station)),'Clicks')       
end

%% Staircase analysis function
function shoreline=Terracem_sfitprof(swat,station,profnum,shorenum,pto,thetaB,Dprof,batchmode)
%=========================================================================
%Sfitprof (staircase analysis)
%Designed by Daniel Melnick and Julius Jara-Munoz
%inputs: swath results
%=========================================================================
resp=shorenum; 
format long
%load swath profile values
x=swat(:,1); %distance along profile
mi=swat(:,2); %min values%
me=swat(:,3); %mean values
ma=swat(:,4); %max values
slp=swat(:,5); %slope values

% Check if flipped, high part must be at right
if ma(1)>ma(end)
   type_profile=2; %identifies if profile is flipped
    
   x=max(x)-x+1;   
else   type_profile=1; %identifies if profile is in correct position
 
    x=x;
end

%check batch mode
if batchmode==1
scrsz = get(0,'ScreenSize'); 
f1=figure(1);
set(f1,'Position',[scrsz(4)/2 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2])
%clf, hold on
%f=f+1;
else
    
 f1=figure(1);   
end
%plot of topographic distributions from swaths
hold on
ax1=subplot(3,1,1:2);
hold on
plot(x,ma,'k')
plot(x,me,'color',[0.3 0.3 0.3])
plot(x,mi,'color',[0.6 0.6 0.6])
hleg1 = legend('Max topography','Mean topography','Min topography','Location','SouthEast');
xlim([min(x) max(x)]);
grid on 
box on
yl1=ylabel('Elevation (m)'); 
xl1=xlabel('Distance along profile (m)');
tit=sprintf('Station %s, Profile %u, Shoreline %u',station,profnum,shorenum);
title(tit);
%plot of slope of maximum topography
ax2=subplot (3,1,3);
hold on
plot(x,slp,'.k') 
plot(x,moving(slp,50,'mean'),'r-');
grid on
box on
leyendita = legend('slope values','Location','NorthEast');
xlim([min(x) max(x)]); 
ylabel('Slope (deg)'); xlabel('Distance along profile (m)'); 
linkaxes([ax1 ax2], 'x'); 
% UNICONTROL RADIO BUTTON for select maximum or minimum topo for analysis
% Create the button group.
ht = uibuttongroup('visible','on','Position',[0 0 1 0.1],'Units','Normalized'); %radio buttons group
% Create two radio buttons inside the button group.
u1 = uicontrol('Style','radiobutton','String','Maximum swath topo','pos',[300 5 130 30],'parent',ht,'HandleVisibility','on');
u2 = uicontrol('Style','radiobutton','String','Minimum swath topo','pos',[100 5 130 30],'parent',ht,'HandleVisibility','on');
% zoom into interest area
if batchmode==1
pause(2)    
else
disp('Zoom in, select area, and press [Enter]')
%pause
end

set(ht, 'Visible','off') %hide buttons when drawing paleocliff and paleo platform
u1p=get(u1,'Value'); %detect max topo radiobutton handles
u2p=get(u2,'Value'); %%detect min topo radiobutton handles
%check which button was pushed and associate with type of topography used
if u1p==1
answer=1; %max topography
else
answer=0; %minumum topography   
end

% CHECK previous clicks routine 
% GUI backup system of results
if batchmode==1
existence=0; 
disp('enter top and bottom of paleocliff extrapolation line');
[clif(1,1),clif(1,2)]=ginputc(1,'Color','r'); 
[clif(2,1),clif(2,2)]=ginputc(1,'Color','r');
disp('enter top and bottom of paleo abrasion platform');
[plat(1,1),plat(1,2)]=ginputc(1,'Color','b'); 
[plat(2,1),plat(2,2)]=ginputc(1,'Color','b');    
else
        
if exist(sprintf ('%s_Clicks.mat',station),'file')==2
load(sprintf ('%s_Clicks.mat',station));  %assumes that Clicks empty matrix was already saved  
else
Clicks=[]; %if clicks matrix not exist creates an empty  matrix for the next condition
end

%check structure existence
doesVarExist = true;
try
    Clicks.(sprintf('prof_%u',profnum)).(sprintf('shore_%u',shorenum)).class
catch
    doesVarExist = false;
end

if doesVarExist == 1
prompt1 = {'would you like to overwrite the data? Yes(0) / No (1)'};

dlg_title1 = 'rewrite or re-read the data?';
num_lines1 = 1;

def1 = {'1'}; %default values
answer1 = str2double(inputdlg(prompt1,dlg_title1,num_lines1,def1));

        if answer1(1,1)==0; 
        existence=0; %parameter for data existence situation, existence=0 no previous data recognized.
        disp('enter top and bottom of paleocliff extrapolation line');
        [clif(1,1),clif(1,2)]=ginputc(1,'Color','r'); 
        [clif(2,1),clif(2,2)]=ginputc(1,'Color','r');
        disp('enter top and bottom of paleo abrasion platform');
        [plat(1,1),plat(1,2)]=ginputc(1,'Color','b'); 
        [plat(2,1),plat(2,2)]=ginputc(1,'Color','b');
    
        else
        existence=1;  %parameter for data existence situation, existence=1 previous data recognized.  
        clif=Clicks.(sprintf('prof_%u',profnum)).(sprintf('shore_%u',shorenum)).class.clif;
        plat=Clicks.(sprintf('prof_%u',profnum)).(sprintf('shore_%u',shorenum)).class.plat;
        end

else
% input clif
existence=0; 
disp('enter top and bottom of paleocliff extrapolation line');
[clif(1,1),clif(1,2)]=ginputc(1,'Color','r'); 
[clif(2,1),clif(2,2)]=ginputc(1,'Color','r');
disp('enter top and bottom of paleo abrasion platform');
[plat(1,1),plat(1,2)]=ginputc(1,'Color','b'); 
[plat(2,1),plat(2,2)]=ginputc(1,'Color','b');
end
end

% select cliff and plat for fitting
%choose maximum or minimum swath topography for analysis
%minimum topography is used when dune depositation is observed above the
%surface / maximum topography is used when erosion processes dominates

if answer == 1
da(:,1)=x; 
da(:,2)=ma; 
md=swat; 
ny=length(md(:,1)); 
else
da(:,1)=x; 
da(:,2)=mi; 
md=swat; 
ny=length(md(:,1));     
end

% delimitate area for interpolation
j=1; k=1;
for i=1:ny
if (da(i,1)>clif(2,1)) && (da(i,1)<clif(1,1))
    daclif(j,:)=da(i,:);
    j=j+1;
end
if (da(i,1)>plat(2,1)) && (da(i,1)<plat(1,1))
    daplat(k,:)=da(i,:);
    k=k+1;
end
end

% fit cliff
xc=daclif(:,1); yc=daclif(:,2);
[pclif,sclif]=polyfit(xc,yc,1);
[p_clif,d_clif]=polyval(pclif,xc,sclif);
% fit plat
xp=daplat(:,1); yp=daplat(:,2);
[pplat,splat]=polyfit(xp,yp,1);
[p_plat,d_plat]=polyval(pplat,xp,splat);
% save points and projected profile
clifx=[clif(1,1),clif(1,2),clif(2,1),clif(2,2)];
platx=[plat(1,1),plat(1,2),plat(2,1),plat(2,2)];

d=xc(2)-xc(1);
d=abs(d);
xx=min(xp):d:max(xc); xx=xx';
[p_clif,d_clif]=polyval(pclif,xx,sclif);
[p_plat,d_plat]=polyval(pplat,xx,splat);
close(f1)
% plot intersection
f2=figure(2);
grid off
hold on
plot(x,ma,'-b'), grid on, box on %ma
plot(x,mi,'-k'), grid on, box on %ma
plot(xc,yc,'ok','MarkerSize',5,'MarkerFaceColor','b')
plot(xx,p_clif,'r-',...
   xx,p_clif+2*d_clif,'r--',xx,p_clif-2*d_clif,'r--') %, axis equal tight
hold on
%plot(da(:,1),da(:,2),'+k')
plot(xp,yp,'ok','MarkerSize',5,'MarkerFaceColor','b')
plot(xx,p_plat,'r-',...
   xx,p_plat+2*d_plat,'r--',xx,p_plat-2*d_plat,'r--') 
axis([min(xx)-15 max(xx)+15 min(yp)-15 max(yc)+15]); 
hleg=legend('Max topography','Min topography','Platform and cliff points','linear fit','2s range','location','northwest');
yl1=ylabel('Elevation (m)'); 
xl1=xlabel('Distance along swath(m)');
% find intersections
mc=pclif(1,1); mp=pplat(1,1); ic=pclif(1,2); ip=pplat(1,2);
% output shoreline angle is sh
sh=[(-ic+ip)/(mc-mp),(ip*mc-ic*mp)/(mc-mp)];
shx=sh(1,1); %distance along ptofile
shz=sh(1,2); %elevation
[shzp,shzep]=polyval(pplat,shx,splat);
[shzc,shzec]=polyval(pclif,shx,sclif);
shze=(shzep+shzec); %error in z
shx1=shx+shze;
str1=num2str(shz,3);
str2=num2str(shze,2);
%text on plot
tx1=text(shx1-2*(max(xx)-min(xx))/3,shz+(max(yp)-min(yp)/2),sprintf('SA: %s +/- %s m',str1,str2),'FontSize',12,'BackgroundColor',[1 1 1]); %Shoreline angle value and location

%text(shx,shz,sprintf('n=%u',shz,'n=+-%u',shze,'FontSize',16)) %fix me
%obtain the UTM coordinate pairs of shoreline
point2utm=Terracem_point2utm(sh,pto,thetaB, Dprof, type_profile); %call to point2utm function
Rx=point2utm(:,1); %UTM E
Ry=point2utm(:,2); %UTM N
h=errorbar(shx,shz,shze,'ok','MarkerFaceColor','k'); %plot shoreline angle with error in z
errorbar_tick(h, 30);
box on
ylabel('Elevation (m)'); xlabel('Distance along swath (m)');
%Figs output
%figure output as a nested structure
%check if zoom adjust is first necesary
if batchmode==1
    batchmode

pause(2)      
tit=sprintf('Station %s, Profile %u, Shoreline %u',station,profnum,shorenum);
titi=title(tit);
v = axis;
xpos=v(1)+(v(2)-v(1))/2;
ypos=v(4)-((v(4)-v(3))/14);
set(titi,'Position',[xpos ypos]);
legend('off')
delete(tx1);
delete(yl1)%delete ylabel

xpos1=v(1)+((v(2)-v(1))/6);
ypos1=v(4)-((v(4)-v(3))/4);

clear v
%tx2=text(xpos1,ypos1,sprintf('SA: %u +/- %u m',shz,shze));
tx2=text(shx1-2*(max(xx)-min(xx))/3,shz+(max(yp)-min(yp)/2),sprintf('SA: %s +/- %s m',str1,str2),'FontSize',12,'BackgroundColor',[1 1 1]); %Shoreline angle value and location
delete(xl1)%delete xlabel
   
else    
msgbox('SAVE FIT: Adjust zoom if necesary and PUSH ENTER when ready')

pause

% title
tit=sprintf('Station %s, Profile %u, Shoreline %u',station,profnum,shorenum);
titi=title(tit);
v = axis;
xpos=v(1)+(v(2)-v(1))/2;
ypos=v(4)-((v(4)-v(3))/14);
set(titi,'Position',[xpos ypos]);
legend('off')
delete(tx1);
delete(yl1)%delete ylabel
delete(xl1)%delete xlabel

tx1=text(shx1-2*(max(xx)-min(xx))/3,shz+(max(yp)-min(yp)/2),sprintf('SA: %s +/- %s m',str1,str2),'FontSize',12,'BackgroundColor',[1 1 1]); %Shoreline angle value and location
end
% The line below converts the current figure handle into a struct.
this_fig = handle2struct(gcf); 
%load &/or combine structures
if exist((sprintf ('%s_fits.mat',station)),'file')==2  
load(sprintf ('%s_fits.mat',station))
Fits.(sprintf('%s_%u_%u_fit',station,profnum,shorenum))=this_fig;
save((sprintf ('%s_fits.mat',station)),'Fits') 
else   
Fits.(sprintf('%s_%u_%u_fit',station,profnum,shorenum))=this_fig;
save((sprintf ('%s_fits.mat',station)),'Fits') 
end
close(findobj('type','figure','name','this_fig'))
close(findobj('type','figure','name','f2'))
close(findobj('type','figure','name','figure(2)'))
close(figure)
close(f2)
%Primary output
shoreline=[shorenum Rx Ry shx shz shze]; %number of shoreline,E, N, distance along profile,z ze
%Clicks output
%check if table exist
if exist((sprintf ('%s_Clicks.mat',station)),'file')==2    
%load table & combine structures
load(sprintf ('%s_Clicks.mat',station))
Clicks.(sprintf('prof_%u',profnum)).(sprintf('shore_%u',shorenum)).class.clif=clif;
Clicks.(sprintf('prof_%u',profnum)).(sprintf('shore_%u',shorenum)).class.plat=plat;
%save the final combination
save((sprintf ('%s_Clicks.mat',station)),'Clicks') 
else
Clicks.(sprintf('prof_%u',profnum)).(sprintf('shore_%u',shorenum)).class.clif=clif;
Clicks.(sprintf('prof_%u',profnum)).(sprintf('shore_%u',shorenum)).class.plat=plat;      
save((sprintf ('%s_Clicks.mat',station)),'Clicks')       
end
clear nny my y_bajo bajo yp y_alto alto yc

%% Stacks analysis function
function [stacks]=Terracem_stacks(swat,station,profnum,pto,thetaB,Dprof,peak)
% Jara-Munoz (2015)
%Function designed to map shorline angles in rocky coasts, characterized by
%stacks and stumps
format long
Peack_search=(peak); %search ratio for peaks
x=swat(:,1); %distance along profile
mi=swat(:,2); %min values%
me=swat(:,3); %mean values
ma=swat(:,4); %max values
slp=swat(:,5); %slope values
% Check if flipped, high part must be at right
if ma(1)>ma(end)
   type_profile=2; %identifies if profile is flipped
    
   x=max(x)-x+1;   
else   type_profile=1; %identifies if profile is in correct position
 
    x=x;
end
%plot swath
f1=figure(1); 
clf, hold on
%plot of topographic curves
ax1=subplot(2,1,1);%new subplot
hold on
plot(x,ma,'k')
plot(x,me,'color',[0.3 0.3 0.3])
plot(x,mi,'color',[0.6 0.6 0.6])
hleg1 = legend('Max topography','Mean topography','Min topography','Location','SouthEast');
xlim([min(x) max(x)]);
box on
ylabel('Elevation (m)'); xlabel('Distance along profile (m)');
tit=sprintf('Station %s, Profile %u',station,profnum);
title(tit);
%plot of slope
ax2=subplot (2,1,2); 
hold on
plot(x,slp,'.k')
plot(x,moving(slp,50,'mean'),'r-');
leyendita = legend('slope values','Location','NorthEast');
xlim([min(x) max(x)]);
ylabel('Slope (deg)'); xlabel('Distance along profile (m)');
grid on 
box on
linkaxes([ax1 ax2], 'x'); %new
da(:,1)=x; %lenght
da(:,2)=ma; %depth
md=swat; 
ny=length(md(:,1)); 

%CHECK previous clicks routine
if exist(sprintf ('%s_Clicks.mat',station),'file')==2
load(sprintf ('%s_Clicks.mat',station));  %assumes that Clicks empty matrix was already saved  
else
Clicks=[]; %if clicks matrix not exist creates an empty  matrix for the next condition
end
doesVarExist = true;
try
    Clicks.(sprintf('prof_%u',profnum)).stack  %check existance of the structure
catch
    doesVarExist = false;
end
if doesVarExist == 1  %if the structure exist
prompt1 = {'would you like to overwrite the data? Yes(0) / No (1)'};
dlg_title1 = 'rewrite or re-read the data?';
num_lines1 = 1;
def1 = {'1'}; %default values
answer1 = str2double(inputdlg(prompt1,dlg_title1,num_lines1,def1));

        if answer1(1,1)==0; 
        existence=0; %parameter for data existence situation, existence=0 no previous data recognized.
        disp('select area for peak detection (1st upper rigth corner then lower left)');
        [quad(1,1),quad(1,2)]=ginputc(1,'color','b'); 
        [quad(2,1),quad(2,2)]=ginputc(1,'color','b');            
        else
        existence=1;  %parameter for data existence situation, existence=1 previous data recognized.          
        quad=Clicks.(sprintf('prof_%u',profnum)).stack.quad;
        
        end    
else
existence=0; 
% input peak detection area
disp('select area for peak detection (1st upper rigth corner then lower left)');
[quad(1,1),quad(1,2)]=ginputc(1,'color','b'); 
[quad(2,1),quad(2,2)]=ginputc(1,'color','b');
end
%Isolate area for peak detection
j=1; k=1;
for i=1:ny
if (da(i,1)>quad(2,1)) && (da(i,1)<quad(1,1)) %reemplazando cliff
    daquad(j,:)=da(i,:);%
    j=j+1;    
end
end
clear i j
% detect peaks
xc=daquad(:,1); 
yc=daquad(:,2);
%Peak det function (Billauer,2009)
[maxtab, mintab] = peakdet(yc, Peack_search, xc); %peak search
%plot maximum peaks detected
hold on
plot(xc,yc,'k');
plot(maxtab(:,1), maxtab(:,2),'.g','MarkerEdgeColor',[1 0 1],'MarkerSize',8);

%Clasify peaks by terrace level 
G=maxtab(:,2);
nim=numel(G);
npts=sqrt(nim);
x1=min(G):((max(G)-min(G))/npts):max(G); %this is the X axis, Rank of values in the X axis (G), the bin rank ins defined in the center, between :: 
y=G;
n_elements = histc(y,x1); %Y axe of cumulative frecuency histogram
c_elements = cumsum(n_elements); %this is for cumulative frecuency plot (not used here)
clear  n_elements c_elements h npts G y x1 nim

%Clasify peaks by terrace level 2
if existence==1
X=Clicks.(sprintf('prof_%u',profnum)).stack.XY(:,1);    
Y=Clicks.(sprintf('prof_%u',profnum)).stack.XY(:,2);  
fringes=numel(X)-1;
else
    
%Input DIALOG
pause(2) 
prompt = {'how many elevation fringes would you like to define? (maximum 3 fringes)'};
dlg_title = 'Input for peaks function';
num_lines = 1;
def = {'3'}; %default values
answer = str2double(inputdlg(prompt,dlg_title,num_lines,def));
fringes=answer(1,1);
disp('mark only the upper bound of each fringe')   
[X,Y]=ginputc(fringes,'color','g'); 
Y=vertcat(Y,zeros(1,1)); %add zeros file (for loop) 
X=vertcat(X,zeros(1,1));
end

maxi=max(x);
% clasify the matrix maxtab in the fringe ranges 
hold off
nh=length(maxtab(:,1)); 
for i=1:(fringes) 
    
    for h=1:nh    
    if (maxtab(h,2)<Y(i,1)) && (maxtab(h,2)>Y(i+1,1))    
        maxtab(h,3)=i;     
    end
    end
    
    line([0 maxi],[Y(i,:) Y(i,:)], 'Color', [0 1-(i/(fringes)) 0],'LineStyle','--')   
end
clear i h j maxi mintab md scrsz tit daquad def nh numlines prompt ax1 ax2 k nh ny yc xc f da 
% shoreline angle on terrace remains on the cliffs 
if existence==1
            if isempty(Clicks.(sprintf('prof_%u',profnum)).stack.XSYS)==1
            shore=[];
            XS=[];
            YS=[];
            else              
            XS=Clicks.(sprintf('prof_%u',profnum)).stack.XSYS(:,1); 
            YS=Clicks.(sprintf('prof_%u',profnum)).stack.XSYS(:,2);                
            shore(:,1)=XS;
            shore(:,2)=YS;
            ns=numel(YS); 
            end
else       
% mark the shoreline point in the screen according with the Fringe
disp ('mark shoreline angles, if not push ENTER')
[XS,YS]=ginputc(fringes,'color','k'); 
            if isempty(YS)==1
            shore=[];
            else
            shore(:,1)=XS;
            shore(:,2)=YS;
            ns=numel(YS); 
            end
end

if size(shore)==0
    maxtab=maxtab;
else   
% loop to automaticaly identify elevation fringe for each shoreline angle
%maximum number is 3 fringes
    for h=1:ns %1:number of shoreline angles selcted could be different that number of fringes in some cases    
        for i=1:fringes %number of fringes   
         if (shore(h,2)<Y(i,1)) && (shore(h,2)>Y(i+1,1))    
        shore(h,3)=i;                
        % Condition: if two shorelines angles have the same fringe range
        % move the minor one to next fringe
        %new condition to check for repeated fringe values in shorelines   
        numb=sum(shore(:,3)==i);             
         if numb==1 || numb==0            
          %shore(h,3)=i; 
          intplat(i).shore=shore(h,:);            
         else   % Move repeated to lower level fringe    
            sortshore=sort(shore,2,'descend');            
            %is not possible to use arbitrary asinged fringes
            %so, we detect repeated fringes to evaluate the next alternatives            
            tabl=tabulate(sortshore(:,3));
            c=max(tabl); %statistics of sortshore
            rept=c(:,1); %identified repeated value              
            %Alternative1 if repeated is 2 use this  (working good)
               if rept==2                   
            intplat(i).shore=sortshore(1,:);
            intplat(i).shore=sortshore(2,:);
            j=i+1;   
            intplat(j).shore=sortshore(3,:);
            intplat(j).shore(:,3)=j;       
               end                 
            %Alternative2 if repeated is 1 use this (working good)
               if rept==1         
            intplat(i).shore=sortshore(1,:);
            k2=i+1;    
            intplat(k2).shore=sortshore(2,:);
            intplat(k2).shore(:,3)=k2;                                      
               end          
         end   
         end  
        end 
    end   
clear i h 
hold on
plot(shore(:,1),shore(:,2),'.k','MarkerSize',10);
% create outputs table
% add shoreline angle data to peaks point matrix
maxtab=vertcat(intplat(1,:).shore,maxtab);
end
% clasify and dismember maxtab on each elevation fringe
for i = 1:fringes
plat(i).level = maxtab(maxtab(:,3) == i,:);
end
clear i sortshore rept c tabl numb shore ns 
% separated vectors and fit stacks in profile
% plot and function interplations
close(f1)
f2=figure(2);
hold on  
for i=1:fringes
   xpeak(:,1)=plat(1,i).level(:,1);
   ypeak(:,1)=plat(1,i).level(:,2);   
   % stack 
   [psurf,ssurf]=polyfit(xpeak,ypeak,1);
   dxm=1; %mean(diff(xpeak));%prepare to extend fit
   xxx=min(xpeak):dxm:max(x);
   xxx=xxx'; %range added to extend fit line   
   [p_surf,d_surf]=polyval(psurf,xxx,ssurf);
% save fiting and other data to table structure
   intplat(i).slope=psurf; 
   intplat(i).p_surf=p_surf;
   intplat(i).d_surf=d_surf;
   intplat(i).xpeak=xpeak(:,1); %x value for peaks
   intplat(i).ypeak=ypeak(:,1); %y value for peaks
   intplat(i).xxx=xxx; %x values for interpolation      
%extract fit-maxtopo intersection
m_rect=intplat(i).slope(1,1); %slope of interpolation rect (tand of rect angle)
%c value
c_rect=intplat(i).slope(1,2); %second coeficient of the rect ecuation 'c'
% check this part because replace x and ma por swat values
swati(:,1)=x;
swati(:,2)=ma;
for j=1:(numel(swati(:,1)));
if fix(m_rect*(swati(j,1))+(c_rect)) == nearest(swati(j,2)); %aproximation to topography    
    shore_int(:,1)=swati(j,1);
    shore_int(:,2)=swati(j,2);   
end
end 
   %take error from interpolation and store in the third row
   %take number of shorelines from the number of fringes in loop
   shore_int(:,3)=min(intplat(i).d_surf); %erro in Y
   shore_int(:,4)=i; %shorenum   
   %store level of each shoreline angle in the fourth row  
   %create structure table
   intplat(i).shore_int=shore_int; %projected shoreline point 
   %FIGURE FINAL    
   plot(x,ma,'-k'), 
   plot(x,me,'Color', [0.5 0.5 0.5]), 
   plot(x,mi,'Color', [0.7 0.7 0.7]), 
   xlim([min(x) max(x)]); 
   ylim([0 max(ma)])
   box on; 
   ylabel('Elevation (m)'); xlabel('Distance along profile (m)');
   plot(intplat(i).xxx,intplat(i).p_surf,'Color', [0+(i/(fringes)) 1-(i/(fringes)) 0],'LineStyle','--')
   plot(intplat(i).xpeak,intplat(i).ypeak,'.','Color', [0+(i/(fringes)) 1-(i/(fringes)) 0],'MarkerSize',20)
   errorbar(intplat(i).shore_int(:,1),intplat(i).shore_int(:,2),intplat(i).shore_int(:,3),'dk','markerfacecolor','k')   
   clear xpeak ypeak xxx dxm p_surf d_surf psurf ssurf shore_int     
end  
   %INSET
   G=maxtab(:,2);  
   nim=numel(G);
   npts=sqrt(nim);
   x1=min(G):((max(G)-min(G))/npts):max(G); %this is the X axis, Rank of values in the X axis (G), the bin rank ins defined in the center, between :: 
   y=G;
   n_elements = histc(y,x1); %Y axe of cumulative frecuency histogram
   c_elements = cumsum(n_elements);
    
ax3=axes('position',[0.15 0.7 0.15 0.2]);
   %%histogram of height
    hist((G),(npts));
    ylim([0 nim/1.5])    
    h = findobj(gca,'Type','patch');
    set(h,'EdgeColor',[0.9 0.9 0.9])
    set(gca, 'YTick', []);
    set(gca,'yaxislocation','right','tickdir','in')
    text((nim-(nim/8)),npts+(npts/1),sprintf('n=%u',nim))
    xlabel('Elevation (m)'); ylabel('Peaks frecuency')
    set(gcf,'PaperUnits','centimeters'); set(gcf, 'PaperType','A4');
    fout = sprintf('%s_%u_stack.pdf',station,profnum);
    saveas(gcf,fout,'pdf');
% outputs
shoreline_angles=vertcat(intplat(1,:).shore_int); %main output shorelineangles with error 
x_sh=shoreline_angles(:,1); %distance along profile
y_sh=shoreline_angles(:,2); %elevation
sh1(:,1)=x_sh; 
sh1(:,2)=y_sh;
e_sh=shoreline_angles(:,3);
num_sh=shoreline_angles(:,4);
all_peaks1=maxtab;
%obtaining the UTM coordinate pairs of shoreline angle
nimb=numel(shoreline_angles(:,4));
for q=1:nimb  
point2utm(q,:)=Terracem_point2utm(sh1(q,:),pto,thetaB, Dprof, type_profile); %call to point2utm function
Rx(q,1)=point2utm(q,1); %UTM E
Ry(q,1)=point2utm(q,2); %UTM N
end
%obtaining the UTM coordinate pairs for stack points (maxtab)
sh2(:,1)=all_peaks1(:,1);%peak distance along profile
sh2(:,2)=all_peaks1(:,2);%peak elevation
sh2(:,3)=all_peaks1(:,3);%peak code
nimc=numel(all_peaks1(:,1)); %number of peaks
for r=1:nimc   
point2utm(r,:)=Terracem_point2utm(sh2(r,:),pto,thetaB, Dprof, type_profile); %call to point2utm function
end
%Primary output table
stacks=[num_sh Rx Ry x_sh y_sh e_sh]; 
%Figures as structure output
msgbox('SAVE FIT: Adjust zoom and push ENTER when ready')
pause
% The line below converts the current figure handle into a struct.
this_fig = handle2struct(gcf); 
%load &/or combine structures
if exist((sprintf ('%s_fits.mat',station)),'file')==2  
load(sprintf ('%s_fits.mat',station))
Fits.(sprintf('%s_%u_1_fit',station,profnum))=this_fig;
save((sprintf ('%s_fits.mat',station)),'Fits') 
else   
Fits.(sprintf('%s_%u_1_fit',station,profnum))=this_fig;
save((sprintf ('%s_fits.mat',station)),'Fits') 
end

close(f2)
%Clicks data output
%check if table exist
if exist((sprintf ('%s_Clicks.mat',station)),'file')==2    
%load table & combine structures
load(sprintf ('%s_Clicks.mat',station))
Clicks.(sprintf('prof_%u',profnum)).stack.quad=quad;
Clicks.(sprintf('prof_%u',profnum)).stack.XY=[X,Y];
Clicks.(sprintf('prof_%u',profnum)).stack.XSYS=[XS,YS];
%save the final combination
save((sprintf ('%s_Clicks.mat',station)),'Clicks') 
else    
Clicks.(sprintf('prof_%u',profnum)).stack.quad=quad;
Clicks.(sprintf('prof_%u',profnum)).stack.XY=[X,Y];
Clicks.(sprintf('prof_%u',profnum)).stack.XSYS=[XS,YS];    
%if table not exist, just save    
save((sprintf ('%s_Clicks.mat',station)),'Clicks') 
end

%% Diffussion analysis function
function [shoreline Difusion]=Terracem_terracediff(swat,station,profnum,shorenum,pto,thetaB,Dprof,thetadegrees,k,t1,ti,t2)
%terrace difussion script By Julius J. and D. Melnick
%inputs: Terracem swat formatted input.
resp=shorenum; %A) modify shorenum to input list
format long
x=swat(:,1); %distance along profile
mi=swat(:,2); %min values%
me=swat(:,3); %mean values
ma=swat(:,4); %max values
slp=swat(:,5); %slope values
% Check if flipped, high part must be at right
if ma(1)>ma(end)
   type_profile=2; %identifies if profile is flipped    
   x=max(x)-x+1;   
else   type_profile=1; %identifies if profile is in correct position 
    x=x;
end

figure(1)
clf, hold on
ax1=subplot(3,1,1:2);
hold on
plot(x,ma,'k')
plot(x,me,'color',[0.3 0.3 0.3])
plot(x,mi,'color',[0.6 0.6 0.6])
hleg1 = legend('Max topography','Mean topography','Min topography','Location','SouthEast');
xlim([min(x) max(x)]);
ylabel('Elevation (m)'); 
xlabel('Distance along profile (m)');
tit=sprintf('Station %s, Profile %u, Shoreline %u',station,profnum,shorenum);
title(tit);
box on
ax2=subplot (3,1,3); 
hold on
plot(x,slp,'.k')
plot(x,moving(slp,50,'mean'),'r-');
leyendita = legend('slope values','Location','NorthEast');
xlim([min(x) max(x)]);
ylabel('Slope (deg)'); xlabel('Distance along profile (m)');
grid on 
box on
linkaxes([ax1 ax2], 'x');
%disp('Zoom in, select area, and press [Enter]')
%pause
if exist(sprintf ('%s_Clicks.mat',station),'file')==2
load(sprintf ('%s_Clicks.mat',station));  %assumes that Clicks empty matrix was already saved  
else
Clicks=[]; %if clicks matrix not exist creates an empty  matrix for the next condition
end
%check structure existence
doesVarExist = true;
try
    Clicks.(sprintf('prof_%u',profnum)).(sprintf('shore_%u',shorenum)).diff
catch
    doesVarExist = false;
end

if doesVarExist == 1
prompt1 = {'would you like to overwrite the data? Yes(0) / No (1)'};

dlg_title1 = 'rewrite or re-read the data?';
num_lines1 = 1;

def1 = {'1'}; %default values
answer1 = str2double(inputdlg(prompt1,dlg_title1,num_lines1,def1));

        if answer1(1,1)==0; 
        existence=0; %parameter for data existence situation, existence=0 no previous data recognized.
        %top point
        disp('...enter inflection point of far-field slope at top surface...');
        [tops(1,1),tops(1,2)]=ginputc(1,'Color','r');       
        %base points
        disp('...enter platform surface...');
        [bots(1,1),bots(1,2)]=ginputc(1,'Color','b'); 
        [bots(2,1),bots(2,2)]=ginputc(1,'Color','b');     
        else
        existence=1;  %parameter for data existence situation, existence=1 previous data recognized.  
        tops=Clicks.(sprintf('prof_%u',profnum)).(sprintf('shore_%u',shorenum)).diff.tops;
        bots=Clicks.(sprintf('prof_%u',profnum)).(sprintf('shore_%u',shorenum)).diff.bots;
       
        end
else
% MANUAL INPUTS
%top point
existence=0; 
disp('...enter inflection point of far-field slope at top surface...');
[tops(1,1),tops(1,2)]=ginputc(1,'Color','r');
%base points
disp('...enter bottom surface...');
[bots(1,1),bots(1,2)]=ginputc(1,'Color','b');
[bots(2,1),bots(2,2)]=ginputc(1,'Color','b');
end
 plot(tops(1,1),tops(1,2),'ob','markerfacecolor','k');
 plot(bots(1,1),bots(1,2),'ob','markerfacecolor','b');
 plot(bots(2,1),bots(2,2),'ob','markerfacecolor','b'); 
% new part for center point
center=bots(1,2)+((tops(1,2)-bots(1,2))/2);
plot([min(x) max(x)],[center center],'k--')

if existence==1
scarp=Clicks.(sprintf('prof_%u',profnum)).(sprintf('shore_%u',shorenum)).diff.scarp;
else
% center point
disp('...enter scarp center...');
[scarp(1,1),scarp(1,2)]=ginputc(1,'Color','g'); 
end

plot(scarp(1,1),scarp(1,2),'sb','markerfacecolor','g');
% fitting line at base with maximum topo
da(:,1)=x; %x axis now defined as da(:,1)
da(:,2)=ma; %fit respect to maximum topo values
da=sortrows(da,1); 
ny=length(da(:,1)); 
daplat1=da(da(:,1)<bots(1,1),:); %only at base
daplat=daplat1(daplat1(:,1)>bots(2,1),:);
clear daplat1
% linear regressions
% platform
xp=daplat(:,1); 
yp=daplat(:,2);
[pplat,splat]=polyfit(xp,yp,1);
% extrapolate platform
d=mean(diff(xp));
xx=min(xp):d:max(xp); xx=xx';
[p_plat,d_plat]=polyval(pplat,xx,splat);
%find x and y values
nr=numel(p_plat);
x1=xx(1,1); 
x2=xx(nr,1); 
y1=p_plat(1,1); 
y2=p_plat(nr,1); %rect eq.
%central point
x_scarp=scarp(1,1);
y_scarp=(((x_scarp-x1)*(y2-y1))/(x2-x1))+y1;
%top points(tops)
x_tops=tops(1,1);
y_tops=(((x_tops-x1)*(y2-y1))/(x2-x1))+y1;
%min value of x in profile
x_minx=min(da(:,1));
y_minx=(((x_minx-x1)*(y2-y1))/(x2-x1))+y1;
%table with basal line elements
bline(1,1)=x_scarp;
bline(1,2)=y_scarp;
bline(2,1)=x1;
bline(2,2)=y1;
bline(3,1)=x_minx;
bline(3,2)=y_minx;
%slope of far field lines %%%%% (b value comes from here)
ml=(y2-y1)/(x2-x1);
%line at top: with slope m1 passing by tops and x_max
%ec en la forma y_max-yt=m(x_max-xt) donde y_max es incognita
xt=tops(1,1); yt=tops(1,2);
x_max=max(da(:,1)); %topline at max x value
y_max=(ml*(x_max-xt))+yt;
x_scarp=scarp(1,1);
yt_scarp=(ml*(x_scarp-xt))+yt;
%table with topline elements
tline(1,1)=x_scarp;
tline(1,2)=yt_scarp;
tline(2,1)=xt;
tline(2,2)=yt;
tline(3,1)=x_max;
tline(3,2)=y_max;
%central vertical line that pass trough scarp 
%values of scarp in y at top and baseline 
v_line(1,1)=x_scarp;
v_line(1,2)=yt_scarp;
v_line(2,1)=x_scarp;
v_line(2,2)=y_scarp;
% precentered data plot %%%%%%%%%%%%%%%%%%
hold on
plot(bline(:,1),bline(:,2),'r--');
plot(tline(:,1),tline(:,2),'r--');
plot(v_line(:,1),v_line(:,2),'k','LineWidth',1);
text((v_line(1,1)),((v_line(1,2)+scarp(1,2))/2),'a','FontSize',9,'BackgroundColor','white')
text((v_line(1,1)),((v_line(2,2)+scarp(1,2))/2),'a','FontSize',9,'BackgroundColor','white')
%values of a height of scarp/2
a=(yt_scarp-y_scarp)/2;
%value of b far field slope
b=atand(ml);
text(bots(2,1),tops(1,2),['2a=',(num2str(2*a,'%0.1f')),' m'],'FontSize',12)
text(bots(2,1),((tops(1,2)+scarp(1,2))/2),['b=',(num2str(b,'%0.1f')),' ?'],'FontSize',12)
%save plot
rect = [5,6,13,18];
set(gcf,'PaperUnits','centimeters'); set(gcf, 'PaperType','A4'); set(gcf,'paperposition',rect);
fout = sprintf('%s_prof_%u_shore_%u_diffuse_picks.pdf',station,profnum,shorenum);
saveas(gcf,fout,'pdf');

% CENTERING THE ENTIRE PLOT DATA
%defining plot in cartesian plane
offx=x_scarp; %centerpoint x
offy=yt_scarp-a; %centerpoint y
%centering topografic data
nim=numel(da(:,1));
xn(:,1)=x(:,1)-offx;    
dan(:,1)=da(:,1)-offx; %for max topo analysis
dan(:,2)=da(:,2)-offy;
%centering far field slope data
num=numel(tline(:,1)); %numel for bline and tline
tlinen(:,1)=tline(:,1)-offx;
tlinen(:,2)=tline(:,2)-offy;
blinen(:,1)=bline(:,1)-offx;
blinen(:,2)=bline(:,2)-offy;
vlinen(:,1)=v_line(:,1)-offx;
vlinen(:,2)=v_line(:,2)-offy;
%centering bots and tops points of ginputs
topsn(1,1)=tops(1,1)-offx;
topsn(1,2)=tops(1,2)-offy;
botsn(1,1)=bots(2,1)-offx;
botsn(1,2)=bots(2,2)-offy;
scarpn(1,1)=scarp(1,1)-offx;
scarpn(1,2)=scarp(1,2)-offy;
%generating database "parameters"
parameters.offset.uperline=tlinen;
parameters.offset.lowerline=blinen;
parameters.offset.vline=vlinen;
parameters.profile.x=dan(:,1);
parameters.profile.maxtopo=dan;
parameters.profile.slope(:,1)=xn;
parameters.profile.slope(:,2)=slp;
parameters.inputs.a=a;
parameters.inputs.b=b;
parameters.profile.nim=nim; %nim=number of elements for dx calculation
%add ginputs point to paramaters to constrain the extension of rms analysis
parameters.inputs.topsn=topsn;
parameters.inputs.botsn=botsn;
parameters.inputs.scarp=scarpn;
%cut plot area using botsn and topsn
parameters.cut.maxtopo=parameters.profile.maxtopo; %vectorial form: first create a matrix equal that the one that i want to modify
parameters.cut.maxtopo(parameters.cut.maxtopo(:,1)>(parameters.inputs.topsn(1,1)+15),:)=[]; %vectorial form, excluve values bigger than x in topsn
parameters.cut.maxtopo(parameters.cut.maxtopo(:,1)<(parameters.inputs.botsn(1,1)-0),:)=[]; %vectorial form, excluve values bigger than x in botsn
nm=numel(parameters.cut.maxtopo(:,1));
parameters.cut.nim=nm;
clear x
%generate a sintetic distribution of the same size of the input topography
xmin=(parameters.cut.maxtopo(1,1));
dx=((parameters.cut.maxtopo(parameters.cut.nim,1)-parameters.cut.maxtopo(1,1))/(parameters.cut.nim-1));
xmax=(parameters.cut.maxtopo(parameters.cut.nim,1));
x=(xmin:dx:xmax);
a=parameters.inputs.a; %half oh the scarp elevation obtained from manual picking at the begining
bdegrees=parameters.inputs.b; %far field slope value, obtained from manual picking at the begining
%difussion loop
thetadegrees
k
t=t1:ti:t2; %range in time
ktrange=t.*k;
rms=zeros(numel(ktrange),1); %matris de zeros igual a ktrange
topopred=(parameters.cut.maxtopo(:,2))';

% RUN MODELS
j=1;
for kt=ktrange
model=finitescarp(thetadegrees,a,x,kt,bdegrees);
% RMS
rms(j)=sqrt(sum((topopred(:)-(-model(:))).^2)/numel(topopred));
j=j+1;
end
%plot RMS
rr=horzcat(ktrange',rms);
ir=sortrows(rr,2);
parameters.KT_RMS=ir;
best=parameters.KT_RMS(1,:)
% Shoreline angle extraction and error %%%%%%%
fface=[thetadegrees 5]; %angular rangle, mean=60 error =5
da2=parameters.cut.maxtopo;
daxclif=da2(da2(:,1)>=0,:); %change scarp por la inversa del shoreline angle para calcular variacion angular
xc=-daxclif(:,1); 
yc=(-tand(fface(1,1)).*abs(daxclif(:,1)-daxclif(1,1)))+daxclif(1,2);
yc1=daxclif(1,2)-((daxclif(:,1)-daxclif(1,1))./tand(90-fface(1,1)+fface(1,2)));
yc2=daxclif(1,2)-((daxclif(:,1)-daxclif(1,1))./tand(90-fface(1,1)-fface(1,2)));    
[pclif,sclif]=polyfit(xc,yc,1);
[pclif1,sclif1]=polyfit(xc,yc1,1);
[pclif2,sclif2]=polyfit(xc,yc2,1);
%platform
xg=daplat(:,1)-offx; %valores alcanzados por la interpolacion
yg=daplat(:,2)-offy;
[pplat,ssplat]=polyfit(xg,yg,1);
% extrapolate
dxm=mean(diff(xg));
%xxx=min(xg):d:(xxa+offx); 
xxx=min(xg):dxm:max(xc); 
xxx=xxx';
[pp_plat,dd_plat]=polyval(pplat,xxx,ssplat);
[p_clif,d_clif]=polyval(pclif,xxx,sclif);
% intersect
mc=pclif(1,1); 
ic=pclif(1,2); 
mp=pplat(1,1); 
ip=pplat(1,2);
sh=[(-ic+ip)/(mc-mp),(ip*mc-ic*mp)/(mc-mp)];
eshx=sh(1,1); 
eshz=sh(1,2);
[shzp,shzep]=polyval(pplat,eshx,ssplat);
% propagate 2s error by intersecting 2s regressions
% intersect max
mc=pclif1(1,1); 
ic=pclif1(1,2);
sh1=[(-ic+ip)/(mc-mp),(ip*mc-ic*mp)/(mc-mp)];
shx1=sh1(1,1); 
shz1=sh1(1,2);
[shzc1,shzec1]=polyval(pplat,shx1,splat);
% intersect min
mc=pclif2(1,1); 
ic=pclif2(1,2);
sh2=[(-ic+ip)/(mc-mp),(ip*mc-ic*mp)/(mc-mp)];
shx2=sh2(1,1); 
shz2=sh2(1,2);
[shzc2,shzec2]=polyval(pplat,shx2,splat);
% calculate error
shze1=shz1-shzec1;
shze2=shz2+shzec2;
shze=(shze2-shze1);
%transform xc and yc to original coordinates
xcyc(:,1)=xc; xcyc(:,2)=yc;
xcyc=xcyc(xcyc(:,2)>eshz,:);
xcyc1(:,1)=xc; xcyc1(:,2)=yc1;
xcyc1=xcyc1(xcyc1(:,2)>eshz,:);
xcyc2(:,1)=xc; xcyc2(:,2)=yc2;
xcyc2=xcyc2(xcyc2(:,2)>eshz,:);
%outputs
shoreline_angle_x=eshx+offx; %distance along profile
shoreline_angle_y=eshz+offy; %elevation
error_shoreline=2*shze; %error in elevation
%plot results
f2=figure(2);
clf
hinitial=finitescarp(thetadegrees,a, x, 0.0001,bdegrees);
bestfit=finitescarp(thetadegrees,a,x,best(1,1),bdegrees);
tit1=sprintf('Station %s, Profile %u, Shoreline %u',station,profnum,shorenum);
ax1=axes('position',[0.2 0.1 0.63 0.65]);
title(tit1)
hold on
box on
plot(x,-hinitial,'b--')
plot(parameters.cut.maxtopo(:,1),parameters.cut.maxtopo(:,2),'k','LineWidth',1.5)
plot(x,-bestfit,'r-')
hleg4 = legend('Initial condition','Topography','Best fit','location','SouthEast');
set(gca,'yaxislocation','left','tickdir','in')
set(gca,'xaxislocation','bottom','tickdir','in');
xlabel('Distance (m)')
ylabel('Elevation (m)')
xlim([(parameters.cut.maxtopo(1,1)) (parameters.cut.maxtopo(parameters.cut.nim,1))])
plot(eshx,eshz,'ok')
errorbar(eshx,eshz,shze,'dk','markerfacecolor','k')
% miniplot rms      
ax2=axes('position',[0.275 0.48 0.17 0.15]);
hold on
set(gca,'xaxislocation','bottom','tickdir','in');
plot(ktrange,rms,'b')
plot(best(1,1),best(1,2),'ok','markerfacecolor','k')
line([0 best(1,1)],[best(1,2) best(1,2)],'Linestyle','--','color','k');
line([best(1,1) best(1,1)],[0 best(1,2)],'Linestyle','--','color','k');
best_kt=best(:,1);
best_rms=best(:,2);
ak=round(best(1,1)+(best(1,1)/1));
ck=round(best(1,1)-(best(1,1)/1));
bk=round((best(1,2))+(best(1,2)/1));
lk=round(best(1,2)-(best(1,2)/1));
ylim([0 (best(1,2)*2)])
%xlim('auto')  
xlim([0 (best(1,1)*2)])
set(gca,'XTick',best(1,1))
set(gca,'YTick',best(1,2))
tit2=title('RMS misfit');
xlabel('Kt (m^2)')
ylabel('RMS (m)')
%sizing
v = axis;
xpos=v(1)+(v(2)-v(1))/2;
ypos=v(4)-((v(4)-v(3))/10);
set(tit2,'Position',[xpos ypos]);
%save figure as pdf
rect = [2,8,17,20];
set(gcf,'PaperUnits','centimeters'); set(gcf, 'PaperType','A4'); set(gcf,'paperposition',rect);
fout = sprintf('%s_prof_%u_shoreline_%u_diff.pdf',station,profnum,shorenum);
saveas(gcf,fout,'pdf');
clear p_clif p_plat pclif1 pclif2 pp_plat 
%obtaining the UTM coordinate pairs of shoreline
p1(:,1)=shoreline_angle_x; %formatting iputs for point2utm
p1(:,2)=shoreline_angle_y; %formatting iputs for point2utm
point2utm=Terracem_point2utm(p1,pto,thetaB, Dprof, type_profile); %call to point2utm function
Rx=point2utm(:,1); %UTM E
Ry=point2utm(:,2); %UTM N
%Primary output
shoreline=[shorenum Rx Ry shoreline_angle_x shoreline_angle_y shze];
%secondary outputs
Difusion=horzcat(profnum,shorenum, Rx, Ry, shoreline_angle_y,shze,best_kt, best_rms, a, b);
%Clicks output
%check if table exist
if exist((sprintf ('%s_Clicks.mat',station)),'file')==2
%load table & combine structures
load(sprintf ('%s_Clicks.mat',station))
Clicks.(sprintf('prof_%u',profnum)).(sprintf('shore_%u',shorenum)).diff.tops=tops;
Clicks.(sprintf('prof_%u',profnum)).(sprintf('shore_%u',shorenum)).diff.bots=bots;
Clicks.(sprintf('prof_%u',profnum)).(sprintf('shore_%u',shorenum)).diff.scarp=scarp;
save((sprintf ('%s_Clicks.mat',station)),'Clicks') 
else
Clicks.(sprintf('prof_%u',profnum)).(sprintf('shore_%u',shorenum)).diff.tops=tops;
Clicks.(sprintf('prof_%u',profnum)).(sprintf('shore_%u',shorenum)).diff.bots=bots;
Clicks.(sprintf('prof_%u',profnum)).(sprintf('shore_%u',shorenum)).diff.scarp=scarp;
save((sprintf ('%s_Clicks.mat',station)),'Clicks') 
end


%% Nested functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Function 1: XYZ to ASCII
function xyz=ascii2xyz(varargin)
%ASCII2XYZ- convert ARC ASCII text file to xyz
%
%   ASCII2XYZ reads in a raster text file in 
%       ARC ASCII format and converts values to 
%       a matrix of x, y, a z values. No data 
%       values are ommited from the output matrix.
% Andrew Stevens, 10/16/2008
% astevens@usgs.gov

%process optional input arg.
if nargin>0;
    fname=varargin{1};
    if exist(fname,'file')==0
        error('File not found, check the file name and try again.');
    end
else
    [filename, pathname] = uigetfile( ...
        {'*.txt','TXT Files (*.txt)'; ...
        '*.asc','ASC Files (*.asc)';...
        '*.*',  'All Files (*.*)'}, ...
        'Pick a file');
    fname=fullfile(pathname,filename);
end

fid=fopen(fname,'r');
%read header
format='%s %f';
hdr=cell(6,2);
try
    for i=1:6;
        [hdr(i,1),hdr{i,2}]=...
            strread(fgetl(fid),format);
    end    
    %try to determine from header if the file uses
    %the corner of the grid or the pixel center. If 
    if findstr(hdr{3,1},'center')~=0
        offset=0;     
    elseif findstr(hdr{3,1},'corner')~=0
        offset=hdr{5,2}/2;
    else %in case the header is poorly formatted
        ansr=questdlg({'Is the grid referenced to the pixel center';...
            'or to the corner of the grid?'}, ...
            'Spatial Reference?', ...
            'Center', 'Corner', 'Center');
        if strcmpi(ansr,'center')
            offset=0;
        else
            offset=hdr{5,2}/2;
        end
    end        
catch %#ok
    fclose(fid)
    errorstr.identifier='myToolbox:ascii2xyz:badHeader';
    errorstr.message=['Error reading header information.',...
        ' See required format in the help section.']
    error(errorstr)
end
%vector of x, y, and z positions
xv=(hdr{3,2}:hdr{5,2}:hdr{3,2}+...
    ((hdr{1,2}-1)*hdr{5,2}))+offset;
yv=(fliplr(hdr{4,2}:hdr{5,2}:hdr{4,2}+...
    ((hdr{2,2}-1)*hdr{5,2})))+offset;
xvec=repmat(xv,[1 hdr{2,2}])';
yvec=cell2mat(cellfun(@(x)(repmat(x,[hdr{1,2} 1])),...
num2cell(yv),'uni',0));
%read data
zvec=fscanf(fid,'%f',hdr{1,2}*hdr{2,2});
fclose(fid);
%get rid of values with no data
xyz=[xvec(zvec~=hdr{6,2}),...
yvec(zvec~=hdr{6,2}),...
zvec(zvec~=hdr{6,2})];

%% Function 2: myTable
function myTable(file,cnames) 
%Inputs:
%table shorelines in each stations directory
%stations name

h = figure('Position',[600 400 902 600],'numbertitle','off','MenuBar','none','ToolBar','none');
%Data input, read shorelines table
%defaultData = rand(90,10); %test default
%cd(handles.stationsdir)
%file=('shorelines.txt');
if exist(file,'file')==1 || exist(file,'file')==2    
%load shorelines table
fid=fopen(file,'r');
i=1;
while ~feof(fid)
    %tline = fgetl(fid);
    A=textscan(fid,'%s %u %u %f %f %f %f %f %s %f/','delimiter',' ','HeaderLines',1);
    %works with this weird / 
    %check headerlines
    
    for j=1:10 %number of columns (default)
        data(i,j) = A(1,j);
    end
    i=i+1;
end

fclose(fid);
data=sortrows(data,10);
nim=numel(data(:,10));
u=1;

for k=1:nim
t1(u,:)=cellstr(data{k,1}); %station
t2(u,:)=(data(k,2));
t3(u,:)=data(k,3);
t4(u,:)=data(k,4);
t5(u,:)=data(k,5);
t6(u,:)=(data(k,6));
t7(u,:)=data(k,7);
t8(u,:)=data(k,8);
t9(u,:)=cellstr(data{k,9}); % type analysis
t10(u,:)=(data(k,10));
u=1+u;
end

defaultData=horzcat(t1,t2,t3,t4,t5,t6,t7,t8,t9,t10);

% Digital table creation
%cnames = {'Station', 'Profile_number', 'shoreline_number', 'east', 'north', 'distance_along_swath', 'shoreline_elevation', 'error', 'analysis_type', 'time'};%colun names

uitable(h,'Units','normalized','Position',[0 0 1 1],...
              'Data', defaultData,... 
              'Tag','myTable',...    
              'ColumnName', cnames,'RowName',[],...
              'CellSelectionCallback',@cellSelect);

% create pushbutton to delete selected rows
tb = uitoolbar(h);

%icon colors for delete button
a = [.20:.05:0.95];
b(:,:,1) = repmat(a,16,1)';
b(:,:,2) = repmat(a,16,1);
b(:,:,3) = repmat(flipdim(a,2),16,1);

uipushtool(tb,'CData',b,'tooltipstring','delete','ClickedCallback',@deleteRow);


end

%% Function 3: select cell
function cellSelect(src,evt)
% get indices of selected rows and make them available for other callbacks
index = evt.Indices;
if any(index)             %loop necessary to surpress unimportant errors.
    rows = index(:,1);
    set(src,'UserData',rows);
end

%% Function 4: celldelete row
function deleteRow(~,~)
th = findobj('Tag','myTable');
% get current data
data1 = get(th,'Data');
% get indices of selected rows
rows = get(th,'UserData');
% create mask containing rows to keep
mask = (1:size(data1,1))';
mask(rows) = [];
% delete selected rows and re-write data
data1 = data1(mask,:);
set(th,'Data',data1);

%save data1 to .txt table


  dlmcell('shorelines.txt',data1,'delimiter',' '); 


%% Function 5: dlmcell
function dlmcell(file,cell_array,varargin) 
% dlmcell - Write Cell Array to Text File     
%                                                 Version:    01.06.2010 
%                                                     (c) Roland Pfister 
%                                             roland_pfister@t-online.de 
%                        ...with many thanks to George Papazafeiropoulos 
%                        for his corrections and improvements.           
% 1. Synopsis                                                            
%                                                                        
% A single cell array is written to an output file. Cells may consist of 
% any combination of (a) numbers, (b) letters, or (c) words. The inputs  
% are as follows:                                                        
%                                                                        
%       - file       The output filename (string).                       
%       - cell_array The cell array to be written.                       
%       - delimiter  Delimiter symbol, e.g. ',' (optional;               
%                    default: tab ('\t'}).                               
%       - append     '-a' for appending the content to the               
%                    output file (optional).                             
% Check input arguments
if nargin < 2
    disp('Error - Give at least two input arguments!');
    Commands=('Error- Give at leat three input arguments');
    T = evalc('Commands');
    Tcell = regexp(T, '\n', 'split');
    set(handles.text38,'String',Tcell);
    
    return;
elseif nargin > 4
    disp('Error - Do not give more than 4 input arguments!');
    
    Commands=('Error -  Dont give more than four input arguments');
    T = evalc('Commands');
    Tcell = regexp(T, '\n', 'split');
    set(handles.text38,'String',Tcell);
    
    return;
end
if ~ischar(file)
    disp(['Error - File input has to be a string (e.g. ' ...
        char(39) 'output.txt' char(39) '!']);
    
    Commands=('Error - File input has to be a string');
    T = evalc('Commands');
    Tcell = regexp(T, '\n', 'split');
    set(handles.text38,'String',Tcell);
    
    return;
end;
if ~iscell(cell_array)
    disp('Error - Input cell_array not of the type "cell"!');
    return;
end;
delimiter = '\t';
append = 'w';
if nargin > 2
    for i = 1:size(varargin,2)
        if strcmp('-a',varargin{1,i}) == 1
            append = 'a';
        else
            delimiter = varargin{1,i};
        end;
    end;
end

%Open output file and prepare output array.
output_file = fopen(file,append);
output = cell(size(cell_array,1),size(cell_array,2));

%Evaluate and write input array.
for i = 1:size(cell_array,1)
    for j = 1:size(cell_array,2)
        if numel(cell_array{i,j}) == 0
            output{i,j} = '';
            % Check whether the content of cell i,j is
            % numeric and convert numbers to strings.
        elseif isnumeric(cell_array{i,j}) || islogical(cell_array{i,j})
            output{i,j} = num2str(cell_array{i,j}(1,1));
            
            % Check whether the content of cell i,j is another cell (e.g. a
            % string of length > 1 that was stored as cell. If cell sizes
            % equal [1,1], convert numbers and char-cells to strings.
            % Note that any other cells-within-the-cell will produce errors
            % or wrong results.
        elseif iscell(cell_array{i,j})
            if size(cell_array{i,j},1) == 1 && size(cell_array{i,j},1) == 1
                if isnumeric(cell_array{i,j}{1,1})
                    output{i,j} = num2str(cell_array{i,j}{1,1}(1,1));
                elseif ischar(cell_array{i,j}{1,1})
                    output{i,j} = cell_array{i,j}{1,1};
                end;
            end;
            
            % If the cell already contains a string, nothing has to be done.
        elseif ischar(cell_array{i,j})
            output{i,j} = cell_array{i,j};
        end;        
        % Cell i,j is written to the output file. A delimiter is appended for
        % all but the last element of each row. At the end of a row, a newline
        % is written to the output file.
        if j < size(cell_array,2)
            fprintf(output_file,['%s',delimiter],output{i,j});
        else
            fprintf(output_file,'%s\r\n',output{i,j});
        end
    end;
end;

% Close output file.
fclose(output_file);


%% Function 6: UTM 2 latlong

function  [Lat,Lon] = utm2deg(xx,yy,utmzone)
%
% Description: Function to convert vectors of UTM coordinates into Lat/Lon vectors (WGS84).
% Some code has been extracted from UTMIP.m function by Gabriel Ruiz Martinez.
% Author: 
%   Rafael Palacios
%   Universidad Pontificia Comillas
%   Madrid, Spain
% Version: Apr/06, Jun/06, Aug/06
% Aug/06: corrected m-Lint warnings

error(nargchk(3, 3, nargin)); %3 arguments required
n1=length(xx);
n2=length(yy);
n3=size(utmzone,1);
if (n1~=n2 || n1~=n3)
   error('x,y and utmzone vectors should have the same number or rows');
end
c=size(utmzone,2);
if (c~=4)
   error('utmzone should be a vector of strings like "30 T"');
end
%Memory pre-allocation
Lat=zeros(n1,1);
Lon=zeros(n1,1);
% Main Loop
for i=1:n1
   if (utmzone(i,4)>'X' || utmzone(i,4)<'C')
      fprintf('utm2deg: Warning utmzone should be a vector of strings like "30 T", not "30 t"\n');
   end
   if (utmzone(i,4)>'M')
      hemis='N';   % Northern hemisphere
   else
      hemis='S';
   end

   x=xx(i);
   y=yy(i);
   zone=str2double(utmzone(i,1:2));

   sa = 6378137.000000 ; sb = 6356752.314245;  
   e2 = ( ( ( sa ^ 2 ) - ( sb ^ 2 ) ) ^ 0.5 ) / sb;
   e2cuadrada = e2 ^ 2;
   c = ( sa ^ 2 ) / sb;
   X = x - 500000;
   
   if hemis == 'S' || hemis == 's'
       Y = y - 10000000;
   else
       Y = y;
   end    
   S = ( ( zone * 6 ) - 183 ); 
   lat =  Y / ( 6366197.724 * 0.9996 );                                    
   v = ( c / ( ( 1 + ( e2cuadrada * ( cos(lat) ) ^ 2 ) ) ) ^ 0.5 ) * 0.9996;
   a = X / v;
   a1 = sin( 2 * lat );
   a2 = a1 * ( cos(lat) ) ^ 2;
   j2 = lat + ( a1 / 2 );
   j4 = ( ( 3 * j2 ) + a2 ) / 4;
   j6 = ( ( 5 * j4 ) + ( a2 * ( cos(lat) ) ^ 2) ) / 3;
   alfa = ( 3 / 4 ) * e2cuadrada;
   beta = ( 5 / 3 ) * alfa ^ 2;
   gama = ( 35 / 27 ) * alfa ^ 3;
   Bm = 0.9996 * c * ( lat - alfa * j2 + beta * j4 - gama * j6 );
   b = ( Y - Bm ) / v;
   Epsi = ( ( e2cuadrada * a^ 2 ) / 2 ) * ( cos(lat) )^ 2;
   Eps = a * ( 1 - ( Epsi / 3 ) );
   nab = ( b * ( 1 - Epsi ) ) + lat;
   senoheps = ( exp(Eps) - exp(-Eps) ) / 2;
   Delt = atan(senoheps / (cos(nab) ) );
   TaO = atan(cos(Delt) * tan(nab));
   longitude = (Delt *(180 / pi ) ) + S;
   latitude = ( lat + ( 1 + e2cuadrada* (cos(lat)^ 2) - ( 3 / 2 ) * e2cuadrada * sin(lat) * cos(lat) * ( TaO - lat ) ) * ( TaO - lat ) ) * ...
                    (180 / pi);
   
   Lat(i)=latitude;
   Lon(i)=longitude;
   
end

%% Function 7: Vertical line
function hhh=vline(x,in1,in2)
% Draws a vertical line on the current axes at the location specified by 'x'.  
% By Brandon Kuczenski for Kensington Labs.
% brandon_kuczenski@kensingtonlabs.com
% November 2001

if length(x)>1  % vector input
    for I=1:length(x)
        switch nargin
        case 1
            linetype='r:';
            label='';
        case 2
            if ~iscell(in1)
                in1={in1};
            end
            if I>length(in1)
                linetype=in1{end};
            else
                linetype=in1{I};
            end
            label='';
        case 3
            if ~iscell(in1)
                in1={in1};
            end
            if ~iscell(in2)
                in2={in2};
            end
            if I>length(in1)
                linetype=in1{end};
            else
                linetype=in1{I};
            end
            if I>length(in2)
                label=in2{end};
            else
                label=in2{I};
            end
        end
        h(I)=vline(x(I),linetype,label);
    end
else
    switch nargin
    case 1
        linetype='r:';
        label='';
    case 2
        linetype=in1;
        label='';
    case 3
        linetype=in1;
        label=in2;
    end
    g=ishold(gca);
    hold on
    y=get(gca,'ylim');
    h=plot([x x],y,linetype);
    if length(label)
        xx=get(gca,'xlim');
        xrange=xx(2)-xx(1);
        xunit=(x-xx(1))/xrange;
        if xunit<0.8
            text(x+0.01*xrange,y(1)+0.1*(y(2)-y(1)),label,'color',get(h,'color'))
        else
            text(x-.05*xrange,y(1)+0.1*(y(2)-y(1)),label,'color',get(h,'color'))
        end
    end     
    if g==0
    hold off
    end
    set(h,'tag','vline','handlevisibility','off')
end % else
if nargout
    hhh=h;
end

%% Function 8: plot google maps
function varargout = plot_google_map(varargin)
% function h = plot_google_map(varargin)
% Plots a google map on the current axes using the Google Static Maps API
% USAGE:
% h = plot_google_map(Property, Value,...)
% Plots the map on the given axes. Used also if no output is specified
% References:
% http://www.mathworks.com/matlabcentral/fileexchange/24113
% http://www.maptiler.org/google-maps-coordinates-tile-bounds-projection/
% http://developers.google.com/maps/documentation/staticmaps/
%  
%Author:
%  Zohar Bar-Yehuda
%  Version 1.3 - 06/10/2013

% store parameters in global variable (used for auto-refresh)
global inputParams
persistent apiKey
if isnumeric(apiKey)
    % first run, check if API key file exists
    if exist('api_key.mat','file')
        load api_key
    else
        apiKey = '';
    end
end
axHandle = gca;
inputParams.(['ax' num2str(axHandle*1e6,'%.0f')]) = varargin;
% Handle input arguments
height = 640;
width = 640;
scale = 2;
maptype = 'roadmap';
alphaData = 1;
autoRferesh = 1;
autoAxis = 1;
ShowLabels = 1;
hold on

markeridx = 1;
markerlist = {};
if nargin >= 2
    for idx = 1:2:length(varargin)
        switch lower(varargin{idx})
            case 'height'
                height = varargin{idx+1};
            case 'width'
                width = varargin{idx+1};
            case 'maptype'
                maptype = varargin{idx+1};
            case 'alpha'
                alphaData = varargin{idx+1};
            case 'refresh'
                autoRferesh = varargin{idx+1};
            case 'showlabels'
                ShowLabels = varargin{idx+1};
            case 'marker'
                markerlist{markeridx} = varargin{idx+1};
                markeridx = markeridx + 1;
            case 'autoaxis'
                autoAxis = varargin{idx+1};
            case 'style'
                styleParams = varargin{idx+1};                
            case 'apikey'
                apiKey = varargin{idx+1}; % set new key
                % save key to file
                funcFile = which('plot_google_map.m');
                pth = fileparts(funcFile);
                keyFile = fullfile(pth,'api_key.mat');
                save(keyFile,'apiKey')
            otherwise
                error(['Unrecognized variable: ' varargin{idx}])
        end
    end
end
if height > 640
    height = 640;
end
if width > 640
    width = 640;
end



curAxis = axis;
% Enforce Latitude constraints of EPSG:900913 
if curAxis(3) < -85
    curAxis(3) = -85;
end
if curAxis(4) > 85
    curAxis(4) = 85;
end
% Enforce longitude constrains
if curAxis(1) < -180
    curAxis(1) = -180;
end
if curAxis(1) > 180
    curAxis(1) = 0;
end
if curAxis(2) > 180
    curAxis(2) = 180;
end
if curAxis(2) < -180
    curAxis(2) = 0;
end

if isequal(curAxis,[0 1 0 1]) % probably an empty figure
    % display world map
    curAxis = [-200 200 -85 85];
    axis(curAxis)
end
if autoAxis
    % adjust current axis limit to avoid strectched maps
    [xExtent,yExtent] = latLonToMeters(curAxis(3:4), curAxis(1:2) );
    xExtent = diff(xExtent); % just the size of the span
    yExtent = diff(yExtent); 
    % get axes aspect ratio
    drawnow
    org_units = get(axHandle,'Units');
    set(axHandle,'Units','Pixels')
    ax_position = get(axHandle,'position');        
    set(axHandle,'Units',org_units)
    aspect_ratio = ax_position(4) / ax_position(3);
    
    if xExtent*aspect_ratio > yExtent        
        centerX = mean(curAxis(1:2));
        centerY = mean(curAxis(3:4));
        spanX = (curAxis(2)-curAxis(1))/2;
        spanY = (curAxis(4)-curAxis(3))/2;
       
        % enlarge the Y extent
        spanY = spanY*xExtent*aspect_ratio/yExtent; % new span
        if spanY > 85
            spanX = spanX * 85 / spanY;
            spanY = spanY * 85 / spanY;
        end
        curAxis(1) = centerX-spanX;
        curAxis(2) = centerX+spanX;
        curAxis(3) = centerY-spanY;
        curAxis(4) = centerY+spanY;
    elseif yExtent > xExtent*aspect_ratio
        
        centerX = mean(curAxis(1:2));
        centerY = mean(curAxis(3:4));
        spanX = (curAxis(2)-curAxis(1))/2;
        spanY = (curAxis(4)-curAxis(3))/2;
        % enlarge the X extent
        spanX = spanX*yExtent/(xExtent*aspect_ratio); % new span
        if spanX > 180
            spanY = spanY * 180 / spanX;
            spanX = spanX * 180 / spanX;
        end
        
        curAxis(1) = centerX-spanX;
        curAxis(2) = centerX+spanX;
        curAxis(3) = centerY-spanY;
        curAxis(4) = centerY+spanY;
    end            
    % Enforce Latitude constraints of EPSG:900913
    if curAxis(3) < -85
        curAxis(3:4) = curAxis(3:4) + (-85 - curAxis(3));
    end
    if curAxis(4) > 85
        curAxis(3:4) = curAxis(3:4) + (85 - curAxis(4));
    end
    axis(curAxis) % update axis as quickly as possible, before downloading new image
    drawnow
end
% Delete previous map from plot (if exists)
if nargout <= 1 % only if in plotting mode
    curChildren = get(axHandle,'children');
    map_objs = findobj(curChildren,'tag','gmap');
    bd_callback = [];
    for idx = 1:length(map_objs)
        if ~isempty(get(map_objs(idx),'ButtonDownFcn'))
            % copy callback properties from current map
            bd_callback = get(map_objs(idx),'ButtonDownFcn');
        end
    end
    delete(map_objs)
    
end
% Calculate zoom level for current axis limits
[xExtent,yExtent] = latLonToMeters(curAxis(3:4), curAxis(1:2) );
minResX = diff(xExtent) / width;
minResY = diff(yExtent) / height;
minRes = max([minResX minResY]);
tileSize = 256;
initialResolution = 2 * pi * 6378137 / tileSize; % 156543.03392804062 for tileSize 256 pixels
zoomlevel = floor(log2(initialResolution/minRes));
% Enforce valid zoom levels
if zoomlevel < 0 
    zoomlevel = 0;
end
if zoomlevel > 19 
    zoomlevel = 19;
end
% Calculate center coordinate in WGS1984
lat = (curAxis(3)+curAxis(4))/2;
lon = (curAxis(1)+curAxis(2))/2;
% CONSTRUCT QUERY URL
preamble = 'http://maps.googleapis.com/maps/api/staticmap';
location = ['?center=' num2str(lat,10) ',' num2str(lon,10)];
zoomStr = ['&zoom=' num2str(zoomlevel)];
sizeStr = ['&scale=' num2str(scale) '&size=' num2str(width) 'x' num2str(height)];
maptypeStr = ['&maptype=' maptype ];
styleStr = ['&style=' styleParams];

if ~isempty(apiKey)
    keyStr = ['&key=' apiKey];
else
    keyStr = '';
end
markers = '&markers=';
for idx = 1:length(markerlist)
    if idx < length(markerlist)
        markers = [markers markerlist{idx} '%7C'];
    else
        markers = [markers markerlist{idx}];
    end
end
if ShowLabels == 0
    labelsStr = '&style=feature:all|element:labels|visibility:off';
else
    labelsStr = '';
end
if ismember(maptype,{'satellite','hybrid'})
    filename = 'tmp.jpg';
    format = '&format=jpg';
    convertNeeded = 0;
else
    filename = 'tmp.png';
    format = '&format=png';
    convertNeeded = 1;
end
sensor = '&sensor=false';
url = [preamble location zoomStr sizeStr maptypeStr format markers labelsStr styleStr sensor keyStr];

%disp(url)
% Get the image
try
    urlwrite(url,filename);
catch % error downloading map
    warning(sprintf(['Unable to download map form Google Servers.\n' ...
        'Possible reasons: no network connection, or quota exceeded.\n' ...
        'Consider using an API key if quota problems persist.']));
    varargout{1} = [];
    varargout{2} = [];
    varargout{3} = [];
    return
end
[M Mcolor] = imread(filename);
M = cast(M,'double');
delete(filename); % delete temp file
width = size(M,2);
height = size(M,1);
% Calculate a meshgrid of pixel coordinates in EPSG:900913
centerPixelY = round(height/2);
centerPixelX = round(width/2);
[centerX,centerY] = latLonToMeters(lat, lon ); % center coordinates in EPSG:900913
curResolution = initialResolution / 2^zoomlevel/scale; % meters/pixel (EPSG:900913)
xVec = centerX + ((1:width)-centerPixelX) * curResolution; % x vector
yVec = centerY + ((height:-1:1)-centerPixelY) * curResolution; % y vector
[xMesh,yMesh] = meshgrid(xVec,yVec); % construct meshgrid 
% convert meshgrid to WGS1984
[lonMesh,latMesh] = metersToLatLon(xMesh,yMesh);
% Convert image from colormap type to RGB truecolor (if PNG is used)
if convertNeeded
    imag = zeros(height,width,3);
    for idx = 1:3
        imag(:,:,idx) = reshape(Mcolor(M(:)+1+(idx-1)*size(Mcolor,1)),height,width);
    end
else
    imag = M/255;
end

% Next, project the data into a uniform WGS1984 grid
sizeFactor = 1; % factoring of new image
uniHeight = round(height*sizeFactor);
uniWidth = round(width*sizeFactor);
latVect = linspace(latMesh(1,1),latMesh(end,1),uniHeight);
lonVect = linspace(lonMesh(1,1),lonMesh(1,end),uniWidth);
[uniLonMesh,uniLatMesh] = meshgrid(lonVect,latVect);
uniImag = zeros(uniHeight,uniWidth,3);

uniImag =  myTurboInterp2(lonMesh,latMesh,imag,uniLonMesh,uniLatMesh);

if nargout <= 1 % plot map
    % display image
    h = image(lonVect,latVect,uniImag);    
    set(gca,'YDir','Normal')
    set(h,'tag','gmap')
    set(h,'AlphaData',alphaData)
    
    % add a dummy image to allow pan/zoom out to x2 of the image extent
    h_tmp = image(lonVect([1 end]),latVect([1 end]),zeros(2),'Visible','off');
    set(h_tmp,'tag','gmap')
    
    % older version (display without conversion to uniform grid)
    % h =pcolor(lonMesh,latMesh,(M));
    % colormap(Mcolor)
    % caxis([0 255])
    % warning off % to avoid strange rendering warnings
    % shading flat
   
    uistack(h,'bottom') % move map to bottom (so it doesn't hide previously drawn annotations)
    axis(curAxis) % restore original zoom
    if nargout == 1
        varargout{1} = h;
    end
    % if auto-refresh mode - override zoom callback to allow autumatic 
    % refresh of map upon zoom actions.
    zoomHandle = zoom;   
    panHandle = pan;    
    if autoRferesh        
        set(zoomHandle,'ActionPostCallback',@update_google_map);          
        set(panHandle, 'ActionPostCallback', @update_google_map);        
    else % disable zoom override
        set(zoomHandle,'ActionPostCallback',[]);
        set(panHandle, 'ActionPostCallback',[]);
    end    
    % set callback for figure resize function, to update extents if figure
    % is streched.
    figHandle = get(axHandle,'Parent');
    set(figHandle, 'ResizeFcn', @update_google_map_fig);        
    % set callback properties 
    set(h,'ButtonDownFcn',bd_callback);
else % don't plot, only return map
    varargout{1} = lonVect;
    varargout{2} = latVect;
    varargout{3} = uniImag;
end
% Coordinate transformation functions
function [lon,lat] = metersToLatLon(x,y)
% Converts XY point from Spherical Mercator EPSG:900913 to lat/lon in WGS84 Datum
originShift = 2 * pi * 6378137 / 2.0; % 20037508.342789244
lon = (x ./ originShift) * 180;
lat = (y ./ originShift) * 180;
lat = 180 / pi * (2 * atan( exp( lat * pi / 180)) - pi / 2);

function [x,y] = latLonToMeters(lat, lon )
% Converts given lat/lon in WGS84 Datum to XY in Spherical Mercator EPSG:900913"
originShift = 2 * pi * 6378137 / 2.0; % 20037508.342789244
x = lon * originShift / 180;
y = log(tan((90 + lat) * pi / 360 )) / (pi / 180);
y = y * originShift / 180;


function ZI = myTurboInterp2(X,Y,Z,XI,YI)
% An extremely fast nearest neighbour 2D interpolation, assuming both input
% and output grids consist only of squares, meaning:
% - uniform X for each column
% - uniform Y for each row
XI = XI(1,:);
X = X(1,:);
YI = YI(:,1);
Y = Y(:,1);

xiPos = nan*ones(size(XI));
xLen = length(X);
yiPos = nan*ones(size(YI));
yLen = length(Y);
% find x conversion
xPos = 1;
for idx = 1:length(xiPos)
    if XI(idx) >= X(1) && XI(idx) <= X(end)
        while xPos < xLen && X(xPos+1)<XI(idx)
            xPos = xPos + 1;
        end
        diffs = abs(X(xPos:xPos+1)-XI(idx));
        if diffs(1) < diffs(2)
            xiPos(idx) = xPos;
        else
            xiPos(idx) = xPos + 1;
        end
    end
end
% find y conversion
yPos = 1;
for idx = 1:length(yiPos)
    if YI(idx) <= Y(1) && YI(idx) >= Y(end)
        while yPos < yLen && Y(yPos+1)>YI(idx)
            yPos = yPos + 1;
        end
        diffs = abs(Y(yPos:yPos+1)-YI(idx));
        if diffs(1) < diffs(2)
            yiPos(idx) = yPos;
        else
            yiPos(idx) = yPos + 1;
        end
    end
end
ZI = Z(yiPos,xiPos,:);


function update_google_map(obj,evd)
% callback function for auto-refresh
drawnow;
global inputParams
if isfield(inputParams,['ax' num2str(gca*1e6,'%.0f')])
    params = inputParams.(['ax' num2str(gca*1e6,'%.0f')]);
    plot_google_map(params{:});
end

function update_google_map_fig(obj,evd)
% callback function for auto-refresh
drawnow;
global inputParams
axes_objs = findobj(get(gcf,'children'),'type','axes');
for idx = 1:length(axes_objs)
    if ~isempty(findobj(get(axes_objs(idx),'children'),'tag','gmap'));
        if isfield(inputParams,['ax' num2str(axes_objs(idx)*1e6,'%.0f')])
            params = inputParams.(['ax' num2str(axes_objs(idx)*1e6,'%.0f')]);
        else
            params = {};
        end
        axes(axes_objs(idx));
        plot_google_map(params{:});
        break;
    end
end

%% Function 9: Moving
function [y]=moving(x,m,fun)
%MOVING will compute moving averages of order n (best taken as odd)
%Usage: y=moving(x,n[,fun])
%where x 	is the input vector (or matrix) to be smoothed. 
%      m 	is number of points to average over (best odd, but even works)
%      y 	is output vector of same length as x
%      fun  (optional) is a custom function rather than moving averages
if m==1
    y=x;
    return
end
if size(x,1)==1
    x=x';
end
if nargin<3
    fun=[];
elseif ischar(fun)
    fun=eval(['@(x)' fun '(x)']);
end
if isempty(fun)
    f=zeros(m,1)+1/m;
    n=size(x,1);
    isodd=bitand(m,1);
    m2=floor(m/2);
    if (size(x,2)==1)
        y=filter(f,1,x);
        y=y([zeros(1,m2-1+isodd)+m,m:n,zeros(1,m2)+n]);
    else
        y=filter2(f,x);
        y(1:(m2-~isodd),:)=y(m2+isodd+zeros(m2-~isodd,1),:);
        y((n-m2+1):end,:)=y(n-m2+zeros(m2,1),:);
    end
else
    y=zeros(size(x));
    sx=size(x,2);
    x=[nan(floor(m*.5),sx);x;nan(floor(m*.5),sx)];
    m1=m-1;
    for ii=1:size(y,1);
        y(ii,:)=fun(x(ii+(0:m1),:));
    end    
end
return

%% Function 10: point2utm
function point2utm=Terracem_point2utm(p1,pto, thetaB, Dprof, type_profile)
%Function to convert profile points to utm coordinates based on the
%geometry of the rectangular profile (swath area)
%pto=red point for orienttion of profile
%p1=calculated shoreline angle in profile
%pto, thetaB, and Dprof comes from extract swats
%type_profile comes from each analysis functions (sfitprof,
%terracediff and stacks)
%Jara & Melnick 2014

dprof=Dprof;
scarp=p1;
%CASE(1) for non-flipped profiles
if type_profile==1 
x1pto=pto(1,1);
y1pto=pto(1,2);
if (thetaB<90) && (thetaB>0) %|| (thetaB>173.7)&& (thetaB<180)   
x2pto=pto(1,1)+((dprof)*cosd(thetaB));
y2pto=pto(1,2)+((dprof)*sind(thetaB));
% rectangle axis pairs for analysis
px1(1,1)=x1pto;
px1(1,2)=x2pto;
py1(1,1)=y1pto;
py1(1,2)=y2pto;     
Rx= (scarp(1,1)/dprof) * (px1(1,2)-px1(1,1)) + px1(1,1); %x for 0 - 90?
Ry= (scarp(1,1)/dprof) * (py1(1,2)-py1(1,1)) + py1(1,1); %y    
elseif (thetaB>173.7) && (thetaB<180)
x2pto=-((dprof)*cosd(thetaB))+pto(1,1);
y2pto=-((dprof)*sind(thetaB))+pto(1,2);
% rectangle axis pairs for analysis
px1(1,1)=x1pto;
px1(1,2)=x2pto;
py1(1,1)=y1pto;
py1(1,2)=y2pto;  
scarp2=scarp(1,1);
Rx= (scarp2(1,1)/dprof) * (px1(1,2)-px1(1,1)) + px1(1,1); %x for 180 - 173.7?
Ry= (scarp2(1,1)/dprof) * (py1(1,2)-py1(1,1)) + py1(1,1); %y       
else %(thetaB<180) && (thetaB>90) 
scarp1=dprof-scarp(1,1);% Modified scarp position for NW oriented profiles
x2pto=pto(1,1)-((dprof)*cosd(thetaB));
y2pto=pto(1,2)-((dprof)*sind(thetaB)); 
% rectangle axis pairs for analysis
px1(1,1)=x1pto; 
px1(1,2)=x2pto;
py1(1,1)=y1pto; py1(1,2)=y2pto;      
Rx= (scarp1(1,1)/dprof) * (px1(1,1)-px1(1,2)) + px1(1,1); %x for 90 - 180?
Ry= (scarp1(1,1)/dprof) * (py1(1,1)-py1(1,2)) + py1(1,1); %y     
end
%CASE(2) for flipped profiles
else
x2pto=pto(1,1);
y2pto=pto(1,2);
scarp_pre=scarp(1,1);
if (thetaB<90) && (thetaB>0)         
x1pto=pto(1,1)+((dprof)*cosd(thetaB));
y1pto=pto(1,2)+((dprof)*sind(thetaB));
% rectangle axis pairs for analysis
px1(1,1)=x1pto;px1(1,2)=x2pto;
py1(1,1)=y1pto;py1(1,2)=y2pto;    
Rx= (scarp_pre(1,1)/dprof) * (px1(1,2)-px1(1,1)) + px1(1,1); %x for 0 - 90?
Ry= (scarp_pre(1,1)/dprof) * (py1(1,2)-py1(1,1)) + py1(1,1); %y        
else %(thetaB<180) && (thetaB>90)     
x1pto=pto(1,1);
y1pto=pto(1,2);         
x2pto=pto(1,1)-((dprof)*cosd(thetaB));
y2pto=pto(1,2)-((dprof)*sind(thetaB)); 
% rectangle axis pairs for analysis
px1(1,1)=x1pto; 
px1(1,2)=x2pto;
py1(1,1)=y1pto; 
py1(1,2)=y2pto;  
scarp1=dprof-scarp(1,1);
Rx= (scarp1(1,1)/dprof) * (px1(1,1)-px1(1,2)) + px1(1,1); %x for 90 - 180?
Ry= (scarp1(1,1)/dprof) * (py1(1,1)-py1(1,2)) + py1(1,1); %y     
end    
end 
%outputs in utm coordinates
Rz=p1(1,2); %z value for input point 
point2utm=[Rx Ry Rz];    

%% Function 11: errorbar_tick
function errorbar_tick(h,w,xtype)
%ERRORBAR_TICK Adjust the width of errorbars
%   ERRORBAR_TICK(H) adjust the width of error bars with handle H.
%      Error bars width is given as a ratio of X axis length (1/80).
%   ERRORBAR_TICK(H,W) adjust the width of error bars with handle H.
%      The input W is given as a ratio of X axis length (1/W). The result 
%      is independent of the x-axis units. A ratio between 20 and 80 is usually fine.
%   ERRORBAR_TICK(H,W,'UNITS') adjust the width of error bars with handle H.
%      The input W is given in the units of the current x-axis.
% Author: Arnaud Laurent
% Creation : Jan 29th 2009
% MATLAB version: R2007a
% Author : Jerome Briot (Dut) 
% Check numbers of arguments
error(nargchk(1,3,nargin))
% Check for the use of V6 flag ( even if it is depreciated ;) )
flagtype = get(h,'type');

% Check number of arguments and provide missing values
if nargin==1
	w = 80;
end
if nargin<3
   xtype = 'ratio';
end

% Calculate width of error bars
if ~strcmpi(xtype,'units')
    dx = diff(get(gca,'XLim'));	% Retrieve x limits from current axis
    w = dx/w;                   % Errorbar width
end

% Plot error bars
if strcmpi(flagtype,'hggroup') % ERRORBAR(...)
    
    hh=get(h,'children');		% Retrieve info from errorbar plot
    x = get(hh(2),'xdata');		% Get xdata from errorbar plot
    
    x(4:9:end) = x(1:9:end)-w/2;	% Change xdata with respect to ratio
    x(7:9:end) = x(1:9:end)-w/2;
    x(5:9:end) = x(1:9:end)+w/2;
    x(8:9:end) = x(1:9:end)+w/2;

    set(hh(2),'xdata',x(:))	% Change error bars on the figure

else  % ERRORBAR('V6',...)    
    x = get(h(1),'xdata');		% Get xdata from errorbar plot    
    x(4:9:end) = x(1:9:end)-w/2;	% Change xdata with respect to the chosen ratio
    x(7:9:end) = x(1:9:end)-w/2;
    x(5:9:end) = x(1:9:end)+w/2;
    x(8:9:end) = x(1:9:end)+w/2;

    set(h(1),'xdata',x(:))	% Change error bars on the figure
    
end

%% Function 12: Peak detection
function [maxtab, mintab]=peakdet(v, delta, x)
%PEAKDET Detect peaks in a vector
%        [MAXTAB, MINTAB] = PEAKDET(V, DELTA) finds the local
%        maxima and minima ("peaks") in the vector V.
%        MAXTAB and MINTAB consists of two columns. Column 1
%        contains indices in V, and column 2 the found values.
%        With [MAXTAB, MINTAB] = PEAKDET(V, DELTA, X) the indices
%        in MAXTAB and MINTAB are replaced with the corresponding
%        X-values.
%        A point is considered a maximum peak if it has the maximal
%        value, and was preceded (to the left) by a value lower by
%        DELTA.

% Eli Billauer, 3.4.05 (Explicitly not copyrighted).
% This function is released to the public domain; Any use is allowed.
maxtab = [];
mintab = [];
v = v(:); % Just in case this wasn't a proper vector
if nargin < 3
  x = (1:length(v))';
else 
  x = x(:);
  if length(v)~= length(x)
    error('Input vectors v and x must have same length');
  end
end
  
if (length(delta(:)))>1
  error('Input argument DELTA must be a scalar');
end

if delta <= 0
  error('Input argument DELTA must be positive');
end
mn = Inf; mx = -Inf;
mnpos = NaN; mxpos = NaN;
lookformax = 1;

for i=1:length(v)
  this = v(i);
  if this > mx, mx = this; mxpos = x(i); end
  if this < mn, mn = this; mnpos = x(i); end
  
  if lookformax
    if this < mx-delta
      maxtab = [maxtab ; mxpos mx];
      mn = this; mnpos = x(i);
      lookformax = 0;
    end  
  else
    if this > mn+delta
      mintab = [mintab ; mnpos mn];
      mx = this; mxpos = x(i);
      lookformax = 1;
    end
  end
end

%% Function 13: Finite scarp diffussion core
function H = finitescarp(thetadegrees,a, x, kt ,bdegrees)
% Solves the finitescarp initial conditions problem for
% thetadegrees is the initial scarp slope in degrees
% a = half scarp initial height
% x = distance along profile
% kt = morphologic age (L^2)
% bdegrees is far field slope in degrees
% from Hanks, T. C. and D. J. Andrews, Effect of far-field slope on morphologic dating 
% of scarplike landforms, Journal of Geophysical Research, 94, 565-573, 1989.
format long
theta=rads(thetadegrees);
b1=rads(bdegrees);
fkt=4.*kt;
sfkt=sqrt(fkt);
ang=theta-b1;
A=a./ang;
X=x+A;
Xm=x-A;
%Difussion core, coded by Arrowsmith, adapted by J. Jara-Munoz for use with swath profiles
    if theta>=rads(75)
   H=(a.*erf(x./(2.*sqrt(kt)))+b1.*x).*-1;%Note: -1 for flip profile
    else
   H=(ang.*sqrt((kt)/pi).*(exp(-((X.^2)./fkt))-exp(-((Xm.^2)./fkt)))+(ang/2).*(X.*erf(X./sfkt)-Xm .*erf(Xm./sfkt))+b1.*x).*-1;
    end
%Small nested function, radian transformation
function outrads = rads(indeg)
format long
outrads = indeg.*pi/180;

%% Function 14: dem
function [h,I,z]=dem(x,y,z,varargin)
%DEM Shaded relief image plot
%
%	DEM(X,Y,Z) plots the Digital Elevation Model defined by X and Y 
%	coordinate vectors and elevation matrix Z, as a lighted image using
%	specific "landcolor" and "seacolor" colormaps. DEM uses IMAGESC 
%	function which is much faster than SURFL when dealing with large 
%	high-resolution DEM. It produces also high-quality and moderate-size 
%	Postscript image adapted for publication.
%
%	[H,I] = DEM(...); returns graphic handle H and illuminated image as I, 
%	an MxNx3 matrix (if Z is MxN and DECIM is 1).
%
%	DEM(X,Y,Z,'Param1',Value1,'Param2',Value2,...) specifies options or
%	parameter/value couple (case insensitive):
%
%
%	--- Lighting options ---
%
%	'Azimuth',A
%		Light azimuth in degrees clockwise relative to North. Default is
%		A = -45 for	a natural northwestern illumination.
%
%	'Contrast',C
%		Light contrast, as the exponent of the gradient value:
%			C = 1 for linear contrast (default),
%			C = 0 to remove lighting,
%			C = 0.5 for moderate lighting,
%			C = 2 or more for strong contrast.
%
%	'LCut',LC
%		Lighting scale saturation cut with a median-style filter in % of 
%	    elements, such as LC% of maximum gradient values are ignored:
%			LC = 0.2 is default, 
%			LC = 0 for full scale gradient.
%
%	'km'
%		Stands that X and Y coordinates are in km instead of m (default).
%		This allows correct lighting. Ignored if LATLON option is used.
%
%
%	--- Elevation colorscale options ---
%
%	'ZLim',[ZMIN,ZMAX]
%		Fixes min and max elevation values for colormap. Use NaN to keep 
%		real min and/or max data values.
%
%	'ZCut',ZC
%		Median-style filter to cut extremes values of Z (in % of elements),
%		such that ZC% of most min/max elevation values are ignored in the
%		colormap application:
%			ZC = 0.5 is default, 
%			ZC = 0 for full scale.
%
%
%	--- "No Value" elevation options ---
%
%	'NoValue',NOVALUE
%		Defines the values that will be replaced by NaN. Note that values 
%		equal to minimum of Z class are automatically detected as NaN 
%		(e.g., -32768 for int16 class).
%
%	'NaNColor',[R,G,B]
%		Sets the RGB color for NaN/NoValue pixels (default is a dark gray).
%		Note that your must specify a valid 3-scalar vector (between 0 and
%		1);	color characters like 'w' or 'k' are not allowed, use [1,1,1]
%		or [0,0,0] instead.
%
%	'Interp'
%		Interpolates linearly all NaN values (fills the gaps using linear 
%		triangulation), using an optimized algorithm.
%
%
%	--- Colormap options ---
%
%	'LandColor',LMAP
%		Uses LMAP colormap instead of default (landcolor, if exists or 
%		jet) for Z > 0 elevations.
%
%	'SeaColor',SMAP
%		Sets the colormap used for Z <= 0 elevations. Default is seacolor 
%		(if exists) or single color [0.7,0.9,1] (a light cyan) to simulate
%		sea color.
%
%	'ColorMap',CMAP
%		Uses CMAP colormap for full range of elevations, instead of default 
%		land/sea. This option overwrites LANDCOLOR/SEACOLOR options.
%
%	'Lake'
%		Detects automaticaly flat areas different from sea level (non-zero 
%		elevations) and colors them as lake surfaces.
%
%	'Watermark',N
%		Makes the whole image lighter by a factor of N.
%
%
%	--- Basemap and scale options ---
%
%	'Legend'
%		Adds legends to the right of graph: elevation scale (colorbar)
%		and a distance scale (in km).
%
%	'Cartesian'
%		Plots staircase basemap-style axis, considering coordinates X and Y 
%		as cartesian in meters. Use parameter "km' for X/Y in km.
%
%	'LatLon'
%		Plots geographic basemap-style axis in deg/min/sec, considering 
%		coordinates X as longitude and Y as latitude. Axis aspect ratio 
%		will be adjusted to approximatively preserve distances (this is  
%		not a real projection!). This overwrites ZRatio option.
%
%	'AxisEqual', 'auto' (default) | 'manual' | 'off'
%		When 'Cartesian' or 'LatLon' option is used, automatic axes scaling
%		is applied to respect data aspect ratio. Default mode is 'auto' and
%		uses AXIS EQUAL and DASPECT functions. The 'manual' mode modifies
%		axes width or height with respect to the paper size in order to
%		produce correct data scaling at print (but not necessarily at 
%		screen). The 'off' mode disables any scaling.
%
%	Additionnal options for basemap CARTESIAN or LATLON:
%
%	'BorderWidth',BW
%		Border width of the basemap axis, in % of axis height. Default is
%		BW = 1%.
%
%	'XTick',DX
%	'YTick',DY
%		X and Y Tick length (same unit as X and Y). Default is automatic.
%		Tick labels are every 2 ticks.
%
%	'FontSize',FS
%		Font size for X and Y tick labels. Default is FS = 10.
%
%	'FontBold'
%		Font weight bold for tick labels.
%
%
%	--- Decimation options ---
%
%	For optimization purpose, DEM will automatically decimate data to limit
%	to a total of 1500x1500 pixels images. To avoid this, use following
%	options, but be aware that large grids may require huge computer 
%	ressources or induce disk swap or memory errors.
%
%	'Decim',N
%		Decimates matrix Z at 1/N times of the original sampling.
%
%	'NoDecim'
%		Forces full resolution of Z, no decimation.
%	--- Informations ---
%
%	Colormaps are Mx3 RGB matrix so it is easy to modify saturation 
%	(CMAP.^N), set darker (CMAP/N), lighter (1 - 1/N + CMAP/N), inverse
%	it (flipud(CMAP)), etc...
%
%	To get free worldwide topographic data (SRTM), see READHGT function.
%
%	For backward compatibility, the former syntax is still accepted:
%	DEM(X,Y,Z,OPT,CMAP,NOVALUE,SEACOLOR) where OPT = [A,C,LC,ZMIN,ZMAX,ZC],
%	also option aliases DEC, DMS and SCALE, but there is no argument 
%	checking. Please prefer the param/value syntax.
%
%	Author: Fran?ois Beauducel <beauducel@ipgp.fr>
%	Updated: 2014-06-06
%
%	Copyright (c) 2014, Fran?ois Beauducel, covered by BSD License.
%	All rights reserved.
%
%	Redistribution and use in source and binary forms, with or without 
%	modification, are permitted provided that the following conditions are 
%	met:
%
%	   * Redistributions of source code must retain the above copyright 
%	     notice, this list of conditions and the following disclaimer.
%	   * Redistributions in binary form must reproduce the above copyright 
%	     notice, this list of conditions and the following disclaimer in 
%	     the documentation and/or other materials provided with the distribution
%	                           
%	THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
%	AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
%	IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
%	ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
%	LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
%	CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
%	SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
%	INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
%	CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
%	ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
%	POSSIBILITY OF SUCH DAMAGE.

if nargin < 3
	error('Not enough input arguments.');
end

degkm = 6378*pi/180; % one latitude degree in km
sea_color = [.7,.9,1]; % default sea color (light cyan)
grey = 0.2*[1,1,1]; % a dark gray


% -------------------------------------------------------------------------
% --- Manage input arguments

% number of arguments param/value
nargs = 0;

if ~isnumeric(x) || ~isnumeric(y) || ~isnumeric(z)
	error('X,Y and Z must be numeric.')
end

if all(size(x) ~= 1) || all(size(y) ~= 1)
	error('X and Y must be vectors, not matrix.')
end

if length(x) ~= size(z,2) || length(y) ~= size(z,1)
	error('If Z has a size of [M,N], X must have a length of N, and Y a length of M.')
end

% OPTIONS and PARAM/VALUE arguments
			
% AZIMUTH param/value
[s,az] = checkparam(varargin,'azimuth',@isscalar);
nargs = nargs + 2;
if s==0
	az = -45; % default
end

% ELEVATION param/value
[s,el] = checkparam(varargin,'elevation',@isscalar);
nargs = nargs + 2;
if s==0
	el = 0; % default
end

% CONTRAST param/value
[s,ct] = checkparam(varargin,'contrast',@isscalar);
nargs = nargs + 2;
if s
	ct = abs(ct);
else
	ct = 1; % default
end

% LCUT param/value
[s,lcut] = checkparam(varargin,'lcut',@isperc);
nargs = nargs + 2;
if s==0
	lcut = .2; % default
end

% NOVALUE param/value
[s,novalue] = checkparam(varargin,'novalue',@isscalar);
nargs = nargs + 2;
if s==0
	% default: min value for integer class / NaN for float
	S = whos('z');
	if strfind(S.class,'int')
		novalue = intmin(S.class);
	else
		novalue = NaN;
	end
end

% NANCOLOR param/value
[s,novalue_color] = checkparam(varargin,'nancolor',@isrgb);
nargs = nargs + 2;
if s==0
	novalue_color = grey; % default
end

% LANDCOLOR param/value
[s,cland] = checkparam(varargin,'landcolor',@isrgb);
nargs = nargs + 2;
if s==0
	% default: landcolor or jet
	if exist('landcolor','file')
		cland = landcolor.^1.3;
	else
		cland = jet(256);
	end
end

% SEACOLOR param/value
[s,csea] = checkparam(varargin,'seacolor',@isrgb);
nargs = nargs + 2;
if s==0
	% default: seacolor or single color
	if exist('seacolor','file')
		csea = seacolor;
	else
		csea = sea_color;
	end
end

% COLORMAP param/value
[s,cmap] = checkparam(varargin,'colormap',@isrgb);
nargs = nargs + 2;
if s
	cland = [];
	csea = [];
else
	% default
	cmap = cland;
end

% ZLIM param/value
[s,zmm] = checkparam(varargin,'zlim',@isvec);
nargs = nargs + 2;
if s
	zmin = min(zmm);
	zmax = max(zmm);
else
	zmin = NaN; % default
	zmax = NaN; % default
end

% ZCUT param/value
[s,zcut] = checkparam(varargin,'zcut',@isperc);
nargs = nargs + 2;
if s==0
	zcut = .5; % default
end

% ZRATIO param/value
[s,zratio] = checkparam(varargin,'zratio',@isscalar);
nargs = nargs + 2;
if s==0
	zratio = 1; % default
end

% WATERMARK param/value
[s,wmark] = checkparam(varargin,'watermark',@isscalar);
nargs = nargs + 2;
if s
	wmark = abs(wmark);
else
	wmark = 0; % default
end

% DECIM param/value and NODECIM option
[s,decim] = checkparam(varargin,'decim',@isscalar);
if s
	decim = round(decim);
	nargs = nargs + 2;
else
	decim = any(strcmpi(varargin,'nodecim')); % default
	nargs = nargs + 1;
end

% FONTSIZE param/value
[s,fs] = checkparam(varargin,'fontsize',@isscalar);
nargs = nargs + 2;
if s==0
	fs = 10; % default
end

% BORDERWIDTH param/value
[s,bw] = checkparam(varargin,'borderwidth',@isperc);
nargs = nargs + 2;
if s==0
	bw = 1; % default
end

% XTICK param/value
[s,ddx] = checkparam(varargin,'xtick',@isscalar);
nargs = nargs + 2;
if s==0
	ddx = 0; % default (automatic)
end

% YTICK param/value
[s,ddy] = checkparam(varargin,'ytick',@isscalar);
nargs = nargs + 2;
if s==0
	ddy = 0; % default (automatic)
end

% AXISEQUAL param/value
[s,axeq] = checkparam(varargin,'axisequal',@ischar,{'auto','manual','off'});
nargs = nargs + 2;
if s==0 || ~any(strcmpi(axeq,{'manual','off'}))
	axeq = 'auto'; % default (automatic)
end

% CROP param/value
[s,crop] = checkparam(varargin,'crop',@isvec,4);
nargs = nargs + 2;

% options without argument value
km = any(strcmpi(varargin,'km'));
dec = any(strcmpi(varargin,'cartesian') | strcmpi(varargin,'dec'));
dms = any(strcmpi(varargin,'latlon') | strcmpi(varargin,'dms'));
scale = any(strcmpi(varargin,'legend') | strcmpi(varargin,'scale'));
inter = any(strcmpi(varargin,'interp'));
lake = any(strcmpi(varargin,'lake'));
fbold = any(strcmpi(varargin,'fontbold'));


% for backward compatibility (former syntax)...
nargs = nargs + dec + dms + scale + lake + inter + km + fbold;

if (nargin - nargs) > 3 && ~isempty(varargin{1})
	opt = varargin{1};
	if ~isnumeric(opt)
		error('OPT = [A,C,S,ZMIN,ZMAX,ZCUT] argument must be numeric.');
	end
	if ~isempty(opt)
		az = opt(1);
	end
	if length(opt) > 1
		ct = opt(2);
	end
	if length(opt) > 2
		lcut = opt(3);
	end
	if length(opt) > 4
		zmin = opt(4);
		zmax = opt(5);
	end
	if length(opt) > 5
		zcut = opt(6);
	end
end

if (nargin - nargs) > 4 && ~isempty(varargin{2})
	cmap = varargin{2};
	csea = [];
end

if (nargin - nargs) > 5 && ~isempty(varargin{3})
	novalue = varargin{3};
end

if (nargin - nargs) > 6 && ~isempty(varargin{4})
	csea = varargin{4};
end


% further test of input arguments
if dms && any(abs(y) > 91)
	error('With LATLON option Y must be in valid latitudes interval (decimal degrees).')
end

if km
	zratio = 1000;
end


% -------------------------------------------------------------------------
% --- Pre-process DEM data

% crops data if needed
if numel(crop)==4
	fprintf('DEM: crops original data from [%g,%g,%g,%g] to [%g,%g,%g,%g]...\n', ...
		min(x(:)),max(x(:)),min(y(:)),max(y(:)),crop);
	kx = find(x >= crop(1) & x <= crop(2));
	ky = find(y >= crop(3) & y <= crop(4));
	x = x(kx);
	y = y(ky);
	z = z(ky,kx);
end

% decimates data to avoid disk swap/out of memory...
nmax = 1500;
if decim
	n = decim;
else
	n = ceil(sqrt(numel(z))/nmax);
end
if n > 1
	x = x(1:n:end);
	y = y(1:n:end);
	z = z(1:n:end,1:n:end);
	fprintf('DEM: on the plot data has been decimated by a factor of %d...\n',n);
end

z = double(z); % necessary for most of the following calculations...
z(z==novalue) = NaN;

if inter
	z = fillgap(x,y,z);
end

if isempty(csea)
	k = (z~=0 & ~isnan(z));
else
	k = ~isnan(z);
end

if isnan(zmin)
	zmin = nmedian(z(k),zcut/100);
end
if isnan(zmax)
	zmax = nmedian(z(k),1 - zcut/100);
end
dz = zmax - zmin;


% -------------------------------------------------------------------------
% --- Process lighting

if dz > 0
	% builds the colormap: concatenates seacolor and landcolor around 0
	if ~isempty(csea)
		l = size(csea,1);
		if zmin < 0 && zmax > 0
			r = size(cland,1)*abs(zmin)/zmax/l;
			cmap = cat(1,interp1(1:l,csea,linspace(1,l,ceil(l*r)),'*linear'),cland);
		elseif zmax <=0
			cmap = csea;
		end
	end
	
	% normalisation of Z using CMAP and convertion to RGB
	I = ind2rgb(uint16(round((z - zmin)*(size(cmap,1) - 1)/dz) + 1),cmap);
	
	if ct > 0
		% computes lighting from elevation gradient
		%[fx,fy] = gradient(z,x,y);
		if dms
			ryz = degkm*1000;
			rxz = degkm*1000*cosd(mean(y));
		else
			rxz = zratio;
			ryz = zratio;
		end
		[xx,yy] = meshgrid(x*rxz,y*ryz);
		[fx,fy,fz] = surfnorm(xx,yy,z);
		[ux,uy,uz] = sph2cart((90-az)*pi/180,el*pi/180,1);
		fxy = fx*ux + fy*uy + fz*uz;
		clear xx yy fx fy fz	% free some memory...
		
		fxy(isnan(fxy)) = 0;

		% computes maximum absolute gradient (median-style), normalizes,
		% saturates and duplicates in 3-D matrix
		li = 1 - abs(sind(el)); % light amplitude (experimental)
		r = repmat(max(min(li*fxy/nmedian(abs(fxy),1 - lcut/100),1),-1),[1,1,3]);
		rp = (1 - abs(r)).^ct;
	
		% applies contrast using exponent
		I = I.*rp;
	
		% lighter for positive gradient
		I(r>0) = I(r>0) + (1 - rp(r>0));
				
	end

	% set novalues / NaN to nancolor
	[i,j] = find(isnan(z));
	if ~isempty(i)
		I(sub2ind(size(I),repmat(i,1,3),repmat(j,1,3),repmat(1:3,size(i,1),1))) = repmat(novalue_color,size(i,1),1);
	end
	
	% lake option
	if lake
		klake = islake(z);
	else
		klake = 0;
	end
	
	% set the seacolor for 0 values
	if ~isempty(csea)
		[i,j] = find(z==0 | klake);
		if ~isempty(i)
			I(sub2ind(size(I),repmat(i,1,3),repmat(j,1,3),repmat(1:3,size(i,1),1))) = repmat(csea(end,:),size(i,1),1);
		end
	end

	if wmark
		I = watermark(I,wmark);
	end
	
	hh = imagesc(x,y,I);
	
else
	
	hh = imagesc(x,y,repmat(shiftdim(sea_color,-1),size(z)));
	text(mean(x),mean(y),'SPLASH!','Color',sea_color/4, ...
		'FontWeight','bold','HorizontalAlignment','center')
	cmap = repmat(sea_color,[256,1]);
	if nargin > 1
		I = ones(size(z));
	end
end

orient tall; axis xy
if strcmpi(axeq,'auto')
	axis equal
end
axis tight
xlim = [min(x),max(x)];
ylim = [min(y),max(y)];
zlim = [min([z(z(:) ~= novalue);zmin]),max([z(z(:) ~= novalue);zmax])];

if dms
	% approximates X-Y aspect ratio for this latitude (< 20-m precision for 1x1? grid)
	xyr = cos(mean(y)*pi/180);
else
	xyr = 1;
end

bw0 = max(diff(xlim)*xyr,diff(ylim))/100;
bwy = bw*bw0; % Y border width = 1%
bwx = bwy/xyr; % border width (in degree of longitude)


% -------------------------------------------------------------------------
% --- Axis basemap style
if dec || dms
	axis off

	if strcmpi(axeq,'manual')
		ppos = get(gcf,'PaperPosition');
		apos = get(gca,'Position');
		xyf = (xyr*diff(xlim)/apos(3)/ppos(3))/(diff(ylim)/apos(4)/ppos(4));
		if xyf >= 1
			set(gca,'Position',[apos(1),apos(2),apos(3),apos(4)/xyf]);
		else
			set(gca,'Position',[apos(1),apos(2),apos(3)*xyf,apos(4)]);
		end
	end
	if strcmpi(axeq,'auto')
		if diff(xlim)*xyr <= diff(ylim)
			set(gca,'DataAspectRatio',[1,xyr,1])
		else
			set(gca,'DataAspectRatio',[1/xyr,1,1])
		end
	end

	if bw > 0
		% transparent borders
		patch([xlim(1)-bwx,xlim(2)+bwx,xlim(2)+bwx,xlim(1)-bwx],ylim(1) - bwy*[0,0,1,1],'k','FaceColor','none','clipping','off')
		patch([xlim(1)-bwx,xlim(2)+bwx,xlim(2)+bwx,xlim(1)-bwx],ylim(2) + bwy*[0,0,1,1],'k','FaceColor','none','clipping','off')
		patch(xlim(1) - bwx*[0,0,1,1],[ylim(1)-bwy,ylim(2)+bwy,ylim(2)+bwy,ylim(1)-bwy],'k','FaceColor','none','clipping','off')
		patch(xlim(2) + bwx*[0,0,1,1],[ylim(1)-bwy,ylim(2)+bwy,ylim(2)+bwy,ylim(1)-bwy],'k','FaceColor','none','clipping','off')
	end
	dlon = {'E','W'};
	dlat = {'N','S'};
	if fbold
		fw = 'bold';
	else
		fw = 'normal';
	end
	
	if ddx == 0
		ddx = dtick(diff(xlim),dms);
	end
	if ddy == 0
		ddy = dtick(diff(ylim),dms);
	end
	xtick = (ddx*ceil(xlim(1)/ddx)):ddx:xlim(2);
	for xt = xtick(1:2:end)
		dt = ddx - max(0,xt + ddx - xlim(2));
		patch(repmat(xt + dt*[0,1,1,0]',[1,2]),[ylim(1) - bwy*[0,0,1,1];ylim(2) + bwy*[0,0,1,1]]','k','clipping','off')
		if fs > 0
			text(xt,ylim(1) - 1.2*bwy,deg2dms(xt,dlon,dec),'FontSize',fs,'FontWeight',fw, ...
				'HorizontalAlignment','center','VerticalAlignment','top');
		end
	end

	ytick = (ddy*ceil(ylim(1)/ddy)):ddy:ylim(2);
	for yt = ytick(1:2:end)
		dt = ddy - max(0,yt + ddy - ylim(2));
		patch([xlim(1) - bwx*[0,0,1,1];xlim(2) + bwx*[0,0,1,1]]',repmat(yt + dt*[0,1,1,0]',[1,2]),'k','clipping','off')
		if fs > 0
			text(xlim(1) - 1.2*bwx,yt,deg2dms(yt,dlat,dec),'FontSize',fs,'FontWeight',fw, ...
				'HorizontalAlignment','center','VerticalAlignment','bottom','rotation',90);
		end
	end
end

% -------------------------------------------------------------------------
% --- Scales legend
if scale
	%wsc = diff(xlim)*0.01;
	wsc = bw0;
	xsc = xlim(2) + wsc*2 + bwx;

	if wmark
		cmap = watermark(cmap,wmark);
	end

	% -- elevation scale (colorbar)
	zscale = linspace(zmin,zmax,length(cmap));
	yscale = linspace(0,diff(ylim)/2,length(cmap));
	ysc = ylim(1);
	ddz = dtick(dz*max(0.5*xyr*diff(xlim)/yscale(end),1));
	ztick = (ddz*ceil(zscale(1)/ddz)):ddz:zscale(end);
	patch(xsc + repmat(wsc*[-1;1;1;-1],[1,length(cmap)]), ...
		ysc + [repmat(yscale,[2,1]);repmat(yscale + diff(yscale(1:2)),[2,1])], ...
		repmat(zscale,[4,1]), ...
		'EdgeColor','flat','LineWidth',.1,'FaceColor','flat','clipping','off')
	colormap(cmap)
	caxis([zmin,zmax]);
	patch(xsc + wsc*[-1,1,1,-1],ysc + yscale(end)*[0,0,1,1],'k','FaceColor','none','Clipping','off')
	text(xsc + 2*wsc + zeros(size(ztick)),ysc + (ztick - zscale(1))*0.5*diff(ylim)/diff(zscale([1,end])),num2str(ztick'), ...
		'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',6)
	% indicates min and max Z values
	text(xsc,ysc - bwy/2,sprintf('%g m',roundsd(zlim(1),3)),'FontWeight','bold', ...
		'HorizontalAlignment','left','VerticalAlignment','top','FontSize',6)
	text(xsc,ysc + .5*diff(ylim) + bwy/2,sprintf('%g m',roundsd(zlim(2),3)),'FontWeight','bold', ...
		'HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',6)
	
	% -- distance scale (in km)
	if dms
		fsc = degkm;
	else
		fsc = zratio/1e3;
	end
	dkm = dtick(diff(ylim)*fsc);
	ysc = ylim(2) - 0.5*dkm/fsc;
	patch(xsc + wsc*[-1,-1,0,0],ysc + dkm*0.5*[-1,1,1,-1]/fsc,'k','FaceColor',grey,'clipping','off')
	if dkm > 1
		skm = sprintf('%g km',dkm);
	else
		skm = sprintf('%g m',dkm*1000);
	end
	text(xsc,ysc,skm,'rotation',-90,'HorizontalAlignment','center','VerticalAlignment','bottom', ...
			'Color',grey,'FontWeight','bold','FontSize',6)
end


if nargout > 0
	h = hh;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = nmedian(x,n)
%NMEDIAN Generalized median filter
%	NMEDIAN(X,N) sorts elemets of X and returns N-th value (N normalized).
%	So:
%	   N = 0 is minimum value
%	   N = 0.5 is median value
%	   N = 1 is maximum value

if nargin < 2
	n = 0.5;
end
y = sort(x(:));
y = interp1(sort(y),n*(length(y)-1) + 1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dd = dtick(dlim,deg)
%DTICK Tick intervals

if nargin < 2
	deg = 0;
end

if deg && dlim <= 2/60
	% less than 2 minutes: base 36
	m = 10^floor(log10(dlim*36))/36;
elseif deg && dlim <= 2
	% less than 2 degrees: base 6
	m = 10^floor(log10(dlim*6))/6;
else
	% more than few degrees or not degrees: decimal rules
	m = 10^floor(log10(dlim));
end
p = ceil(dlim/m);
if p <= 1
	dd = .1*m;
elseif p == 2
	dd = .2*m;
elseif p <= 5
	dd = .5*m;
else
	dd = m;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = deg2dms(x,ll,dec)
%DEG2DMS Degree/minute/second display

if dec
	s = sprintf('%7.7g',x);
else
	xa = abs(x) + 1/360000;
	%sd = sprintf('%d%c',floor(xa),176);	% ASCII char 176 is the degree sign
	sd = sprintf('%d?',floor(xa));
	sm = '';
	ss = '';
	if mod(x,1)
		sm = sprintf('%02d''',floor(mod(60*xa,60)));
		sa = floor(mod(3600*xa,60));
		if sa
			ss = sprintf('%02d"',sa);
		else
			if strcmp(sm,'00''')
				sm = '';
			end
		end
	end
	s = [sd,sm,ss,ll{1+int8(x<0)}];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = fillgap(x,y,z)
% GRIDDATA is not efficient for large arrays, but has great advantage to be
% included in Matlab's core functions! To optimize interpolation, we
% reduce the amount of relevant data by building a mask of all surrounding
% pixels of novalue areas... playing with linear index!

sz = size(z);
k = find(isnan(z));
k(k == 1 | k == numel(z)) = []; % removes first and last index (if exist)
if ~isempty(k)
	[xx,yy] = meshgrid(x,y);
	mask = zeros(sz,'int8');
	k2 = ind90(sz,k); % k2 is linear index in the row order
	% sets to 1 every previous and next index, both in column and row order
	mask([k-1;k+1;ind90(fliplr(sz),[k2-1;k2+1])]) = 1; 
	mask(k) = 0; % removes the novalue index
	kb = find(mask); % keeps only border values
	z(k) = griddata(xx(kb),yy(kb),z(kb),xx(k),yy(k));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function k2 = ind90(sz,k)

[i,j] = ind2sub(sz,k);
k2 = sub2ind(fliplr(sz),j,i); % switched i and j: k2 is linear index in row order


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function k = islake(z)
% ISLAKE mask of zero gradient on 3x3 tiles
% We use diff matrix in row and column directions, and shift it to build
% a single vectorized test of surrounding pixels. To do this we must
% concatenate unit vectors in different combinations...

dx = diff(z,1,2);	% differences in X direction
dy = diff(z,1,1);	% differences in Y direction
u1 = ones(size(z,1),1);	% row unit vector 
u2 = ones(1,size(z,2));	% column unit vector
u2r = u2(2:end);

% index of the tiles center pixel
k = ( ...
	[u2;dy] == 0 & [dy;u2] == 0 & ...
	[u1,dx] == 0 & [dx,u1] == 0 & ...
	[u1,[dx(2:end,:);u2r]] == 0 & [[dx(2:end,:);u2r],u1] == 0 & ...
	[u1,[u2r;dx(1:end-1,:)]] == 0 & [[u2r;dx(1:end-1,:)],u1] == 0 ...
);

% now extends it to surrounding pixels
k(1:end-1,:) = (k(1:end-1,:) | k(2:end,:));
k(2:end,:) = (k(2:end,:) | k(1:end-1,:));
k(:,1:end-1) = (k(:,1:end-1) | k(:,2:end));
k(:,2:end) = (k(:,2:end) | k(:,1:end-1));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = isrgb(x,n)

if nargin < 2
	n = 0;
end
if isnumeric(x) && (n == 1 && all(size(x) == [1,3]) || n == 0 && size(x,2) == 3) ...
		&& all(x(:) >= 0 & x(:) <= 1)
	s = 1;
else
	s = 0;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = isperc(x)

if isnumeric(x) && isscalar(x) && x >= 0 && x <= 100
	s = 1;
else
	s = 0;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = isvec(x,n)

if nargin < 2
	n = 2;
end
if isnumeric(x) && numel(x) == n
	s = 1;
else
	s = 0;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y=roundsd(x,n)

og = 10.^(floor(log10(abs(x)) - n + 1));
y = round(x./og).*og;
y(x==0) = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = watermark(x,n)

if nargin < 2
	n = 2;
end

if n == 0
    y = x;
else
    y = (x/n + 1 - 1/n);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [s,v] = checkparam(arg,nam,func,val)

switch func2str(func)
	case 'isscalar'
		num = 1;
		mes = 'scalar value';
	case 'isperc'
		num = 1;
		mes = 'percentage scalar value';
	case 'isvec'
		num = 1;
		if nargin < 4
			val = 2;
		end
		mes = sprintf('%d-element vector',val);
	case 'isrgb'
		num = 1;
		mes = '[R,G,B] vector with 0.0 to 1.0 values';
	case 'ischar'
		num = 0;
		mes = 'string';
		if nargin > 3
			mes = sprintf('%s (%s)',mes,strjoin(val,' or '));
		end
	otherwise
		num = 1;
		mes = 'value';
end

s = 0;
v = [];
k = find(strcmpi(arg,nam));
if ~isempty(k)
	if (k + 1) <= length(arg) ...
			&& (~num || isnumeric(arg{k+1})) ...
			&& (nargin < 4 && func(arg{k+1}) ...
				|| (nargin > 3 && (strcmp(func2str(func),'ischar') && ismember(arg{k+1},val)) ...
					 || strcmp(func2str(func),'isvec') && func(arg{k+1},val)))
		v = arg{k+1};
		s = 1;
	else
		error('%s option must be followed by a valid %s.',upper(nam),mes)
	end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s=strjoin(c,d)
%STRJOIN Join cell array of strings
%(this is for Matlab versions < 2013a backward compatibility)

if nargin < 2
	d = '';
end
n = numel(c);
ss = cell(2,n);
ss(1,:) = reshape(c,1,n);
ss(2,1:n-1) = {d};
s = [ss{:}];

%% Function 15: Landcolor for dem
function y = landcolor(n)
%LANDCOLOR Land colormap
%
%	Author: Francois Beauducel <beauducel@ipgp.fr>
%	Date: 2012/05/17 11:22:44 $

J = [ ...
0.095678 0.53427 0.21682 
0.15785 0.5979 0.23274 
0.21286 0.64673 0.2514 
0.26411 0.68789 0.27268 
0.32959 0.72416 0.31308 
0.39794 0.75695 0.36038 
0.46153 0.7871 0.40624 
0.52108 0.81516 0.45135 
0.57702 0.84152 0.49547 
0.62973 0.86645 0.53891 
0.67946 0.89016 0.58187 
0.72647 0.91282 0.62427 
0.77095 0.93455 0.66619 
0.81306 0.95546 0.70772 
0.85292 0.97563 0.7489 
0.89066 0.99514 0.78976 
0.88379 0.98595 0.77038 
0.86389 0.96758 0.73236 
0.84615 0.94972 0.69623 
0.8303 0.93233 0.66186 
0.81612 0.91536 0.6291 
0.80341 0.8988 0.59784 
0.79201 0.8826 0.56795 
0.78191 0.86676 0.53946 
0.7729 0.85123 0.51224 
0.76479 0.83602 0.48615 
0.75747 0.8211 0.46111 
0.75084 0.80645 0.43704 
0.74506 0.79206 0.41414 
0.73981 0.77792 0.39211 
0.73501 0.76401 0.37089 
0.73068 0.75033 0.35052 
0.72683 0.73685 0.33106 
0.72042 0.72074 0.31228 
0.71032 0.70085 0.29417 
0.69761 0.67821 0.27694 
0.68489 0.65558 0.26026 
0.67235 0.63313 0.24418 
0.65997 0.61082 0.22889 
0.64775 0.58874 0.21406 
0.63568 0.56689 0.19983 
0.62376 0.54527 0.18622 
0.61197 0.52391 0.17299 
0.60033 0.50283 0.16046 
0.58881 0.48203 0.14832 
0.57742 0.46151 0.13667 
0.56616 0.44133 0.12555 
0.55502 0.4214 0.11472 
0.54398 0.4019 0.10456 
0.53306 0.38266 0.094633 
0.52226 0.36382 0.085242 
0.51155 0.3453 0.076179 
0.50095 0.32714 0.067515 
0.49045 0.30938 0.059259 
0.48005 0.29193 0.051294 
0.46973 0.27495 0.043796 
0.45951 0.25823 0.0365 
0.44938 0.24206 0.029715 
0.43934 0.22609 0.023063 
0.42938 0.21074 0.016949 
0.41951 0.19556 0.010917 
0.40971 0.18105 0.0054326 
0.4 0.16667 0 
];

l = length(J);
if nargin < 1
	n = 256;
end
y = interp1(1:l,J,linspace(1,l,n),'*linear');