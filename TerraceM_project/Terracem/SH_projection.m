function varargout = SH_projection(varargin)
% SH_PROJECTION MATLAB code for SH_projection.fig
%      SH_PROJECTION, by itself, creates a new SH_PROJECTION or raises the existing
%      singleton*.
%
%      H = SH_PROJECTION returns the handle to a new SH_PROJECTION or the handle to
%      the existing singleton*.
%
%      SH_PROJECTION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SH_PROJECTION.M with the given input arguments.
%
%      SH_PROJECTION('Property','Value',...) creates a new SH_PROJECTION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SH_projection_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SH_projection_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SH_projection

% Last Modified by GUIDE v2.5 07-Apr-2015 21:29:31

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SH_projection_OpeningFcn, ...
                   'gui_OutputFcn',  @SH_projection_OutputFcn, ...
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





% --- Executes just before SH_projection is made visible.
function SH_projection_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SH_projection (see VARARGIN)

% Choose default command line output for SH_projection
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SH_projection wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SH_projection_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



% --- Executes on button press in export.
function export_Callback(hObject, eventdata, handles)
% hObject    handle to export (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)





% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

handles = guidata(hObject); %define handles

% Hint: place code in OpeningFcn to populate axes1
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
%extract data from table
nim=numel(data(:,10));
u=1;

for k=1:nim
md1(u,1)=data{k,4}; %x
md1(u,2)=data{k,5}; %Y
md1(u,3)=data{k,7}; %z
md1(u,4)=data{k,8}; %ze
md1(u,5)=data{k,3}; %shore
u=1+u;
end

else
%no shorelines option    
 msgbox('No shorelines have been processed')  
 return
 
end

% set shoreline levels
shore=md1(:,5);

%save handles
handles.md1=md1;
handles.shore=shore;

%plot all shorelines
hold on
scatter(md1(:,1),md1(:,2),75,md1(:,3),'filled'); %plot again shorelines
%g1=scatter(x,y,75,z,'filled');
%col=colorbar;


axis equal
box on
ylabel('UTM   N')
xlabel('UTM   E')

guidata(hObject, handles); %save handles

% --- Executes on button press in radiobutton11.
function radiobutton11_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject); %define handles
% Hint: get(hObject,'Value') returns toggle state of radiobutton11

if (get(hObject,'Value') == get(hObject,'Max'))
[xl,yl]=ginput(2);
l(:,1)=xl;
l(:,2)=yl;

%plot profile
hold on
p1=plot(l(:,1),l(:,2),'-k'); %plot of profile
p2=plot(l(1,1),l(1,2),'xk'); %point1
p3=plot(l(2,1),l(2,2),'xk'); %point2

%profile geometry
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
text(min(l(:,1))+((max(l(:,1))-min(l(:,1)))/2)+((max(l(:,1))-min(l(:,1)))/10), min(l(:,2))+(max(l(:,2))-min(l(:,2)))/2,sprintf('%s km',str1),'FontSize',12)

%call handles
md1=handles.md1;

%rotate and reproject data
for i=1:length(md1(:,1))
    b2=md1(i,2)-m2*md1(i,1);
    Pnx=(b2-b)/(m-m2);
    Pny=m*Pnx+b;
    Pn=[Pnx Pny];
    dp(i,1)=Pn(1,1); %xp
    dp(i,2)=Pn(1,2); %yp
    dp(i,3)=((l(1,2)-dp(i,2))/sin(anr)); % distance in km
    dp(i,4)=md1(i,3);%Z
    dp(i,5)=md1(i,4);%ZE
    dp(i,6)=md1(i,5);%shore
end

%save handles
handles.dp=dp;
handles.dx=dx;
handles.dy=dy;
handles.dkm=dkm;

guidata(hObject, handles); %save handles
end    

% --- Executes on button press in radiobutton12.
function radiobutton12_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject); %define handles
% Hint: get(hObject,'Value') returns toggle state of radiobutton12

if (get(hObject,'Value') == get(hObject,'Max'))
%load shapefile    
[Filename,PathName]=uigetfile('*.shp','Select shapefile');

    if isequal(Filename,0)
    disp('User selection Canceled')
    return
    
    else
    S4 = shaperead(Filename);    
    l(:,1)=[S4.X];
    l(:,2)=[S4.Y];     
    end
    
%plot profile    
hold on
p1=plot(l(:,1),l(:,2),'-k'); %plot of profile
p2=plot(l(1,1),l(1,2),'xk'); %point1
p3=plot(l(2,1),l(2,2),'xk'); %point2

%define profile geometry
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

%call handles
md1=handles.md1;

%rotate and reproject data
for i=1:length(md1(:,1))
    b2=md1(i,2)-m2*md1(i,1);
    Pnx=(b2-b)/(m-m2);
    Pny=m*Pnx+b;
    Pn=[Pnx Pny];
    dp(i,1)=Pn(1,1); %xp
    dp(i,2)=Pn(1,2); %yp
    dp(i,3)=((l(1,2)-dp(i,2))/sin(anr)); % distance in km
    dp(i,4)=md1(i,3); %Z
    dp(i,5)=md1(i,4); %Ze
    dp(i,6)=md1(i,5);%shore
end

%save handles
handles.dp=dp;
handles.dx=dx;
handles.dy=dy;
handles.dkm=dkm;
end
guidata(hObject, handles); %save handles
   
% --- Executes on button press in B1. (plot in axes 2)
function B1_Callback(hObject, eventdata, handles)
% hObject    handle to B1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject); %define handles

a1=handles.axes3;
axes(a1)


if isempty('handles.dp')==0
dp=handles.dp;


%call handles
shore=handles.shore;

dkm=handles.dkm;
dx=handles.dx;
dy=handles.dy;

hold on
box on
%loop shorelines colors
%check if only one shore level
for i=min(shore):max(shore)       
   plot(a1,dp(dp(:,6)==i,3)/1000,dp(dp(:,6)==i,4),'.','color',[1/i 0 1/i],'markersize',20)%color for shoreine levels
   j(i,1)=i;
end

%cartesian plane conditions for profiles
if dx>0 && dy>0   
xlim([-dkm/1000 0])

proj_filt1=dp(dp(:,3)<=0,:);
proj_filt=proj_filt1(proj_filt1(:,3)>=-dkm,:);
a=1; %a is a parameter for data export

elseif dx<0 && dy<0 
xlim([0 dkm/1000]) 

proj_filt1=dp(dp(:,3)>=0,:);
proj_filt=proj_filt1(proj_filt1(:,3)<dkm,:);
a=3;

elseif dx>0 && dy<0  
xlim([-dkm/1000 0]) 

proj_filt1=dp(dp(:,3)<=0,:);
proj_filt=proj_filt1(proj_filt1(:,3)>=-dkm,:);
a=4;

else
xlim([0 dkm/1000]) 
    
proj_filt1=dp(dp(:,3)>=0,:);    
proj_filt=proj_filt1(proj_filt1(:,3)<=dkm,:);  
a=2;

end
%legend
%leg=legend(num2str(j),'location','EastOutside'); %legend specifications
%tit = get(leg,'title');
%set(tit,'string','Shorelines'); %legend title

xg=dp(:,3)./1000; %km scale
errorbar(xg,dp(:,4),(dp(:,5)),'ok') %plot errorbars

ylabel('Elevation (m)')
xlabel('Distance along profile (km)')

handles.proj_filt=proj_filt; %constrained and filtered shorelines
else

msgbox('No profiles load or drawn yet')    
end    
guidata(hObject, handles); %save handles

% --- Executes on button press in B2.(Save projected points)
function B2_Callback(hObject, eventdata, handles)
% hObject    handle to B2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%dp=handles.dp;
proj_filt=handles.proj_filt; %call projected shorelines

proj(:,1)=proj_filt(:,6); %shoreline number
proj(:,2)=proj_filt(:,3); %distance profile
proj(:,3)=proj_filt(:,4); %Z
proj(:,4)=proj_filt(:,5); %Ze

titles={'shoreline_number', 'distance_along_profile', 'shoreline_elevation', 'error'}; %table headers
dlmcell('sh_projected.txt',titles); %add headers for table
proj2=num2cell(proj); %transform in cell type
dlmcell('sh_projected.txt',proj2,'-a'); %save data to table


% --- Executes on button press in B3.
function B3_Callback(hObject, eventdata, handles) %clear plots
% hObject    handle to B3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject); %define handles

a1=handles.axes1;
a3=handles.axes3;
cla(a3)% delete contents
cla(a1)% delete contents


axes(a1)
md1=handles.md1;
scatter(md1(:,1),md1(:,2),75,md1(:,3),'filled'); %plot again shorelines

guidata(hObject, handles); %save handles

%% %%%%% Nested function 1 dlmcell
function dlmcell(file,cell_array,varargin)
%% <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> %%
% <><><><><>     dlmcell - Write Cell Array to Text File      <><><><><> %
% <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> %
%                                                 Version:    01.06.2010 %
%                                                     (c) Roland Pfister %
%                                             roland_pfister@t-online.de %
%                        ...with many thanks to George Papazafeiropoulos %
%                        for his corrections and improvements.           %
% 1. Synopsis                                                            %
%                                                                        %
% A single cell array is written to an output file. Cells may consist of %
% any combination of (a) numbers, (b) letters, or (c) words. The inputs  %
% are as follows:                                                        %
%                                                                        %
%       - file       The output filename (string).                       %
%       - cell_array The cell array to be written.                       %
%       - delimiter  Delimiter symbol, e.g. ',' (optional;               %
%                    default: tab ('\t'}).                               %
%       - append     '-a' for appending the content to the               %
%                    output file (optional).                             %
%                                                                        %
% 2. Example                                                             %
%                                                                        %
%         mycell = {'Numbers', 'Letters', 'Words','More Words'; ...      %
%                    1, 'A', 'Apple', {'Apricot'}; ...                   %
%                    2, 'B', 'Banana', {'Blueberry'}; ...                %
%                    3, 'C', 'Cherry', {'Cranberry'}; };                 %
%         dlmcell('mytext.txt',mycell);                                  %
%                                                                        %
% <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> %



%% Check input arguments
if nargin < 2
    disp('Error - Give at least two input arguments!');
    return;
elseif nargin > 4
    disp('Error - Do not give more than 4 input arguments!');
    return;
end
if ~ischar(file)
    disp(['Error - File input has to be a string (e.g. ' ...
        char(39) 'output.txt' char(39) '!']);
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

%% Open output file and prepare output array.
output_file = fopen(file,append);
output = cell(size(cell_array,1),size(cell_array,2));

%% Evaluate and write input array.
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
            %
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
