function varargout = SH_filter(varargin)
% SH_FILTER MATLAB code for SH_filter.fig
%      SH_FILTER, by itself, creates a new SH_FILTER or raises the existing
%      singleton*.
%
%      H = SH_FILTER returns the handle to a new SH_FILTER or the handle to
%      the existing singleton*.
%
%      SH_FILTER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SH_FILTER.M with the given input arguments.
%
%      SH_FILTER('Property','Value',...) creates a new SH_FILTER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SH_filter_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SH_filter_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SH_filter

% Last Modified by GUIDE v2.5 31-Dec-2014 09:28:51

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SH_filter_OpeningFcn, ...
                   'gui_OutputFcn',  @SH_filter_OutputFcn, ...
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


% --- Executes just before SH_filter is made visible.
function SH_filter_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SH_filter (see VARARGIN)

% Choose default command line output for SH_filter
handles.output = hObject;

file=('sh_projected.txt');
if exist(file,'file')==1 || exist(file,'file')==2
fid=fopen(file,'r');
i=1;
while ~feof(fid)
    A=textscan(fid,'%u %u %u %u/','delimiter',' ','HeaderLines',1);
    
    for j=1:4 %number of columns (default from shorelines.txt)
        data(i,j) = A(1,j);
    end
    i=i+1;
end
fclose(fid);

%extract data from table
mdp(:,1)=cell2mat(data(:,1)); %shoreline
mdp(:,2)=cell2mat(data(:,2)); %X
mdp(:,3)=cell2mat(data(:,3)); %Z
mdp(:,4)=cell2mat(data(:,4)); %Ze

end
 
mdp=double(mdp);

% set shoreline levels
shr=mdp(:,1);%shore level

%save handles
handles.mdp=mdp; %values
handles.shr=shr; %shore
sh_lev=num2cell(unique(shr)); %types of levels

%plot all shorelines
plot(mdp(:,2),mdp(:,3),'o','MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7],'MarkerSize',6)


%populate listbox
set(handles.listbox_level,'String',sh_lev);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SH_filter wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SH_filter_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in PL2.
function PL2_Callback(hObject, eventdata, handles)
% hObject    handle to PL2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


level_selection=str2num(handles.level);
mdp=handles.mdp;
tee=handles.edit1;%%%% Tolerance input

% interpolation eq
x_min=min(mdp(mdp(:,1)==level_selection,2))-tee;
x_max=max(mdp(mdp(:,1)==level_selection,2))+tee;


%C(:,1)=mdp(:,2); %higway terraces X axis
%C(:,2)=mdp(:,3); %%higway terraces Y axis

C(:,1)=mdp(mdp(:,1)==level_selection,2);
C(:,2)=mdp(mdp(:,1)==level_selection,3);


CC=sortrows(C);

NEDX=CC(:,1);
NEDY=CC(:,2);

%tee=1000;
Cix=x_min:tee:x_max; %size of interpolation
handles.Cix=Cix;
try
Ci=interp1(NEDX,NEDY,Cix);%interpolation
Ci=smooth(Ci); %optional optimization
xs1=handles.xs1;%main plot
cla(xs1)
axes(xs1);
yy =(Ci);
handles.Ci=Ci;
catch   
msgbox('There are repeated shoreline values, delete repeated values form the table and try to filter again')    
return
end    

hold on
plot(mdp(:,2),mdp(:,3),'o','MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7],'MarkerSize',6)
plot(mdp(mdp(:,1)==level_selection,2),mdp(mdp(:,1)==level_selection,3),'o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',6)
plot(x_min:tee:x_max,yy,'-r')
clear Ci tee

guidata(hObject, handles);


% --- Executes on button press in PL4.
function PL3_Callback(hObject, eventdata, handles) %set filters
% hObject    handle to PL4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
Ci=handles.Ci;%interpolation from previous routine
yy = (Ci);

tee=handles.edit1;%interpolation spacing
plus=handles.edit2; %standard deviation
mdp=handles.mdp; %shorelines loaded
xs1=handles.xs1; %main plot
level_selection=str2num(handles.level);
% ranges
x_min=min(mdp(mdp(:,1)==level_selection,2))-tee;
x_max=max(mdp(mdp(:,1)==level_selection,2))+tee;
yy_mp=yy+plus; % upper bound
yy_np=yy-plus; % lower bound

cla(xs1)%clear plot area of GUI
axes(xs1);
%plot shorelines and fits
hold on
plot(mdp(:,2),mdp(:,3),'o','MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7],'MarkerSize',6);
plot(mdp(mdp(:,1)==level_selection,2),mdp(mdp(:,1)==level_selection,3),'o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',6);
plot(x_min:tee:x_max,yy,'-b');
plot(x_min:tee:x_max,yy_mp,'--r');
plot(x_min:tee:x_max,yy_np,'--r');

%call handles of miniplot
bord_up(:,1)=x_min:tee:x_max; %xvalues from interpolation upper bound
bord_up(:,2)=yy_mp;
bord_down(:,1)=x_min:tee:x_max; %xvalues from interpolation lower bound
bord_down(:,2)=yy_np;

shores(:,1)=mdp(mdp(:,1)==level_selection,2);
shores(:,2)=mdp(mdp(:,1)==level_selection,3);
kk=1;
%shorelines excluded above upper bound

x=bord_up(:,1);
y=bord_up(:,2);
x2=bord_down(:,1);
y2=bord_down(:,2);

%use denecke function to eaxtract points above or below interp1 matrix (Ci)
for ik=1:numel(shores(:,1));
point=shores(ik,:);        
den(kk,1) = denecke(x,y,point);
den2(kk,1) = denecke(x2,y2,point);
kk=kk+1;
end

%prepare plot inputs
hold on
shores3=horzcat(shores, den2);
shores2=horzcat(shores, den);
%upper shores excluded red
X1=shores2(shores2(:,3)==1,1);
Y1=shores2(shores2(:,3)==1,2);

%lower shores excluded
X2=shores3(shores3(:,3)==2,1);
Y2=shores3(shores3(:,3)==2,2);

%plot excluded and incuded
plot(X1,Y1,'o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',6)
plot(X2,Y2,'o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',6)

%prepare barplot
ax2=handles.xs2; 
axes(ax2);
cla(ax2);%clear plot area of GUI

%mini bar plot
labelsx={'outliers', 'included'};
Yb=[numel(X1)+numel(X2) numel(shores(:,1))-(numel(X1)+numel(X2))];
barh(Yb); 
title('% of shoreline-angles interpolated');
set(gca,'yticklabel',labelsx);

percent1=((numel(X1)+numel(X2))*100)/numel(shores(:,1));
percent2=100-percent1;
per=('%');
t1=sprintf('%.f%s',percent1,per);
t2=sprintf('%.f%s',percent2,per);
hold on
text(numel(shores(:,1)),1,t1);
text(numel(shores(:,1)),2,t2);
%set output handles for saving function
catch
msgbox('try interpolating X befor defining the Y-range')  
return
end    
handles.outliers=vertcat(horzcat(X1, Y1),horzcat(X2, Y2));
guidata(hObject, handles);%save handles


function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
handles.output = hObject;
handles.edit1=str2double(get(hObject,'String')); %returns contents of edit1 as a double

guidata(hObject, handles);%save handles
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


% --- Executes on selection change in listbox_level.
function listbox_level_Callback(hObject, eventdata, handles) %shoreline levels list
% hObject    handle to listbox_level (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_level contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_level
handles = guidata(hObject);
contents = cellstr(get(hObject,'String')); %returns listbox_level contents as cell array
level_selection=(contents{get(hObject,'Value')}); %returns selected item from listbox
handles.level=level_selection;


guidata(hObject, handles);%save handles


% --- Executes during object creation, after setting all properties.
function listbox_level_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_level (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in PL1.
function PL1_Callback(hObject, eventdata, handles)
% hObject    handle to PL1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = guidata(hObject);

level_selection=str2num(handles.level);
mdp=handles.mdp;

xs1=handles.xs1;%main plot

cla(xs1)

axes(xs1)
hold on
plot(mdp(:,2),mdp(:,3),'o','MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7],'MarkerSize',6);
plot(mdp(mdp(:,1)==level_selection,2),mdp(mdp(:,1)==level_selection,3),'o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',6);
%plot(mdp(mdp(:,1)==level_selection,2),mdp(mdp(:,1)==level_selection,3),'ok')

xs2=handles.xs2; %miniplot hist

ht=numel(mdp(mdp(:,1)==level_selection,3));

axes(xs2)

guidata(hObject, handles);

% --- Executes on button press in PL4.
function PL4_Callback(hObject, eventdata, handles)
% hObject    handle to PL4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%dp=handles.dp;
proj_out=handles.outliers; %call filtered outliers
levsel=str2num(handles.level);
yinterp=(handles.Ci);
xinterp=(handles.Cix');

sh_sel(1:numel(proj_out(:,1)))=levsel;

%fil(:,1)=sh_sel; %shoreline number
%proj(:,2)=proj_out(:,3); %distance profile
fil(:,1)=proj_out(:,1); %x
fil(:,2)=proj_out(:,2); %y

filename2=sprintf('Outliers_level_%.f.txt',levsel);

titles={'X', 'Y'}; %table headers
dlmcell(filename2,titles); %add headers for table
proj2=num2cell(fil); %transform in cell type
dlmcell(filename2,proj2,'-a'); %save data to table

filename=sprintf('Interpolation_level_%.f.txt',levsel);
%save( filename,'xinterp','yinterp','-ASCII')
titles2={'X', 'Y'}; %table headers
dlmcell(filename,titles2); %add headers for table

tabf1=horzcat(xinterp,yinterp);
proj3=num2cell(tabf1); %transform in cell type
dlmcell(filename,proj3,'-a'); %save data to table

function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text

handles.output = hObject;
handles.edit2=str2double(get(hObject,'String')); %returns contents of edit2 as a double

guidata(hObject, handles);%save handles


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% utilitary functions
%% Nested Function 1 denecke
function den = denecke(x,y,point)
% detect point above and below an interp1 curve
%Jara-Mu?oz 2015
% used for detect points above or below interp1 outputs
% This function takes as arguments two equal-length vectors x and y. point is 
% a vector of length 2 (any other point). The function determines if the point
%lies above (+1), on (0) or below (-1) the function y = f(x). 
%% tip, interpoation must be a litle bit longer point must lie inside interp1

p=point;
lx = max(x(x < p(1)));
ly = y(x == lx); %y to max x

if ly<p(2) && lx<p(1);
den=1;
elseif ly>p(2); 
den=2;    
else    
 den=0;   
end

%% %%%%% Nested function 2 dlmcell
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
