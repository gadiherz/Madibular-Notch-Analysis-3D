% -------------------------------------------------------------------------
% Mandibular upper-ramus workflow (Mandibles_WF)
%
% Authors: Gadi Herzlinger; Uzy Smilansky
% Copyright (c) 2024 Gadi Herzlinger and Uzy Smilansky. All rights reserved.
% Documentation assistance: ChatGPT (comments and documentation only).
% -------------------------------------------------------------------------


function varargout = SelectSide4Pos(varargin)
% SELECTSIDE4POS MATLAB code for SelectSide4Pos.fig
%      SELECTSIDE4POS, by itself, creates a new SELECTSIDE4POS or raises the existing
%      singleton*.
%
%      H = SELECTSIDE4POS returns the handle to a new SELECTSIDE4POS or the handle to
%      the existing singleton*.
%
%      SELECTSIDE4POS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SELECTSIDE4POS.M with the given input arguments.
%
%      SELECTSIDE4POS('Property','Value',...) creates a new SELECTSIDE4POS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SelectSide4Pos_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SelectSide4Pos_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SelectSide4Pos

% Last Modified by GUIDE v2.5 08-Sep-2020 20:58:52

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SelectSide4Pos_OpeningFcn, ...
                   'gui_OutputFcn',  @SelectSide4Pos_OutputFcn, ...
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


% --- Executes just before SelectSide4Pos is made visible.
function SelectSide4Pos_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SelectSide4Pos (see VARARGIN)

RotPCNeg = varargin{1};
PeaksNeg = varargin{2};
angNeg = varargin{3};
RotPCPos = varargin{4};
PeaksPos = varargin{5};
angPos = varargin{6};
RotPCComb = varargin{7};
PeaksComb = varargin{8};
angComb = varargin{9};

if angNeg==0
    set(handles.Left, 'Enable', 'off');
end
if angPos==0
    set(handles.Right, 'Enable', 'off');
end
if angComb==0
    set(handles.Combined, 'Enable', 'off');
end

RotPCNeg = ([cos(-angNeg), -sin(-angNeg); sin(-angNeg), cos(-angNeg)] *RotPCNeg')';
hold(handles.NegAx,'on');
plot(handles.NegAx, RotPCNeg(:,1), RotPCNeg(:,2));
%axis equal
if ~isempty(PeaksNeg)
    scatter(handles.NegAx, RotPCNeg(PeaksNeg(find(PeaksNeg(:,3)>0),1),1),RotPCNeg(PeaksNeg(find(PeaksNeg(:,3)>0),1),2),...
        'r','o');
    scatter(handles.NegAx, RotPCNeg(PeaksNeg(find(PeaksNeg(:,3)<0),1),1),RotPCNeg(PeaksNeg(find(PeaksNeg(:,3)<0),1),2),...
        'g','o');
end

RotPCComb = ([cos(-angComb), -sin(-angComb); sin(-angComb), cos(-angComb)] *RotPCComb')';
hold(handles.CombAx,'on');
plot(handles.CombAx, RotPCComb(:,1), RotPCComb(:,2));
%axis equal
if ~isempty(PeaksComb)
    scatter(handles.CombAx, RotPCComb(PeaksComb(find(PeaksComb(:,3)>0),1),1),RotPCComb(PeaksComb(find(PeaksComb(:,3)>0),1),2),...
        'r','o');
    scatter(handles.CombAx, RotPCComb(PeaksComb(find(PeaksComb(:,3)<0),1),1),RotPCComb(PeaksComb(find(PeaksComb(:,3)<0),1),2),...
        'g','o');
end

RotPCPos = ([cos(-angPos), -sin(-angPos); sin(-angPos), cos(-angPos)] *RotPCPos')';
hold(handles.PosAx,'on');
plot(handles.PosAx, RotPCPos(:,1), RotPCPos(:,2));
%axis equal
if ~isempty(PeaksPos)
    scatter(handles.PosAx, RotPCPos(PeaksPos(find(PeaksPos(:,3)>0),1),1),RotPCPos(PeaksPos(find(PeaksPos(:,3)>0),1),2),...
        'r','o');
    scatter(handles.PosAx, RotPCPos(PeaksPos(find(PeaksPos(:,3)<0),1),1),RotPCPos(PeaksPos(find(PeaksPos(:,3)<0),1),2),...
        'g','o');
end

XLims = [get(handles.NegAx,'XLim'); get(handles.CombAx,'XLim');...
    get(handles.PosAx,'XLim')];
XLims = [min(XLims(:,1)), max(XLims(:,2))];
YLims = [get(handles.NegAx,'YLim'); get(handles.CombAx,'YLim');...
    get(handles.PosAx,'YLim')];
YLims = [min(YLims(:,1)), max(YLims(:,2))];

set(handles.NegAx,'XLim',XLims);
set(handles.NegAx,'YLim',YLims);
set(handles.CombAx,'XLim',XLims);
set(handles.CombAx,'YLim',YLims);
set(handles.PosAx,'XLim',XLims);
set(handles.PosAx,'YLim',YLims);


% Choose default command line output for SelectSide4Pos
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SelectSide4Pos wait for user response (see UIRESUME)
 uiwait(handles.figure1);
 


% --- Outputs from this function are returned to the command line.
function varargout = SelectSide4Pos_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
delete(hObject)


% --- Executes on button press in Left.
function Left_Callback(hObject, eventdata, handles)
% hObject    handle to Left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.output = 1;
guidata(hObject, handles);
uiresume(handles.figure1)

% --- Executes on button press in Combined.
function Combined_Callback(hObject, eventdata, handles)
% hObject    handle to Combined (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.output = 2;
guidata(hObject, handles);
uiresume(handles.figure1)

% --- Executes on button press in Right.
function Right_Callback(hObject, eventdata, handles)
% hObject    handle to Right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.output = 3;
guidata(hObject, handles);
uiresume(handles.figure1)


% --- Executes on button press in Cancel.
function Cancel_Callback(hObject, eventdata, handles)
% hObject    handle to Cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.output = 4;
guidata(hObject, handles);
uiresume(handles.figure1)
