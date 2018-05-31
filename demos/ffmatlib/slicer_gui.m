%slicer_gui.m Graphical user interface to slice FreeFem++ simulations
%
% Author: Chloros2 <chloros2@gmx.de>
% Created: 2018-05-19
%
% Copyright (C) 2018 Chloros2 <chloros2@gmx.de>
% 
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see
% <https://www.gnu.org/licenses/>.
%
function slicer_gui(bddatafile,tetdatafile)

    %listbox slicing normals preselection
    slicingNormals={'0;0;1';'0;1;0';'1;0;0';'0;0;-1';'0;-1;0';'-1;0;0'; ...
                    '0;1;1';'0;1;-1';'1;1;0';'1;-1;0';'1;0;-1'; ...
                    '-1;1;-1';'-1;1;1';'1;1;1';'1;-0.5;-0.75'};
    %slicing point
    handles.O=[0.5 0.5 0.2];
    %startup slicing plane normal
    handles.N=[0 0 1];

    %slider and colorbar control variables
    handles.Tspread=[0.5 2.0 1.0];
    %initial guess at startup, will be updated later on
    handles.Tshift=[-1.0 1.0 0.0];
    %hold min and max temperature of the acutal slice used for colorbar
    handles.Tmax=1;
    handles.Tmin=0;
    %used for color scale if "FIX" check box is activated in order to freeze
    %the color scale limits
    handles.Tminfrz=handles.Tmin;
    handles.Tmaxfrz=handles.Tmax;

    handles.haveSlicingdata=false;

    winWidth=301;
    winHeight=440;
    D1=75;
    D2=-5;
    handles.f=figure('Name','GUI-Slicer','NumberTitle','off','MenuBar', ...
                     'none','Position', [20 50 winWidth winHeight], ...
                     'Resize','off','Visible','off');

    uicontrol('Style', 'pushbutton', 'String', ...
              'Close All', 'Position', ...
              [199 29 80 23], ...
              'Callback', @pushbuttonExit_Callback);
    uicontrol('Style', 'pushbutton', 'String', ...
              'Draw', 'Position', ...
              [20 29 80 23], ...
              'Callback', @pushbuttonCalc_Callback);
    uicontrol('Style', 'pushbutton', 'String', ...
              'Reset Color', 'Position', ...
              [19 (D1+184) 80 23], ...
              'Callback', @pushbuttonTreset_Callback);

    uicontrol('style','text','String', 'Slicing Point:','position', ...
              [20 (D2+123) 100 21], ...
              'visible','on','HorizontalAlignment','left', ...
              'BackgroundColor',[0.95,0.95,0.95]);
    handles.editx = uicontrol('Style','edit',...
                              'Position',[20 (D2+76) 51 22], ...
                              'String',num2str(handles.O(1)), ...
                              'BackgroundColor',[1,1,1], ...
                              'CallBack',@editX_Callback);
    uicontrol('style','text','String', 'x','position', ...
              [41 (D2+99) 31 21], ...
              'visible','on','HorizontalAlignment','left', ...
              'BackgroundColor',[0.95,0.95,0.95]);
    handles.edity = uicontrol('Style','edit',...
                              'Position',[121 (D2+76) 51 22], ...
                              'String',num2str(handles.O(2)), ...
                              'BackgroundColor',[1,1,1], ...
                              'CallBack',@editY_Callback);
    uicontrol('style','text','String', 'y','position', ...
              [141 (D2+99) 31 21], ...
              'visible','on','HorizontalAlignment','left', ...
              'BackgroundColor',[0.95,0.95,0.95]);
    handles.editz = uicontrol('Style','edit',...
                              'Position',[221 (D2+76) 51 22], ...
                              'String',num2str(handles.O(3)), ...
                              'BackgroundColor',[1,1,1], ...
                              'CallBack',@editZ_Callback);
    uicontrol('style','text','String', 'z','position', ...
              [241 (D2+99) 31 21], ...
              'visible','on','HorizontalAlignment','left', ...
              'BackgroundColor',[0.95,0.95,0.95]);

    uicontrol('Style','listbox',...
              'Position',[190 (D1+132-80) 92 270],...
              'String',slicingNormals,...
              'BackgroundColor',[1,1,1], ...
              'Callback', @listboxNormal_Callback);
    uicontrol('style','text','String', 'Slicing Normal:','position', ...
              [190 (D1+333) 100 21], ...
              'visible','on','HorizontalAlignment','left', ...
              'BackgroundColor',[0.95,0.95,0.95]);

    handles.cbbd = uicontrol('Style','checkbox','String','Boundary', ...
                             'Value',1,'Position',[20 (D1+142) 115 23], ...
                             'BackgroundColor',[0.95,0.95,0.95]);
    handles.cbsd = uicontrol('Style','checkbox','String','Crosssection', ...
                             'Value',1,'Position',[21 (D1+112) 115 23], ...
                             'BackgroundColor',[0.95,0.95,0.95]);
    handles.cbed = uicontrol('Style','checkbox','String','Edge Color', ...
                             'Value',1,'Position',[21 (D1+82) 115 23], ...
                             'BackgroundColor',[0.95,0.95,0.95]);

    handles.cbfixScale = uicontrol('Style','checkbox','String','FIX', ...
                                   'Value',0,'Position',[125 (D1+184) 50 23], ...
                                   'BackgroundColor',[0.95,0.95,0.95]);

    uicontrol('style','text','String', 'Color Shift:','position', ...
              [21 (D1+333) 100 21], ...
              'visible','on','HorizontalAlignment','left', ...
              'BackgroundColor',[0.95,0.95,0.95]);
    uicontrol('style','text','String', 'Color Spread:','position', ...
              [20 (D1+263) 100 21], ...
              'visible','on','HorizontalAlignment','left', ...
              'BackgroundColor',[0.95,0.95,0.95]);
    handles.sptxt=uicontrol('style','text','String', 'x 1.0','position', ...
                            [130 (D1+263) 42 21], ...
                            'visible','on','HorizontalAlignment','right', ...
                            'BackgroundColor',[0.95,0.95,0.95]);
    handles.shtxt=uicontrol('style','text','String', '0.0','position', ...
                            [130 (D1+334) 42 21], ...
                            'visible','on','HorizontalAlignment','right', ...
                            'BackgroundColor',[0.95,0.95,0.95]);

    [handles.bdata,handles.tdata] = ...
                  ffreadfile('File1',bddatafile, ...
                             'File2',tetdatafile, ...
                             'Delimiter',';','Format','%f %f %f %f');

    g=figure('Name','View','Position', [150 75 800 600]);
    handles.ax=axes();
    hold on;
    handles.g=g;
    handles=draw_handler(handles);

    %if we know the min and max temperatures of the object at startup
    %we initialize the color slider accordingly.
    if handles.haveSlicingdata
        b=handles.Tmax-handles.Tmin;
        handles.Tshift=[-1.5*b 1.5*b 0];
        handles.Tminfrz=handles.Tmin;
        handles.Tmaxfrz=handles.Tmax;
        %no need to call setcolbar_handler() once more as long as the preset
        %is spread=1; shift=0 
        %handles=setcolbar_handler(handles);
    end

    handles.SliderTShiftCtrl= ...
          uicontrol(handles.f,'style','slider','position', ...
                    [21 (D1+302) 150 21], ...
                    'min',handles.Tshift(1),'max',handles.Tshift(2), ...
                    'Value',handles.Tshift(3), ...
                    'callback',@sliderShiftCtrl_Callback);
    handles.SliderTSpreadCtrl= ...
          uicontrol(handles.f,'style','slider','position', ...
                    [20 (D1+233) 150 21], ...
                    'min',handles.Tspread(1),'max',handles.Tspread(2), ...
                    'Value',handles.Tspread(3), ...
                    'callback',@sliderSpreadCtrl_Callback);

    set(handles.f,'color',[0.95,0.95,0.95]);
    set(handles.f,'Visible','on');
    guidata(handles.f, handles);
end

function pushbuttonExit_Callback(hObject, eventdata, handles)
    handles = guidata(hObject);
    B = ishghandle(handles.g);
    if B
        close(handles.g);
    end
    close(handles.f);
end

function pushbuttonCalc_Callback(hObject, eventdata, handles)
    handles=guidata(hObject);
    handles=draw_handler(handles);
    guidata(hObject, handles);
end

function sliderShiftCtrl_Callback(hObject, eventdata, handles)
    handles=guidata(hObject);
    handles.Tshift(3)=get(hObject,'Value');
    handles=setcolbar_handler(handles);
    guidata(hObject, handles);
end

function sliderSpreadCtrl_Callback(hObject, eventdata, handles)
    handles=guidata(hObject);
    handles.Tspread(3)=get(hObject,'Value');
    handles=setcolbar_handler(handles);
    guidata(hObject, handles);
end

function editX_Callback(hObject, eventdata, handles)
    handles=guidata(hObject);
    handles.O(1)=str2double(get(hObject,'String'));
    guidata(hObject, handles);
end

function editY_Callback(hObject, eventdata, handles)
    handles=guidata(hObject);
    handles.O(2)=str2double(get(hObject,'String'));
    guidata(hObject, handles);
end

function editZ_Callback(hObject, eventdata, handles)
    handles=guidata(hObject);
    handles.O(3)=str2double(get(hObject,'String'));
    guidata(hObject, handles);
end

%reset the color slider settings. if slicing data is available
%set the Tshift slider borders to a preset value
function pushbuttonTreset_Callback(hObject, eventdata, handles)
    handles=guidata(hObject);
    handles.Tshift(3)=0;
    handles.Tspread(3)=1;
    set(handles.SliderTShiftCtrl, 'value', handles.Tshift(3));
    set(handles.SliderTSpreadCtrl, 'value', handles.Tspread(3));
    if handles.haveSlicingdata
        b=handles.Tmax-handles.Tmin;
        handles.Tshift=[-1.5*b 1.5*b 0];
        set(handles.SliderTShiftCtrl, 'min', handles.Tshift(1));
        set(handles.SliderTShiftCtrl, 'max', handles.Tshift(2));
        handles.Tminfrz=handles.Tmin;
        handles.Tmaxfrz=handles.Tmax;
    end
    handles=setcolbar_handler(handles);
    guidata(hObject, handles);
end

function listboxNormal_Callback(hObject, eventdata, handles)
    handles=guidata(hObject);
    list = get(hObject, 'String'); 
    idx = get(hObject, 'Value'); 
    name = list{idx};  
    handles.N=str2double(strtrim(strsplit(name,';')));
    handles=draw_handler(handles);
    guidata(hObject, handles);
end

%slider controlled colorbar
%take min and max temperatures of the actual slice
%multiply spread and shift by Tshift
%--> colorbar is changing every time if a new slice is made and is
%showing the range existing in the slice
function handles=setcolbar_handler(handles)
    logVal=10^(handles.Tspread(3)-1);
    if get(handles.cbfixScale, 'Value')
        %fixed colorbar boundary
        Tmid=handles.Tminfrz+0.5*(handles.Tmaxfrz-handles.Tminfrz);
        dT=logVal*(handles.Tmaxfrz-handles.Tminfrz);    
    else
        %colorbar boundary will change with actual slice
        Tmid=handles.Tmin+0.5*(handles.Tmax-handles.Tmin);
        dT=logVal*(handles.Tmax-handles.Tmin);
        %update those in case FIX checkbox is activated, then we have the same
        %settings as before activating the checkbox
        handles.Tminfrz=handles.Tmin;
        handles.Tmaxfrz=handles.Tmax;
    end
    Tu=Tmid-0.5*dT;
    To=Tmid+0.5*dT;
    caxis(handles.ax,[(Tu+handles.Tshift(3)) (To+handles.Tshift(3))]);
    str=sprintf('x %1.1f',logVal);
    set(handles.sptxt, 'String', str);
    str=sprintf('%1.2f',handles.Tshift(3));
    set(handles.shtxt, 'String', str);
end

function handles=draw_handler(handles)
    %find a valid orthogonal
    if ~all([-handles.N(3) 0 handles.N(1)] == 0)
        A=[-handles.N(3) 0 handles.N(1)];
    end
    if ~all([0 -handles.N(3) handles.N(2)] == 0)
        A=[0 -handles.N(3) handles.N(2)];
    end
    B=cross(A,handles.N);
    %find three points on the plane
    S1=handles.O';
    S2=handles.O'+A';
    S3=handles.O'+B';

    set(0,'CurrentFigure',handles.g);
    clf;
    handles.ax=axes();
    hold on;

    [BX,BY,BZ,BC] = slicebd2patch(handles.bdata,S1,S2,S3);
    [SX,SY,SZ,SC] = slicetet2patch(handles.tdata,S1,S2,S3);

    if ~(isempty(BX) || isempty(BY) || isempty(BZ))
        if get(handles.cbbd, 'Value')
            if ~get(handles.cbed, 'Value');
                patch(BX,BY,BZ,BC,'EdgeColor','none');
            else
                patch(BX,BY,BZ,BC,'EdgeColor',[0 0 0],'LineWidth',1);
            end
        end
    end
    if ~(isempty(SX) || isempty(SY) || isempty(SZ))
        if get(handles.cbsd, 'Value')
            if ~get(handles.cbed, 'Value');
                patch(SX,SY,SZ,SC,'EdgeColor','none');
            else
                patch(SX,SY,SZ,SC,'EdgeColor',[0 0 0],'LineWidth',1);
            end
        end
    end
    %you can also put everything at once
    %patch([SX BX],[SY BY],[SZ BZ],[SC BC]);

    axis('image');
    view(3);
    zlabel('z');
    ylabel('y');
    xlabel('x');
    %possibly there is no data to draw
    if ~(isempty(SC) && isempty(BC))
        handles.haveSlicingdata=true;
        handles.Tmin=min(min([SC BC]));
        handles.Tmax=max(max([SC BC]));
        colormap(jet(128));
        handles.cbar_axes=colorbar;
        title(handles.cbar_axes,'dT[K]');
        handles=setcolbar_handler(handles);
    else
        handles.haveSlicingdata=false;
    end
end
