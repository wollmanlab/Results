function GUIlookAndAnnotateDatasets

%% make comaptible with multiple microscopes
if isunix
    basepth = '/data4/Images'; 
elseif ispc
    basepth = 'E:\WollmanLabEpiScopeData'; 
end

%% define GUI layout
figure(1)
clf
set(1,'position',[106   150   560   500]); 

Channels = {''}; 
Types = {''};

hDSbotton = uicontrol('String','Choose Dataset','callback',@chooseDS,'BackgroundColor','red','position',[20 450 100 25]);
hDSname = uicontrol('style','text','string','','position',[125 450 300 25]); 

hPosition = uicontrol('String','Choose Position','callback',@choosePos,'BackgroundColor','red','position',[200 400 100 25]);
PosChosen = false; 

uicontrol('style','text','string','Channel','position',[430 405 100 20]);
hChannelmenu = uicontrol('style','popupmenu','string',Channels,'position',[325 400 100 25],'BackgroundColor','red','callback',@chooseChannel); 
ChannelChosen = false; 


uicontrol('style','text','string','Types','position',[430 355 100 20]);
hShowPorperty = uicontrol('style','popupmenu','string',Types,'position',[325 350 100 25],'callback',@showheatmap);

hWatch = uicontrol('String','Watch','BackgroundColor','red','position',[20 380 100 25],'callback',@watch);

hSave = uicontrol('String','Save','position',[500 20 50 20],'BackgroundColor','red','callback',@savefcn);
hComment = uicontrol('Style','Edit','Max',5,'position',[20 50 520 300],'HorizontalAlignment','left');

hExportMD = uicontrol('String','Export MD','position',[20 20 70 20],'BackgroundColor','m','callback',@exportMD);  %#ok<NASGU>

MD = []; 
Pos = {}; 
Scp.Chamber = [];  
pth=''; 
Rslt = [];  %#ok<NASGU>
Plt=[]; 

    function exportMD(~,~,~)
        assignin('base','MD',MD); 
    end

    function chooseDS(~,~,~)
        pth = uigetdir(basepth); 
        MD = Metadata(pth); 
        
        if ismember('PlateType',MD.Types)
            ptype = unique(MD,'PlateType'); 
            Plt = Plate(ptype{1});
        else
            Plt = Plate;
        end
        Plt.x0y0=[0 0];
        Scp.Chamber=Plt; 
%         Rslt = Results(pth);
        Scp.Chamber=Plt; 
        Plt.x0y0=[0 0];
%         set(hComment,'String',Rslt.Conclusions); 
        Channels = unique(MD,'Channel'); 
        set(hChannelmenu,'string',Channels); 
        Types = setdiff(MD.Types,{'Channel'
            'Exposure'
            'PixelSize'
            'Position'
            'TimestampFrame'
            'TimestampImage'
            'XY'
            'frame'
            'Fluorophore'
            'Marker'
            'Skip'
            'dZ'
            'group'});
        set(hShowPorperty,'string',Types); 
        
        set(hDSbotton,'BackgroundColor','green'); 
        [~,ds]=fileparts(pth); 
        set(hDSname,'String',ds);
    end

    function choosePos(~,~,~)
        Wells = unique(MD,'group');
        msk =  double(ismember(Plt.Wells,Wells))*0.5; 
        Pos = markPlatePosition(Scp,'fig',2,'msk',msk);
        Pos = unique(Pos); 
        set(hPosition,'BackgroundColor','green');
        PosChosen = true;
        if ChannelChosen
            set(hWatch,'BackgroundColor','green');
        end
    end

    function chooseChannel(~,~,~)
        set(hChannelmenu,'BackgroundColor','green'); 
        ChannelChosen = true; 
        if PosChosen
            set(hWatch,'BackgroundColor','green'); 
        end
    end

    function watch(~,~,~)
        if PosChosen && ChannelChosen
            
            %% show all chosen stacks
            stk = stkread(MD,'group',Pos,'Channel',Channels{get(hChannelmenu,'value')},...
                             'resize',0.33,'groupby','group','sortby','TimestampImage','montage',true); % all the different flags
            ttl = cellfun(@(p) sprintf('%s-%s',p,Channels{get(hChannelmenu,'value')}),Pos,'uniformoutput',0);             
            stkshow(stk,'title',ttl); 
       
        else
            msgbox('Choose position and channel before pressing watch'); 
        end
    end

    function showheatmap(~,~,~)
        if isempty(MD)
            return
        end
        try 
            Type = Types{get(hShowPorperty,'value')}; 
            figure(999)
            clf
            MD.plotMetadataHeatmap(Type,'removeempty',1); 
        catch e
            close(999)
            warning('Error during show heatmap: %s',e.message); 
        end
        
    end

    function savefcn(~,~,~)
        Conclusions = get(hComment,'String'); 
        save(fullfile(pth,'Conclusions.mat'),'Conclusions', '-v7.3')
        if ~exist(fullfile(pth,'Results.mat'),'file')
            R = Results; 
            R.pth=pth; 
            R.Conclusions = Conclusions; 
            R.saveResults; 
        end
        set(hSave,'backgroundcolor','green'); 
    end

end
