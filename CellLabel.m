classdef CellLabel < handle
    
    properties
        
        T % Time, this should be either empty, or set BEFORE the labels are being added.
        Reg % Registration object 
        
        saveToFile = false; 
        posname
        pth
        cell2exclude % there is also a cell2exclude in the MultiPositionSingleCellResults class, this one is for backing the cells that excluded by sizes. added by Yanfei
        spacedims = 2; 
        tracked=false; 
        useCC = false; 
    end
    
    properties (SetAccess = protected)
        
        Lbl = zeros(0,0,'uint16');
        CC = struct('Connectivity',{},'ImageSize',{},'NumObjects',{},'PixelIdxList',{})
        Msk
        
        Bnd %boundaries
        
        Regions
        sz
        Nt % number of time points for frames
        Nc % max unique number of cells 
        XY
        
        OldLbl %Lbl before track
    end
    
    properties (Transient = true)
       verbose = true; 
    end
    
    properties (Dependent = true)
       dummped 
       filename
    end
    
    methods
        
        function preAllocateMemory(Lbl,sz,nFrames,Types)
            Lbl.Lbl=zeros(sz(1),sz(2),nFrames,'uint16');
            sz3=size(Lbl.Lbl); 
            for i=1:numel(Types)
                Lbl.Msk.(Types{i})=false(sz3); 
            end
        end
        
        function createCellBoundaries(Lbl)
            %%
            if ~Lbl.tracked
                error('Can only create boundaries for tracked cells, this might change in the future')
            end
            if isempty(Lbl.Lbl) && Lbl.saveToFile
                Lbl.verbose && fprintf('loading CellLabel from drive\n'); %#ok<VUNUS>
                s=load(fullfile(Lbl.pth,sprintf('CellLabel_%s.mat',Lbl.posname)));
                if Lbl.useCC
                    Lbl.CC=s.L;
                else
                    Lbl.Lbl=s.L;
                end
                Lbl.Msk=s.M;
            end
            for j=1:numel(Lbl.Regions)
                %% 
                
                lbl=Lbl.Lbl; 
                lbl(~Lbl.Msk.(Lbl.Regions{j}))=0; 
                prm = bwperim(lbl,4);
                lbl(~prm)=0; 
                for i=1:Lbl.Nt
                    [y,x,id]=find(lbl(:,:,i)); 
                    [bnd,unqid]=grp2cell([x y],id); 
                    Lbl.Bnd.(Lbl.Regions{j}){i}=cell(Lbl.Nc,1); 
                    Lbl.Bnd.(Lbl.Regions{j}){i}(unqid)=bnd; 
                end
            end
        end
        
        function filename = get.filename(Lbl)
            filename = fullfile(Lbl.pth,sprintf('CellLabel_%s.mat',Lbl.posname)); 
        end
        
        function dummped = get.dummped(Lbl)
            dummped = Lbl.saveToFile && isempty(Lbl.Lbl) && isempty(Lbl.Msk) && exist(Lbl.filename,'file');
        end
        
        function dumpToDrive(Lbl,pth,posname)
            if nargin==1 && ~Lbl.saveToFile || Lbl.dummped
                return
            end
            if nargin>1
                Lbl.pth=pth; 
            end
            if nargin>2
                Lbl.posname=posname; 
            end
            
            % method will save current Lbl and Msk variables into a CellLabel.mat file
            assert(~isempty(Lbl.pth),'Please specify pth & posname properties before calling dumpToDrive');
            assert(~isempty(Lbl.posname),'Please specify pth & posname properties before calling dumpToDrive');

            if Lbl.useCC
                L=Lbl.CC; 
            else
                L=Lbl.Lbl;  %#ok<*NASGU>
            end
            M=Lbl.Msk; 
            B=Lbl.Bnd; 
            OL=Lbl.OldLbl; 
            save(Lbl.filename,'L','M','OL','B', '-v7.3'); 
            
            Lbl.Lbl=[]; 
            Lbl.Msk=[]; 
            Lbl.OldLbl=[];
            Lbl.Bnd=[]; 
            Lbl.saveToFile=true; 
        end
        
        
        
        function tf = isempty(Lbl)
            tf = isempty(Lbl.Lbl); 
        end
        
        function sz = get.sz(Lbl)
            sz = [size(Lbl.Lbl,1) size(Lbl.Lbl,2)]; 
        end
        
        function Nt = get.Nt(Lbl)
            Nt = numel(Lbl.T); 
        end
        
        function Nc = get.Nc(Lbl) 
            Nc = Lbl.Nc;
%             Nc=max(max(max(Lbl.Lbl))); 
        end
        

        function Regions = getRegions(Lbl) %create byYanfei
            Regions = Lbl.Regions;
        end
        
        function Bnd = getCellBoundaries(Lbl,cellid,T)
            Bnd = cell(numel(T),1);
            for i=1:numel(T)
                lbl = Lbl.getLbls('base',T(i)); 
                bnd = bwboundaries(lbl==cellid); 
                Bnd{i} = bnd{1}(:,[2 1]);  
            end
                       
        end
        
        
        function R = getRegionMatrix(Lbl,T)
            
            % get index based on timepoint T (1 if empty)
            if nargin==1 || isempty(T)
                ix = 1; 
            else
                [~,ix] = min(abs(Lbl.T-T)); 
            end
            
            % create the empty matrix
            R = zeros(Lbl.sz); 
            for i=1:numel(Lbl.Regions)
                R(Lbl.Msk(:,:,ix).(Lbl.Regions{i}))=i; 
            end
            
        end
              
        function xy= getXY(Lbl,T,varargin)
            arg.type = 'base';
            arg = parseVarargin(varargin,arg); 
            
            if nargin==1 
                T=[]; 
            end
            lbl = getLbls(Lbl,arg.type,T);
            [Xgrd,Ygrd]=meshgrid(1:size(lbl,2),1:size(lbl,1)); 
            [~,~,x]=grp2cell(Xgrd(:),lbl(:)); 
            [~,id,y]=grp2cell(Ygrd(:),lbl(:));
            xy=nan(Lbl.Nc,2); 
            xy(id,1)=x; 
            xy(id,2)=y; 
%             cntr = regionprops(lbl,'Centroid');
%             xy = cat(1,cntr.Centroid); 
        end
        
        function addLbl(Lbl,newlbl,type,T,varargin)
            
            % to save memory & harddrive space move to uint16
            newlbl = uint16(newlbl); 
            
            % defults
            arg.relabel = 'none'; 
            arg.maxcelldistance = 25; % distance for tracking, 25 is default
            arg = parseVarargin(varargin,arg); 
            
            
            %% allow a "data dump" of a whole matrix of Labels in the case where type=base and relabel=none and we got accurate number of frames. 
            if strcmp(arg.relabel,'none') && strcmp(type,'base') && isempty(Lbl.T) && numel(T)==size(newlbl,3)
                Lbl.T=T; 
                Lbl.Lbl=newlbl; 
                Lbl.Msk.base = newlbl>0; 
                return; 
            end
            
            %% otherwise....
              
            
            % default value for type is base
            if ~exist('type','var') || nargin==2 || isempty(type)
                type='base'; 
            end
                        
            % if T is not specified change the first (and only!) layer 
            % if T is specified, either get the index for that timepoint or
            % add another timepoint in the end and get its index. 
            if nargin==3 || isempty(T)
                assert(isempty(Lbl.T),'Can only add a base label without a timepoint when Label doesn''t have a Time vector!'); 
                ix = 1; 
                T=[]; 
            else % find the existing exact timepoint, or add a new one in the end
                ix = find(ismember(Lbl.T,T));
                assignin('base', 'addLblT', ix); %assign addLblT inthe  base workspace
                if isempty(ix)
                    Lbl.T(end+1)=T; 
                    ix = numel(Lbl.T); 
                end
            end
               
            
            if strcmp(type,'base') 
                if  ~isempty(Lbl.Lbl)
                    switch arg.relabel
                        case 'none'
                        case 'nearest'
                            Lbl.tracked=1; 
                            oldlbl = newlbl; % call the newlbl oldlbl so I can then rename it into newlbl and then keep the same variable name
                            prps = regionprops(oldlbl,{'Centroid','PixelIdxList'});
                            xy = cat(1,prps.Centroid);
                            PxlIdx = {prps.PixelIdxList};
                            
                            % relabel the base according to distances to the closest
                            possT = Lbl.T(setdiff(1:numel(Lbl.T),ix));
                            [~,ix2] = min(abs(possT-T));
                            XY = getXY(Lbl,Lbl.T(ix2));
                            
                            M=annquerysingle(xy',XY', arg.maxcelldistance); 
                            
                            newlbl=zeros(Lbl.sz);
                            for i=1:size(M,1)
                                newlbl(PxlIdx{M(i,2)})=M(i,1);
                            end
                            
                        case 'lapjv' % pretty slow, but optimal assignment.
                            Lbl.tracked=1; 
                            oldlbl = newlbl; % call the newlbl oldlbl so I can then rename it into newlbl and then keep the same variable name
                            prps = regionprops(oldlbl,{'Centroid','PixelIdxList'});
                            xy = cat(1,prps.Centroid);
                            PxlIdx = {prps.PixelIdxList};
                            
                            % relabel the base according to distances to the closest
                            possT = Lbl.T(setdiff(1:Lbl.Nt,ix));
                            [~,ix2] = min(abs(possT-T));
                            XY = getXY(Lbl,Lbl.T(ix2));
                            
                            D=distance(XY',xy');
                            
                            n1=size(D,1); % label already in Lbl
                            n2=size(D,2); % new label
                            ass = lapjv(D,0.1);
                            if n1>n2
                                id = ass;
                            else
                                id = zeros(n2,1);
                                id(ass)=1:n1;
                                id(id==0)=(n1+1):n2;
                            end
                            newlbl=zeros(Lbl.sz);
                            for i=1:n2
                                newlbl(PxlIdx{i})=id(i);
                            end
                    end
                end
                if ~isa(newlbl,'uint16')
                    newlbl = uint16(newlbl); 
                end
                if Lbl.useCC
                    PxlIds = regionprops(newlbl,'PixelIdxList');
                    Lbl.CC(ix).PixelIdxList=PxlIds;
                    Lbl.CC(ix).NumObjects=numel(PxlIds);
                    Lbl.CC(ix).Connectivity
                    Lbl.CC(ix).ImageSize=size(newlbl);
                else
                    try Lbl.Lbl(:,:,ix)=newlbl;     
                    catch
                        error('This Lbl has two timestamps and only one image')
                    end
                    
                end
            end
            % set up the type in the Msk datastructure. 
            Lbl.Msk.(type)(:,:,ix)=newlbl>0;  %this line is mainly what it does.
            
            
            if isempty(Lbl.Nc)
                if Lbl.useCC
                    Lbl.Nc = Lbl.CC(ix).NumObjects;
                else
                    Lbl.Nc =max(newlbl(:)); 
                end
            else
                if Lbl.useCC
                    Lbl.Nc = [max(Lbl.Nc) Lbl.CC(ix).NumObjects];
                else
                    Lbl.Nc =max([Lbl.Nc max(newlbl(:))]); 
                end
            end
            Lbl.Regions = union(Lbl.Regions,{type}); 
            
        end
        
        function Traj = getTrajectories(Lbl)
            Traj=cellfun(@(~) nan(Lbl.Nt,2),cell(Lbl.Nc,1),'uniformoutput',0);
            for i=1:numel(Lbl.T)
                xy=Lbl.getXY(Lbl.T(i));
                for j=1:size(xy,1)
                    Traj{j}(i,:)=xy(j,:); 
                end
                
            end
        end
                
        function Traj = track(Lbl,varargin)
            arg.maxdisp=25; 
            arg.mem=3; 
            arg.minlength=0.75; 
            arg.reset=false; % can do it many times
            arg.jitter='none'; 
            arg.pairrule='fwdnn';
            arg.subsetofcells = {}; 
            arg = parseVarargin(varargin,arg); 
            
            if arg.reset
                if isempty(Lbl.OldLbl)
                    Lbl.verbose && fprintf('loading CellLabel from drive\n');  %#ok<VUNUS>
                    s=load(fullfile(Lbl.pth,sprintf('CellLabel_%s.mat',Lbl.posname))); 
                    Lbl.OldLbl=s.OL;
                end
                Lbl.Lbl=Lbl.OldLbl; 
                return
            end
            
            %% first create the cell array with in each frame the XY and ID for each cell
            Pnts=cell(numel(Lbl.T),1);
            for i=1:numel(Lbl.T)
                Pnts{i}=Lbl.getXY(Lbl.T(i));%looks like all regions are used
                ix=find(~isnan(Pnts{i}(:,1)));
                Pnts{i}=[Pnts{i}(ix,:) ix(:)]; %last col is the cell ID
            end
            %% now track the cells
            Traj=ultTrack(Pnts,'dim',1:2,'maxdisp',arg.maxdisp,'mem',arg.mem,'minlength',arg.minlength);
            
            %% relabel Lbl into cell id that are based on indexes of Traj
            % first create a conversion matrix that says which new ids from
            % lbl belong to which new cell id (Traj index) for each
            % timepoint.
            Conv=nan(numel(Lbl.T),numel(Traj));
            for i=1:numel(Traj),
                Conv(Traj{i}(:,4),i)=Traj{i}(:,3);
            end
            
            %% 
            Lbl.OldLbl=Lbl.Lbl; %save the OldLbl
            
            
            %% use conversion matrix and the pixel idx (regionprops)
            % to assign the pixel idx in newlbl that belong to the old id
            % the new id that now exist in Conv matrix
            L=Lbl.Lbl; 
            
            parfor i=1:Lbl.Nt
                oldlbl=L(:,:,i);
                %%
                Pix=regionprops(oldlbl,'PixelIdxList');
                ix=find(~isnan(Conv(i,:)));
                newlbl=zeros(size(oldlbl),'uint16');
                for j=1:numel(ix)
                    % the pixel idx should reference the old labeling which
                    % is the actual values in Conv. It should get the new
                    % label which are the collum idx in Conv. Do this for
                    % each row
                    newlbl(Pix(Conv(i,ix(j))).PixelIdxList)=ix(j); %#ok<PFBNS>
                end
                L(:,:,i)=newlbl;
            end
            Lbl.Lbl=L; 
            
            %% marked as labeled
            Lbl.tracked=true;
            Lbl.Nc=numel(Traj); 
        end

        function lbl = get.Lbl(Lbl)
            if Lbl.useCC
                lbl = zeros([Lbl.CC(1).ImageSize numel(Lbl.CC)],'uint16');
                for i=1:numel(Lbl.CC)
                    lbl(:,:,i)=labelmatrix(Lbl.CC(i));
                end
            else
                lbl = Lbl.Lbl;
            end
        end
        
        
        function imshow(Lbl,type,T)
            if nargin<2
                type=[]; 
                T=[]; 
            elseif nargin<3
                T=[]; 
            end
            lbl = Lbl.getLbls(type,T); 
            imshow(label2rgb(lbl,'jet','k','shuffle'));
        end
        
       
        function Movie = showVideoWithOutlines(Lbl,stk,lvl,response,varargin)
            %% 
            % input parameters: 
            % Lbl - The calling CellLabel object
            % stk - a stack of images usually from stkread(MD,...) possibly
            %       after registration (with Lbl.Reg.register(stk,Lbl.T) etc. 
            % lvl - [min max] of the contrast in 2^16 bit space. 
            % response- a matrix of responses per cell. Or a cell array of
            %           up to three of such matirces to show up multiple
            %           plots. Must have the same # of cols as cells. This
            %           is what measureStackIntensities returns. 
            % 
            % Other options: called as property/value pairs: 
            % cell_to_track : chose the id of the cell to zoom into
            % fig : use a specific figure, if empty will create one. 
            % ylim: a matrix of the ylim for the response plots, a row per
            %       plot (if resposne is a cell array)
            % position: the poisition of the figure on the screen
            % 
            
            arg.cell_to_track=[]; 
            arg.zoomed=false; 
            arg.fig=[]; 
            arg.movie=false;
            arg.region='base'; 
            arg.ylim=[0 1.5]; 
            arg.position=[0.05 0.05 0.85 0.85]; 
            arg.t=Lbl.T; 
            arg.width=100; 
            arg = parseVarargin(varargin,arg);
                   
            
            %% get all XY for all frames
            if isempty(Lbl.XY)
                XY=cell(Lbl.Nt,1);  %#ok<*PROP>
                h=waitbar(0,'calculating centroids');
                for i=1:Lbl.Nt
                    XY{i}=Lbl.getXY(Lbl.T(i));
                    waitbar(i/Lbl.Nt);
                end
                delete(h);
                Lbl.XY=XY; 
            else
                XY=Lbl.XY; 
            end
            if isempty(Lbl.Bnd) && Lbl.saveToFile
                Lbl.verbose && fprintf('loading CellLabel boundaries from drive\n'); %#ok<VUNUS>
                s=load(fullfile(Lbl.pth,sprintf('CellLabel_%s.mat',Lbl.posname)),'B');
                Lbl.Bnd=s.B; 
            end
                
            if isempty(arg.fig)
                arg.fig=figure; 
            end
            
            XY=cellfun(@ceil,XY,'uniformoutput',0); 
            
            if nargin<4 
                response = {}; 
            elseif ~iscell(response)
                response = {response};
            end
       
            
            %%
            if max(lvl)>1
                lvl=lvl/2^16;
            end
            ttl='video';
            if ~isempty(inputname(2))
                ttl=inputname(2);
            end
            if isempty(Lbl.Bnd)
                Lbl.createCellBoundaries; 
            end
            image_axis=[]; 
                      
            %% define colors
            clr=jet(double(Lbl.Nc)); 
            clr=clr(randperm(size(clr,1)),:); 
            
            %%
            cell_to_track=arg.cell_to_track;
            if ~isempty(cell_to_track)
                arg.zoomed=true;
            end
            
            %% start video player
            play_fps=3; 
            [fig_handle, image_axis,~,play_func,plot_handles] = videofig(size(stk,3),@redraw,play_fps,[],[],numel(response),'position',arg.position,'fig',arg.fig,'Name',ttl,'numbertitle','off');
            
            % setup initial zoom factor
            zoomed = arg.zoomed; 
            set(fig_handle,'userdata',zoomed);
            
            % determine if the output should be a avi file or a 
            if ~arg.movie
                Movie=[]; 
                play_func(1./play_fps);
            else
                Movie=struct('cdata',{},'colormap',{}); 
                for i=1:size(stk,3)
                    redraw(i); 
                    Movie(i) = getframe(fig_handle); 
                end
            end
            
            %% nested drawing function 
            function redraw(frame)
                %%
                h = findobj(fig_handle,'tag','plot_axes');
                zoomed = get(fig_handle,'userdata'); 
                
                % find the closets frame in Lbl by Time
                [~,closestFrameInLbl] = min((Lbl.T-arg.t(frame)).^2);
                
                bnd=Lbl.Bnd.(arg.region){closestFrameInLbl}; 
                bnd(cellfun(@isempty,bnd))={[nan nan]};
                
                if zoomed
                    set(fig_handle,'Units','norm')
                    if isempty(cell_to_track)
                        pnt = get(fig_handle, 'CurrentPoint');
                        pnt=[pnt(1,1) 1-pnt(1,2)];
                        if numel(response)==0
                            pnt(1)=ceil(pnt(1)*size(stk,2));
                        else
                            pnt(1)=ceil(pnt(1)*size(stk,2)*2);
                        end
                        pnt(2)=ceil(pnt(2)*size(stk,1));
                        d=distance(pnt(:),XY{closestFrameInLbl}'); 
                        [~,cell_to_track]=min(d);
                    end
                    % update the "click point" to the the center of the
                    % closets cell. If XY do not exist for this cell at
                    % this time, find the XY from the next non NAN frame
                    % set the "clicked" pnt as the center of the cell
                    pnt = XY{closestFrameInLbl}(cell_to_track,:);
                    f=closestFrameInLbl;
                    while isnan(pnt(1))
                        if f>=numel(XY)
                            f=0;
                        end
                        f=f+1; 
                        pnt = XY{f}(cell_to_track,:);
                    end
                                    
                    set(fig_handle,'CurrentAxes',image_axis);
                    hold off
                    xlm=max(1,pnt(1)-arg.width):min(size(stk,2),pnt(1)+arg.width);
                    ylm=max(1,pnt(2)-arg.width):min(size(stk,1),pnt(2)+arg.width);
                    imshow(stk(ylm,xlm,frame),lvl,'XData',xlm,'YData',ylm)
                    hold on
                    if ~isempty(bnd{cell_to_track})
                        plot(bnd{cell_to_track}(:,1),bnd{cell_to_track}(:,2),'.','color',clr(cell_to_track,:));
                    end
                    for ii=1:numel(response)
                        set(fig_handle,'CurrentAxes',plot_handles(ii));
                        hold off
                        plot(plot_handles(ii),response{ii}(:,cell_to_track))
                        ylim(arg.ylim(ii,:))
                        hold on
                        plot(plot_handles(ii),[frame frame],arg.ylim(ii,:))
                    end
                    
                else % show all cells
                    cell_to_track=[];
                    set(h,'Visible','off');
                    set(fig_handle,'CurrentAxes',image_axis);
                    hold off
                    imshow(stk(1:2:end,1:2:end,frame),lvl)
                    hold all
                    set(image_axis,'ColorOrder',clr)
                    id = cellfun(@(m,cnt) cnt*ones(size(m,1),1),bnd,num2cell((1:numel(bnd))'),'uniformoutput',0);
                    id=cat(1,id{:});
                    bnd=cat(1,bnd{:});
                    plot(bnd(:,1)/2,bnd(:,2)/2,'c.','markersize',1)
                    
                    %                     for ii=1:numel(Lbl.Bnd.base{frame})
                    %                         plot(bnd{ii}(:,1)/2,bnd{ii}(:2)/2,'.');
                    %                     end
                    %
                    for ii=1:numel(response)
                        set(fig_handle,'CurrentAxes',plot_handles(ii));
                        hold off
                        plot(nanmean(response{ii},2),'-k')
                        ylim(arg.ylim(ii,:))
                        hold on
                        plot([frame frame],arg.ylim(ii,:))
                    end
                end
            end
        end
        
        function exploreCellsVideosWithScatterPlot(Lbl,Stk,X,Y,M,varargin)
            
            arg.scatter_fig=1;
            arg.video_fig=2;
            arg.region='ring';
            arg.lvl=[100 400];
            arg.ylim=[0 1]; 
            arg = parseVarargin(varargin,arg); 
            scatter_fig=arg.scatter_fig; 
            video_fig=arg.video_fig;
            region=arg.region; 
            lvl=arg.lvl; 
            
            
            figure(scatter_fig)
            clf
            set(scatter_fig,'WindowButtonDownFcn',@button_down,'position', [48   314   610   554]);
            plot(X,Y,'.');
            hold on
            unplot set
            scatter_axes=gca;
            
            
            function button_down(~,~)
                pnt = get(scatter_axes, 'CurrentPoint');
                pnt=pnt(1,1:2);
                d=distance(pnt(:),[X(:) Y(:)]');
                [~,cell_to_track]=min(d);
                unplot revert
                plot(X(cell_to_track),Y(cell_to_track),'ro');
                if ismember(video_fig,get(0,'children'))
                    delete(video_fig);
                end
                showVideoWithOutlines(Lbl,Stk,lvl,M,'cell_to_track',cell_to_track,'fig',video_fig,'position',[0.3589    0.2285    0.6042    0.6327],'region',region,'ylim',arg.ylim)
            end
            
            
            
        end
        
        function time = get.T(Lbl)
            time = Lbl.T;
        end
        
        
        
        function lbl = getLbls(Lbl,type,T)
            if nargin==1 || isempty(type)
                type='base';
                T=[];
            end
            
            % decide which slice to get. 
            if nargin==2 || isempty(T)
                ix = 1; 
            else
                ix = nan(size(numel(T))); 
                for i=1:numel(T)
                    [~,ix(i)] = min(abs(Lbl.T-T(i)));
                end
            end
            
            if ~Lbl.useCC
                if isempty(Lbl.Lbl) && Lbl.saveToFile
                    Lbl.verbose && fprintf('loading CellLabel from drive\n'); %#ok<VUNUS>
                    s=load(fullfile(Lbl.pth,sprintf('CellLabel_%s.mat',Lbl.posname)));
                    Lbl.Lbl=s.L;
                    Lbl.Msk=s.M;
                end
            else
                if isempty(Lbl.CC) && Lbl.saveToFile
                    Lbl.verbose && fprintf('loading CellLabel from drive\n'); %#ok<VUNUS>
                    s=load(fullfile(Lbl.pth,sprintf('CellLabel_%s.mat',Lbl.posname)));
                    Lbl.CC=s.L;
                    Lbl.Msk=s.M;
                end
            end
            if Lbl.useCC
                lbl = labelmatrix(Lbl.CC(ix)); 
            else
                lbl = Lbl.Lbl(:,:,ix);
            end
            
            % zero out everything that is not in the required region type 
            lbl(~Lbl.Msk.(type)(:,:,ix))=0; 
            
            % if there are cells that don't have cyto or rings, then using
            % this function gives different number of cyto or rings.
        end
        
        function lbl = getOldLbls(Lbl,type,T)
            %added by Yanfei, similar to getLbls method, just get the old
            %labels before tracking, first used to debugging purpose.
            if nargin==1 || isempty(type)
                type='base';
                T=[];
            end
            
            % decide which slice to get. 
            if nargin==2 || isempty(T)
                ix = 1; 
            else
                ix = nan(size(numel(T))); 
                for i=1:numel(T)
                    [~,ix(i)] = min(abs(Lbl.T-T(i)));
                end
            end
            
            if isempty(Lbl.OldLbl) && Lbl.saveToFile
                Lbl.verbose && fprintf('loading CellLabel from drive\n'); %#ok<VUNUS>
                s=load(fullfile(Lbl.pth,sprintf('CellLabel_%s.mat',Lbl.posname))); 
                Lbl.OldLbl=s.OL; 
                Lbl.Msk=s.M; 
            end
            lbl = Lbl.OldLbl(:,:,ix); 
            
            % zero out everything that is not in the required region type 
            lbl(~Lbl.Msk.(type)(:,:,ix))=0;             
        end
        
        function lbl = loadLbls(Lbl, posname)
            % added by Yanfei, similar to getLbls method, just load the
            % whole labels, use when labels have been dump to drive
            Lbl.verbose && fprintf('loading CellLabel from drive\n'); %#ok<VUNUS>
            s=load(fullfile(Lbl.pth,sprintf('CellLabel_%s.mat',posname)));
            
            lbl.Lbl=s.L;
            lbl.Msk=s.M;
            lbl.OldLbl=s.OL;
            lbl.Bnd=s.B;
        end
        
        function undoTracking(Lbl)
            % added by Yanfei, change Lbl to OldLbl, used for debugging
            % purpose. It becomes useless when I find the track method has a reset
            % option. 
            warning('This function only reset Lbl to OldLbl and set OldLbl to empty, beware of other changes you have made');
            if Lbl.tracked                
                if ~isempty(Lbl.OldLbl)
                    Lbl.Lbl = Lbl.OldLbl;
                    Lbl.OldLbl = [];
                    Lbl.tracked = false;
                else
                    fprintf('OldLbl is empty, use loadLbls to load full label\n');
                end
            else
                fprintf('OldLbl is empty, tracking hasnot been performed\n');
            end
        end
        
        
        function M = meanIntensityPerLabel(Lbl,Stk,T,varargin)
            % method will calculate a function (mean by defualt) per id per
            % time. It will use the label matrix that is closest in time to
            % the image measured in. 
            arg.type = 'base'; 
            arg.func = 'mean'; 
            arg.interp = []; 
            arg = parseVarargin(varargin,arg); 
                        
            if isempty(Lbl.T) || nargin==2 || isempty(T)
                lbl = getLbls(Lbl,arg.type);
                M = meanIntensityOverTime(Stk,lbl,arg.func);
            else
                assert(numel(T) == size(Stk,3),'Must provide a timepoint for each slice in Stack');
                if numel(T)>2
                    assert(min(diff(T))>0,'Timeseries must be monotonically increasing');  
                end
                M=nan(numel(T),Lbl.Nc); 
                dt = abs(repmat(Lbl.T(:),1,numel(T))-repmat(T(:)',numel(Lbl.T(:)),1));
                [~,ix]=min(dt,[],1);
                unq = unique(ix);

                % Made changes to meanIntensityOverTime to result in missing labels
                % becoming a column of NaN.
                for i=1:numel(unq)
                    lbl = getLbls(Lbl,arg.type,Lbl.T(unq(i)));
                    s = Stk(:,:,ix==unq(i));
                    % there could be a case where max(lbl(:))<Lbl.Nc in
                    % those cases we need to match M till max of lbl
                    % in most cases size(m,2)==Lbl.Nc. In any case the col
                    % number for a cell should match its label number and
                    % the XY position it has etc. 
                    m = meanIntensityOverTime(s,lbl,arg.func); 
                    M(ix==unq(i),1:size(m,2)) = m; 
                end
            end
            if ~isempty(arg.interp)
                M=interp1(T,M,arg.Tinterp); 
            end
            
        end

        
        function stkshowLabel(Lbl,varargin)
            arg.shuffle = true;
            
            
        end
        

        
        function allRegionProps = cellRegionProps(Lbl, varargin)
            % this is an substitue of the matlab regionpros function,            
            for j = 1:numel(Lbl.Regions)
                region = Lbl.Regions{j};
                for i = 1:numel(Lbl.T)
                    % get Labels first
                    Time = Lbl.T(i);
                    lbl = Lbl.getLbls(region,Time);
                    props{i} = regionprops(lbl, varargin);                 
                end
                %allRegionProps(i,:) = props;
                allRegionProps.region = props;
            end
        end                

        
        
%         %% the next function is added by Yanfei
%         function excludeSize(Lbl, varargin)
%             % Summary of this function goes here
%             
%             % varargin in format like Type1, minsize1, maxsize1, Type2, minsize2,
%             % maxsize2,.....
%             
%             % this function does not delete any cells, it only gives an index of cells
%             % that meet your exclusion criteria
%             
%             %% check whether varargin has right number and correct datatype
%             if mod(nargin-1,3) == 0
%             else
%                 error('please provide  varargin in format like Type1, minsize1, maxsize1, Type2, minsize2, maxsize2,.....');
%             end
%             
%             Types = varargin(1:3:end);
%             Mins = varargin(2:3:end);
%             Maxes = varargin(3:3:end);
%             
%             if any(~isstr(Types))
%             else
%                 error('Types should string datatype');
%             end
%             
%             if any(~isnumeric(Mins))
%             else
%                 error('min should be numeric datatype');
%             end
%             
%             if any(~isnumeric(Maxes))
%             else
%                 error('max should be numeric datatype');
%             end
%             %% check whether types (like 'nuc', 'cyt') are all present
%             if all(ismember(Types, Lbl.Regions))
%             else
%                 error('Types not found in the CellLabel');
%             end
%             
%             %%
%             
%             %% get labels
%             for i = numel(Types)
%                 regionLbl{i} = getLbls(Lbl, Types{i}, Lbl.T);
%                 % getLbls methods has the danger of get cells that don't
%                 % have rings or cytos, 
%             end
%             
%             %% get the area property of each cells 
%             %replace whole loop using allRegionProps method
%             for i = 1:numel(Lbl.T)
%                 for j = 1:numel(Types)
%                     regionArea{j}{i} = regionprops(regionLbl{j}(:,:,i), 'Area');
%                 end
%             end
%             
%             %% get the indexes that meet the size requirements
%             cellnum = cellfun(@numel,regionArea{1}); % has cell numbers of each frame
%  
% 
%             cell2exclude = [];
%             
%             for i = 1:numel(Lbl.T)
%                 allcells = 1:cellnum(i);
%                 for j = 1:numel(Types)
%                     cell2exclude_CurrentType = allcells([regionArea{j}{i}.Area]>cell2mat(varargin(3*j)) | [regionArea{j}{i}.Area]<cell2mat(varargin(3*j-1)));
%                     cell2exclude = union(cell2exclude,cell2exclude_CurrentType); %these number need to be change if you change ring_width etc
%                 end
%             end
%             Lbl.cell2exclude=union(cell2exclude,Lbl.cell2exclude); %Lbl.cell2exclude can be changed by other methods
%         end


    end
    
   

    
    
end
        