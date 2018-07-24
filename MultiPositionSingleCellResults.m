classdef MultiPositionSingleCellResults < MultiPositionResults
    
    properties % only use the specialized setLbl / getLbl 
        Lbl = CellLabel.empty(1,0); % an empty CellLabel array. It will have an element for each potision.
        TimeVecs
        
        SubPlotData=struct('nrows',{},'ncols',{},'title',{}); 
        cell2exclude  %added by Yanfei
        dataBeforeExclusion %added by Yanfei
        oldDataStored = false; %added by Yanfei
    end
    
    
    methods (Static)
        function R = loadobj(S)
            R = MultiPositionSingleCellResults; 
            R = reload(R,S); 
        end
    end
    
    methods
        
        function S = saveobj(R)
            S = toStruct(R); 
        end
        
        function R = reload(R,S)
            R = reload@MultiPositionResults(R,S); 
            if isfield(S,'Lbl') % is the fields exists load them, if not, 
                                % they will be filled with default values
                                % effectivly upcasting an object. 
                R.Lbl = S.Lbl;
                R.TimeVecs = S.TimeVecs;
            end
            if isfield(S,'cell2exclude') %added by Yanfei
                R.cell2exclude = S.cell2exclude;
            end
        end
        
        function plotByPosition(R,H,func_plot,varargin)
            
            assert(~isempty(R.SubPlotData),'please set SubPlotData property before you call plotByPosition')
            
            arg.fig=[]; 
            arg.hold='off'; 
            arg.ylim=[]; 
            arg.xlim=[]; 
            arg.xlabel=R.SubPlotData.xlabel; 
            arg.ylabel=R.SubPlotData.ylabel;
            arg.title=R.SubPlotData.title; 
            arg.ncols=R.SubPlotData.ncols; 
            arg.nrows=R.SubPlotData.nrows; 
            arg.positions = R.PosNames; 
            arg.order=R.SubPlotData.order; 
            
            arg = parseVarargin(varargin,arg); 
            
            %% fill up inputs arguments to match the arg.position size. 
            assert(numel(arg.positions)<=arg.ncols*arg.nrows,'Mismatch between number of requested positions and subplots');
            if isempty(arg.fig)
                fig_handle=figure; 
            else
                fig_handle=arg.fig; 
            end
            figure(fig_handle);
            if size(arg.ylim,1)==1
                arg.ylim=repmat(arg.ylim,numel(arg.positions),1);
            elseif size(arg.ylim,1)>numel(arg.positions)
                arg.ylim(~ismember(R.PosNames,arg.positions),:)=[]; 
            end
            if size(arg.xlim,1)==1
                arg.xlim=repmat(arg.xlim,numel(arg.positions),1);
            elseif size(arg.xlim,1)>numel(arg.positions)
                arg.xlim(~ismember(R.PosNames,arg.positions),:)=[]; 
            end
            if ischar(arg.xlabel) 
                arg.xlabel=repmat({arg.xlabel},numel(arg.positions),1);
            elseif numel(arg.xlabel)>numel(arg.positions)
                arg.xlabel(~ismember(R.PosNames,arg.positions),:)=[]; 
            end
            if ischar(arg.ylabel) 
                arg.ylabel=repmat({arg.ylabel},numel(arg.positions),1);
            elseif numel(arg.ylabel)>numel(arg.positions)
                arg.ylabel(~ismember(R.PosNames,arg.positions),:)=[]; 
            end
            if ischar(arg.title) 
                arg.title=repmat({arg.title},numel(arg.positions),1);
            elseif numel(arg.title)>numel(arg.positions)
                arg.title(~ismember(R.PosNames,arg.positions),:)=[]; 
            end
            if isempty(arg.order)
                arg.order=1:numel(arg.positions); 
            end
                        
            if strcmp(arg.hold,'off')
                clf
            end
            for i=1:numel(arg.positions)
                ix=arg.order(i);
                subplot(arg.nrows,arg.ncols,i)
                hold(arg.hold); 
                if nargin(func_plot)==1
                    D=R.getData(H,arg.positions{ix});
                    func_plot(D);
                else
                    [D,T]=R.getTimeseriesData(H,arg.positions{ix}); 
                    func_plot(T,D);
                end
                if ~isempty(arg.ylim)
                    ylim(arg.ylim(ix,:))
                end
                if ~isempty(arg.xlim)
                    xlim(arg.xlim(ix,:))
                end
                title(arg.title{ix})
                xlabel(arg.xlabel{ix}); 
                ylabel(arg.ylabel{ix}); 
            end
        end

        function cellid = showSegmentationAndResponseForRandomCell(R,MD,pos,channel,response,varargin)
            % pos: a position that must be in R.PosNames 
            % response: a timeseries resposse 
            % MD: a Metadata object
            % channel: the channel to use to read that stack
            
            arg.k=1;
            arg.movie = false; 
            arg.prctile=[6 99]; % to deal with the 5% dark corners
            arg.region = 'base'; 
            arg.ylim = [];
            arg.tags = {}; 
            arg = parseVarargin(varargin,arg);
            
            assert(ismember(pos,R.PosNames),'Position %s doesn''t exist',pos);
            assert(ismember(response,R.Header),'Response %s doesn''t exist',response);
            
            %% read Stack
            Stk = stkread(MD,'Position',pos,'Channel',channel); 
            Marker = unique(MD,'Marker','Channel',channel); 
            Marker = Marker{1}; 
            lbl = R.getLbl(pos); 
            Tchnl = MD.getSpecificMetadata('TimestampFrame','Position',pos,'Channel',channel); 
            Tchnl = cat(1,Tchnl{:}); 
            Stk = lbl.Reg.register(Stk,Tchnl); 
            lvl = prctile(Stk(randi(numel(Stk),10000,1)),arg.prctile); 
            
            %% add tags
            arg.tags = [arg.tags {'cellplayer',Marker}];
            
            resData=R.getTimeseriesData(response,pos); 
            
            cellid = randi(size(resData,2),arg.k);
            for i=1:arg.k
                if arg.movie
                    Vid = showVideoWithOutlines(lbl,Stk,lvl,resData,'cell_to_track',cellid(i),'T',Tchnl,'movie',true,...
                                                'Position',[0.0500 0.4903 0.4786 0.4104],'region',arg.region,'ylim',arg.ylim); 
                    name = sprintf('%s_%s_%s_%g.avi',pos,channel,response,cellid(i)); 
                    Vid=cat(4,Vid(2:end).cdata); 
                    R.stack2movie(Vid,name,'tags',arg.tags); 
                else
                    showVideoWithOutlines(lbl,Stk,lvl,resData,'cell_to_track',cellid(i),'T',Tchnl); 
                end
            end
                
            
        end
        
        function setCellLabelToSaveToDrive(R)
            for i=1:numel(R.PosNames)
                R.Lbl(i).saveToFile=true;
                R.Lbl(i).pth = R.pth;
                R.Lbl(i).posname = R.PosNames{i};
            end
        end
        
        function S = toStruct(R)
            % call the superclass method to start the transition to a
            % struct
            S = toStruct@MultiPositionResults(R);
            
            % add all the new fields
            for i=1:numel(R.Lbl) 
                if R.Lbl(i).saveToFile
                    R.Lbl(i).dumpToDrive;
                end
            end
            S.Lbl = R.Lbl;
            S.TimeVecs = R.TimeVecs;
            S.cell2exclude = R.cell2exclude;
        end
        
        function R = MultiPositionSingleCellResults(pth,reset) %constructor
            if nargin==0
                pth=''; 
                reset=false; 
            end
            if nargin==1
                reset=false; 
            end
            R@MultiPositionResults(pth,reset);
        end
        
        function [D,T] = getTimeseriesData(R,H,pos)
            if nargin==2
                pos=R.PosNames;  
            end
            D = getData(R,H,pos); 
            T = getT(R,H,pos); 
        end
        
        function [D,T] = getTimeseriesDataByProperties(R,H,varargin) 
            props = varargin(1:2:end);
            values = varargin(2:2:end);
            pos = getPosByProperty(R,props,values); 
            [D,T] = getTimeseriesData(R,H,pos);
        end
        
        function addTimeSeries(R,H,M,T,pos)
            % T can be a array of absolute times, a cell array of absolute
            % times for each position or a char to copy from
            if nargin<5 || isempty(pos)
                pos = R.PosNames; 
            end
            if ischar(T) % can be header names
                T=R.getTabs(T,pos);
            end
            add(R,H,M,pos);%inherited method
            ix = find(ismember(R.PosNames,pos));
            
            if numel(ix)==1 % single position
                R.TimeVecs(ix).(H)=T;
            else % there are multiple positions
                assert(numel(ix)==numel(T),'size of matrix must match the time vecotr');
                for i=1:numel(ix)
                    R.TimeVecs(ix(i)).(H)=T{i};
                end
            end
        end
        
        
        function T = getT(R,type,pos)
            if nargin==2 || isempty(pos)
                pos = R.PosNames; 
            end
            assertPositionExist(R,pos);
            T = getTabs(R,type,pos); 
            % temporarly make T/pos a cell if its not already one
            if ~iscell(pos), pos={pos}; end
            if ~iscell(T), T={T}; end
            fn = fieldnames(R.TimeVecs);
            for j=1:numel(pos)
                ix = ismember(R.PosNames,pos{j});
                minT=Inf;
                for i=1:numel(fn)
                    minT = min([minT; R.TimeVecs(ix).(fn{i})]);
                end
                T{j}=T{j}-minT;
                T{j}=T{j}*24*3600;
            end
            if numel(T)==1
                T=T{1}; 
            end
        end
        
        function T = getTabs(R,type,pos)
            if nargin==2
                pos = R.PosNames; 
            end
            ix = find(ismember(R.PosNames,pos));
            if numel(ix)==1
                T = R.TimeVecs(ix).(type);
            else
                T={R.TimeVecs(ix).(type);}; 
            end
        end
        
        function XY = getXY(R,pos,T)
            L = getLbl(R,pos); 
            if nargin==2, T=[]; end
            XY = L.getXY(T); %getXY is function of L		
        end
        
        
        function Lbl = setLbl(R,Lbl,pos,varargin)
            % accespts either a CellLabel (single or array) or a function
            % handle to calcualte the Lbl and then add it
            
            arg.redo = false; 
            arg = parseVarargin(varargin,arg); 
            
            
            if isa(Lbl,'function_handle') 
                assert(~iscell(pos),'Can only provide function handle for one position at a time'); 
                assertPositionExist(R,pos);
                % check to see if already exist, only if missing redo
                % calculations
                psname = R.PosNames; 
                if isnumeric(psname{1})
                    psname = cat(1,psname{:}); 
                end
                ix = find(ismember(psname,pos)); 
                if arg.redo || numel(R.Lbl) < ix || isempty(R.Lbl(ix)) 
                    Lbl=Lbl(); 
                else
                    Lbl=R.Lbl(ix); 
                end
            end
                    
            % do some checks
            assert(isa(Lbl,'CellLabel'),'second variable Lbl must be of Class CellLabel'); 
            assertPositionExist(R,pos); 
            if iscell(pos)
                assert(numel(Lbl)==numel(pos),'if position list is a cell array, Lbl must be of same size')
            else
                pos={pos}; 
            end
            psname = R.PosNames; 
            if isnumeric(psname{1})
                psname = cat(1,psname{:}); 
            end
            ix = ismember(psname,pos);
            
            R.Lbl(ix)=Lbl;
            
        end
        
        function lbl = getLbl(R,pos)
            assertPositionExist(R,pos); 
            ix = ismember(R.PosNames,pos); 
            lbl = R.Lbl(ix); 
        end
        
        function R = merge(Rvec,varargin)
            % TODO: extend support for merging object with different
            % TimeVec structures, for now I'm assuming they are the same!
            arg.prefix='%g_'; 
            arg = parseVarargin(varargin,arg); 
            arg.Results = MultiPositionSingleCellResults; 
            R = merge@MultiPositionResults(Rvec,arg); 
            
            R.Lbl=CellLabel.empty(1,0);
            for i=1:numel(Rvec)
                R.Lbl = [R.Lbl Rvec(i).Lbl]; 
            end
            
            allfields = {};
            for i=1:numel(Rvec)
                allfields = union(allfields,fieldnames(Rvec(i).TimeVecs)); 
            end
            for i=1:numel(Rvec)
                T=Rvec(i).TimeVecs; 
                missingflds = setdiff(allfields,fieldnames(T));
                for j=1:numel(missingflds)
                    T(1).(missingflds{j})=[]; 
                end
                T=orderfields(T,allfields); 
                if i==1
                    TimeVec=T; 
                else
                    TimeVec=[TimeVec T];  %#ok<AGROW>
                end
            end
            R.TimeVecs=TimeVec;               
        end
                
        
        function graphByPosition(R,Hx,Hy,func_plot,varargin)
            % function wrote by Yanfei. compared to plotByPosition, Plot not just Time traces, Hx, Hy
            % can be chosen to be something else. But data of Hx Hy has to
            % be numeric and aligned well. So when using addData method to add
            % data, make sure they are numeric and are in same orientation , plottable
            
            assert(~isempty(R.SubPlotData),'please set SubPlotData property before you call graphByPosition')
            
            arg.fig=[];
            arg.hold='off';
            arg.ylim=[];
            arg.xlim=[];
            arg.xlabel=R.SubPlotData.xlabel;
            arg.ylabel=R.SubPlotData.ylabel;
            arg.title=R.SubPlotData.title;
            arg.ncols=R.SubPlotData.ncols;
            arg.nrows=R.SubPlotData.nrows;
            arg.positions = R.PosNames;
            arg.order=R.SubPlotData.order;
            
            arg = parseVarargin(varargin,arg);
            
            %% fill up inputs arguments to match the arg.position size.
            assert(numel(arg.positions)<=arg.ncols*arg.nrows,'Mismatch between number of requested positions and subplots');
            if isempty(arg.fig)
                fig_handle=figure;
            else
                fig_handle=arg.fig;
            end
            figure(fig_handle);
            if size(arg.ylim,1)==1
                arg.ylim=repmat(arg.ylim,numel(arg.positions),1);
            elseif size(arg.ylim,1)>numel(arg.positions)
                arg.ylim(~ismember(R.PosNames,arg.positions),:)=[];
            end
            if size(arg.xlim,1)==1
                arg.xlim=repmat(arg.xlim,numel(arg.positions),1);
            elseif size(arg.xlim,1)>numel(arg.positions)
                arg.xlim(~ismember(R.PosNames,arg.positions),:)=[];
            end
            if ischar(arg.xlabel)
                arg.xlabel=repmat({arg.xlabel},numel(arg.positions),1);
            elseif numel(arg.xlabel)>numel(arg.positions)
                arg.xlabel(~ismember(R.PosNames,arg.positions),:)=[];
            end
            if ischar(arg.ylabel)
                arg.ylabel=repmat({arg.ylabel},numel(arg.positions),1);
            elseif numel(arg.ylabel)>numel(arg.positions)
                arg.ylabel(~ismember(R.PosNames,arg.positions),:)=[];
            end
            if ischar(arg.title)
                arg.title=repmat({arg.title},numel(arg.positions),1);
            elseif numel(arg.title)>numel(arg.positions)
                arg.title(~ismember(R.PosNames,arg.positions),:)=[];
            end
            if isempty(arg.order)
                arg.order=1:numel(arg.positions);
            end
            
            if strcmp(arg.hold,'off')
                clf
            end
            for i=1:numel(arg.positions)
                ix=arg.order(i);
                subplot(arg.nrows,arg.ncols,i)
                hold(arg.hold);
                Xs = R.getData(Hx,arg.positions{ix});
                Ys = R.getData(Hy,arg.positions{ix});
                func_plot(Xs, Ys);
                
                %there are lot of cases, like the layout of Xs, Ys; is Xs just an array
                %or a same size matrix; do Xs Ys have same size? are they numeric datatype, you need to
                %payattention to how you add new data to Results using addData. Make
                %sure you can plot them or if you want you can change the script here
                %to make Xs and Ys plottable, just go ahead
                
                
                if ~isempty(arg.ylim)
                    ylim(arg.ylim(ix,:))
                end
                if ~isempty(arg.xlim)
                    xlim(arg.xlim(ix,:))
                end
                title(arg.title{ix})
                xlabel(arg.xlabel{ix});
                ylabel(arg.ylabel{ix});
            end
        end
        
        function dataBeforeExclusion = storeOldData(R)
            if ~oldDataStored                
                dataBeforeExclusion = R.Data;
                R.dataBeforeExclusion = dataBeforeExclusion;
            end
        end

        
    end
    
end
