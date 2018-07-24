classdef MultiPositionResults < Results
    
    properties
       PosNames = {}; 
       ByPosition = []; % Flag that says which of the headers is by positon.
       
       PosProperties
       
    end
    
    properties (Dependent = true)%If false, property value is stored in object. If true, property value is not stored in object. The set and get functions cannot access the property by indexing into the object using the property name.
       Np %number of positions
    end
    
    
    methods (Static) %I guess loadobj has to be static method
        function R = loadobj(S)
            R = MultiPositionResults; %this is why it is overloaded, it was R = Results
            R = reload(R,S); 
        end
    end
    
    methods
        

        
        function R = merge(Rvec,varargin)%merge a vector of Rs
            
            arg.prefix='%g_'; %default setting
            arg = parseVarargin(varargin,arg); % function arg=parseVarargin(var,arg)
% parses any additional input argument that came in property/value pair
% arg is a struct where field names ar eproperty names and values ar
% edefault values. 
            
            %% deal with position names and Labels
            if isfield(arg,'Results') && ~isempty(arg.Results)
                R = arg.Results; 
            else
                R = MultiPositionResults;
            end
            Pos={}; 
            for i=1:numel(Rvec)
                oldposnames = cellfun(@(nm) sprintf([arg.prefix nm],i),Rvec(i).PosNames,'uniformoutput',0); 
                Pos = [Pos oldposnames(:)']; %#ok<AGROW>
            end
            R.PosNames=Pos; 
            
            %% deal with data that is NOT by position 
            % basically just add a prefix for each header and thats it. 
            % no real point in merging data that is not by position.
            NonPosHeaders = {}; 
            NonPosData = {}; 
            for i=1:numel(Rvec)
                PossHeaders = Rvec(i).Header(~(Rvec(i).ByPosition)); 
                PossData = Rvec(i).Data(~(Rvec(i).ByPosition));
                PossHeaders = cellfun(@(nm) sprintf([arg.prefix nm],i),PossHeaders,'uniformoutput',0);
                NonPosHeaders=[NonPosHeaders(:)' PossHeaders(:)']; 
                NonPosData=[NonPosData(:)' PossData(:)']; 
            end
            R.ByPosition = false(1,numel(NonPosData)); 
            
            %% Deal with data that IS by position
            SharedPosHeader={}; 
            for i=1:numel(Rvec)
                SharedPosHeader=union(SharedPosHeader,Rvec(i).Header(Rvec(i).ByPosition>0)); 
            end
            
            for j=1:numel(SharedPosHeader)
                R.Header=[R.Header SharedPosHeader{j}]; 
                Data={}; 
                for i=1:numel(Rvec)
                    if ismember(SharedPosHeader{j},Rvec(i).Header)
                        D=Rvec(i).getData(SharedPosHeader{j}); 
                    else
                        D=cell(Rvec(i).Np,1); 
                    end
                    Data=[Data(:); D];  
                end
                R.Data{end+1}=Data; 
            end
            R.ByPosition=[R.ByPosition true(1,numel(SharedPosHeader))];
            
            %% merge PosProperties
            
            allfields = {}; 
            for i=1:numel(Rvec)
                allfields = union(allfields,fieldnames(Rvec(i).PosProperties));
            end
            for i=1:numel(allfields)
                newPP.(allfields{i})=[]; 
                for j=1:numel(Rvec)
                    if isfield(Rvec(j).PosProperties,allfields{i})
                        newPP.(allfields{i})=[newPP.(allfields{i}) Rvec(j).PosProperties.(allfields{i})]; 
                    else
                        newPP.(allfields{i})=[newPP.(allfields{i}) nan(1,Rvec(j).Np)]; 
                    end
                end
            end
            R.PosProperties=newPP;  
        end
        
        function S = saveobj(R)
            S = toStruct(R); 
        end
        
        function R = reload(R,S)
            R = reload@Results(R,S); 
            if isfield(S,'PosNames') % is the fields exists load them, if not, 
                                     % they will be filled with default values
                                     % effectivly upcasting an object. 
                R.PosNames = S.PosNames;
                R.ByPosition = S.ByPosition;
                R.PosProperties = S.PosProperties;
            end
        end
        
        function S = toStruct(R)
            % call the superclass method to start the transition to a
            % struct
            S = toStruct@Results(R);
            
            % add all the new fields
            S.PosNames = R.PosNames;
            S.ByPosition = R.ByPosition;
            S.PosProperties = R.PosProperties;
        end
        
        function R = MultiPositionResults(pth,reset)%constructor
            if nargin==0
                pth='';
                reset=false; 
            end
            if nargin==1
                reset=false; 
            end
            R@Results(pth,reset);
        end
        
        function Np = get.Np(R)%number of positions
            Np = numel(R.PosNames); 
        end
        
        function P = getProperty(R,prop,pos)
            assert(isfield(R.PosProperties,prop),'Requested property %s not found',prop);
            P = R.PosProperties.(prop);
            if nargin==3
                assertPositionExist(R,pos)
                P = P(ismember(R.PosNames,pos));
            end
            
        end
        
        function setProperty(R,prop,value,pos)
            if nargin==3
                pos = R.PosNames; 
            end
            if ischar(value)
                value={value}; 
            end
            ixpos = ismember(R.PosNames,pos); 
                if numel(value) == 1;
                    R.PosProperties.(prop)(ixpos)=value; 
                end
            
        end
        
        function setAllNonStandardProperties(R,MD,pos)
            %% add all properties to all positions
            if nargin==2
                pos=R.PosNames; 
            end
            if iscell(pos)
                for i=1:numel(pos)
                    setAllNonStandardProperties(R,MD,pos{i}); %recursive, if there are two arguments, it means multiple positions
                end
                return
            end
            % add all new properties to single position         
            Props =  MD.NewTypes;%New types that are not in the default types list
            for j=1:numel(Props)
                try
                    tmp=unique(MD,Props{j},'group',pos);
                catch
                    tmp = MD.getSpecificMetadata(Props{j},'group',pos);
                    if iscell(tmp{1}) 
                        % remove empties 
                        tmp(cellfun(@isempty,tmp))=[]; 
                        assert(all(cellfun(@(c) numel(c)<=1,tmp)),'Something wrong with the properties, please check'); 
                        tmp = cellfun(@(c) c{1},tmp,'uniformoutput',0); 
                    end
                    if isnumeric(tmp{1})
                        tmp = unique(cat(1,tmp{:}));
                    end
                end
                if numel(tmp)==1;%only proprties that are common between all frames in position
                    if isnumeric(tmp(1)) || islogical(tmp(1))
                        setProperty(R,Props{j},tmp(1),pos);
                    elseif iscell(tmp{1})
                        setProperty(R,Props{j},tmp{1},pos);
                    elseif iscell(tmp)
                        setProperty(R,Props{j},tmp,pos);
                    end
                end
            end%end of for loop
        end%end of setAllNonStandardProperties function

function cleanEmptyProps(R,MD)
        Props =  MD.NewTypes;%New types that are not in the default types list
        for j=1:numel(Props)
            if isempty(R.PosProperties.(Props{j}))
                R.PosProperties = rmfield(R.PosProperties,Props{j})
            end
        end
end
        
        function pos = getPosByProperty(R,varargin)

            if numel(varargin)==2 && iscell(varargin{1}) && iscell(varargin{2})
                props = varargin{1};
                values = varargin{2}; 
            else
                props = varargin(1:2:end);
                values = varargin(2:2:end);
            end
            
            assert(numel(props)==numel(values),'Must suplly pairs of props / value'); 
            
            tf = true(R.Np,1); 
            for i=1:numel(props)
                P = getProperty(R,props{i});
                P=P(:); 
                tf = tf & ismember(P,values{i}); 
            end
            pos = R.PosNames(tf); 
        end
        
        function D = getDataByProperty(R,H,varargin)
            props = varargin(1:2:end);
            values = varargin(2:2:end);
            pos = getPosByProperty(R,props,values); 
            D = getData(R,H,pos); 
        end
        
        function set.PosNames(R,PN)
            R.PosNames=PN; 
        end
        
        %@OVERLOAD getData from Results
        function D = getData(R,H,pos)
            
            % make sure requested H exist
            assert(any(ismember(R.Header,H)),'Requested data type: %s doesn''t exist',H); 
            
            % get data from Results
            D = getData@Results(R,H);
            
            % if user asked for specific position and that position has
            % data by label
            if nargin==3 
                if ~R.ByPosition(ismember(R.Header,H))
                    warning('Asked for data of type %s which is not set to be by position!',H)
                    return
                end
                assertPositionExist(R,pos);
                if iscell(pos)
                    D=D(ismember(R.PosNames,pos));
                else
                    D=D{ismember(R.PosNames,pos)}; 
                end
            end
            
        end
        
        %@OVERLOAD setData from Results
        function setData(R,H,D,pos)
            
            if isempty(R.PosNames)
                error('Cann''t set/add data before setting PosNames!'); 
            end
            
            
            ix = find(ismember(R.Header,H), 1);
            if isempty(ix)
                warning('Header %s does not exist, preforming an ADD operation!',H);
                if nargin==3
                    add(R,H,D); 
                else
                    add(R,H,D,pos); 
                end
                return
            end
            
            % if position is not supplied add whole data and possibly mark
            % as by position dependeing on the data shape
            if nargin==3 
                setData@Results(R,H,D);
                % if data is a cell array of the same size as PosName
                % assumes its a position based result and mark it that way
                if numel(D)==numel(R.PosNames) && iscell(D)
                    R.ByPosition(ix)=true;
                end
                
            % if position was provided, see if data is byposition    
            else
                % make sure the data type is by position and that the
                % positions requested all exist. 
                assert(R.ByPosition(ix)==1,...
                    'Data to set: %s is not by position, yet a position %s was supplied!',H,pos); %added condition ==1. AOY
                assertPositionExist(R,pos); 
                ixpos = ismember(R.PosNames,pos); 
                if iscell(pos) 
                    assert(iscell(D) && numel(D) ==numel(pos),...
                           'if pos is a cell array than D must be a cell array with the same number of items!'); 
                    R.Data{ix}(ixpos)=D;
                else
                    R.Data{ix}{ixpos}=D; 
                end
            end
        end
        
        %@OVERLOAD add from Results, add pos 
        function add(R,H,D,pos,varargin)
            arg.redo = true; % not used here, it is not used in this method, but in the Results add method
            arg = parseVarargin(varargin,arg); 
            
            %% decide whether its a regular add or a position based add
            % position based add is either if pos is specified of if D
            % "looks" like it is ment to be all positions. 
            
            
            if (nargin>3 && ~isempty(pos)) || (iscell(D) && numel(D)==R.Np)
            % If this is a position based data
                if ismember(H,R.Header) 
                    % get the whole cell array and replace specific
                    % positions if needed or just replace the while thing
                    if (nargin>3 && ~isempty(pos)) % pos exists as input
                        assertPositionExist(R,pos); 
                        DD = getData(R,H); % get the whole cell array; 
                        if iscell(pos)
                            assert(numel(pos)==numel(D),'you must supply an array of data same size of position cell'); 
                            DD(ismember(R.PosNames,pos))=D;
                        elseif ischar(pos)
                            DD{ismember(R.PosNames,pos)}=D; 
                        end
                    else % pos doesn't exist, but D seems to be for all positions, just ignore current data 
                        DD=D; 
                    end
                    % add to existing header
                    R.Data{ismember(R.Header,H)}=DD;
                    R.ByPosition(ismember(R.Header,H))=true;
                else
                % H is not in the dataset, need to add new, first make an empty cell array
                    DD = cell(R.Np,1); 
                    % now either populate just a few cells based on pos
                    if (nargin>3 && ~isempty(pos)) % pos existed
                        if iscell(pos)
                            assert(numel(pos)==numel(D),'you must supply an array of data same size of position cell'); 
                            DD(ismember(R.PosNames,pos))=D;
                        elseif ischar(pos)
                            DD{ismember(R.PosNames,pos)}=D; 
                        end
                    else % or if pos is empty, but D is a cell array size of Np (based on if statement above)
                        DD=D; 
                    end
                    R.Header{end+1}=H; 
                    R.Data{end+1}=DD; 
                    R.ByPosition(end+1)=true; 
                end
            else
                % we are not adding a position based data, just use Results
                add@Results(R,H,D,arg);
            end
        end
        
        function deleteData(R,Headers) %added by Yanfei, use with caution
            % Headers should be in a cell array of strings
            % example: {'JNK_cyto', 'JNK_nuc', 'JNK_ring'}
            ix = ismember(R.Header, Headers);
            R.Data(ix) = [];
            R.Header(ix) = [];
            R.ByPosition(ix) = []; % the only property that requires this overloading
            R.TimeVecs = rmfield(R.TimeVecs,Headers);
        end
    end
    
    
    methods (Access = protected)
        function assertPositionExist(R,pos)
            if ~ischar(pos) && ~iscell(pos)
                warning('not checking position name - make sure you got it right!'); 
                return
            end
            df = setdiff(pos,R.PosNames);
            if ~isempty(df)
                posstr = df{1};
                for i=2:numel(df)
                    posstr=[posstr ' / ' df{i}];  %#ok<AGROW>
                end
                error('Position requested: %s are not in PosNames!',pos);
            end
        end
    end
    
end
