classdef ResultsCollection < MultiPositionSingleCellResults
    
    properties
		PossibleBasePaths = {'/synology/data/Images/','/data2/Images/','/data3/Images/'}; % will default to saving in the last entry. 
		                                                                                  % will use the rest to find Results objects. 
		Username
		Project
		Datasets ={};
		
		CollectionName 
	end
    
	methods (Static)
        function R = loadobj(S)
            R = ResultsCollection; 
            R = reload(R,S); 
        end
    end
	
    methods
		% overload add method to throw an error since it shouldn't be called ever....
		function add(R) %#ok<MANU>
			error('Cann''t add data directly to a collection. Please add the data to any of the Results in the collection')
		end
		
		% overload addTimeSeries method to throw an error since it shouldn't be called ever....
		function addTimeSeries(R) %#ok<MANU>
			error('Cann''t add data directly to a collection. Please add the data to any of the Results in the collection')
		end
	
		function R = ResultsCollection(pth,varargin)
			arg.project='';
			arg.user='';
			arg.datasets={};
			arg.name='';
			
			arg = parseVarargin(varargin,arg);
			if nargin==0
				pth='';
			end
			% call the base result. If result is empty it will do nothing,
			% it it is not it will create the object. 
			R@MultiPositionSingleCellResults(pth);
			% if pth is empty, it means that we are recreating the object
			% from scratch. 
			if isempty(pth)
				R = createFromScratch(R,arg.username,arg.project,arg.datasets,arg.name,'type',arg.type);
			end
		end
		
		function R = createFromScratch(R,username,project,datasets,name,varargin)
			arg.date=datestr(now,'yyyymmmdd');
			arg.basepth=R.PossibleBasePaths{end};
			arg = parseVarargin(varargin,arg); 
			
			% make sure all inputs make sense
			assert(~isempty(username),'Please provide username');
			assert(~isempty(project),'Please provide project name');
			assert(~isempty(datasets) || ~iscell(datasets),'Please provide dataset list as cell array');
			assert(~isempty(name),'Please provide name for the collection');
			pth=fullfile(arg.basepth,username,project,sprintf('%s_%s',name,arg.date));
			if ~exist(pth,'dir')
				mkdir(pth);
			end
			R.pth=pth;
			for i=1:numel(arg.datasets)
				% find in which of the possible base paths the actual
				% dataset really is. 
				possiblePaths=cellfun(@(b) fullfile(b,arg.user,arg.project,arg.datasets{i}),R.PossibleBasePaths,'uniformoutput',0);
				correctPath=cellfun(@(p) exist(p,'dir'),possiblePaths);
				assert(nnz(correctPath)==1,'Did not find a single correct path for user: %s, project: %s, dataset: %s', arg.user,arg.project,arg.datasets{i});
				datasetpth=possiblePaths(possiblePaths);
				Rvec(i)=MultiPositionSingleCellResults(datasetpth); %#ok<AGROW>
			end
			Rmerge=merge(Rvec);
			R.Header=Rmerge.Header;
			R.Data=Rmerge.Data;
			R.PosNames=Rmerge.PosNames;
			R.PosProperties=Rmerge.PosProperties;
			R.ByPosition=Rmerge.ByPosition;
			R.TimeVec=Rmerge.TimeVec;
		end
		
		function S = saveobj(R)
            S = toStruct(R); 
        end
			
		function S = toStruct(R)
            % Unlike all other Results subclasses we are NOT calling the
            % base toStruct method.  (S = toStruct@Results(R));
            
            % add fields information only for the stuff needed to reacreate
            % from scratch
            S.Username = R.Username;
            S.Project = R.Project;
            S.Datasets = R.Datasets;
			S.CollectionName=R.CollectionName; 
		end
		
		function R = reload(R,S)
			assert(exist(R.pth,'dir'),'Call for reload for ResultsCollection must happen when the pth is set')
			splt = strsplit(R.pth,'_');
			date = splt{end};
			date = regexprep(date,'/','');
			basepth='';
			
			% this will create the object from all the other datasets
			R = createFromScratch(R,R.Username,R.Project,R.Datasets,R.CollectionName,'date',date,'basepth',basepth); 
			
			% manulally get fields that usually Results would get since we
			% are not calling the superclass reload; 
			R.Conclusions = S.Conclusions;
			R.analysisScript = S.analysisScript;
			R.reportScript = S.reportScript;
			R.reportPth = S.reportPth;
		end
    end
end