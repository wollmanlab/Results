function DatasetList = findAllDatasets(basepth,varargin)
% creates a struct array for DatasetList from filesystem 
% it has the following fields: Dataset / Date / User / Project / Results 

arg.report = ''; % options: 'only','missing'
arg.startdate = datenum('March-12-2013'); 
arg = parseVarargin(varargin,arg); 



DatasetList= struct('User',{},...
                        'Project',{},...
                        'Dataset',{},...
                        'Date',{},...
                        'Datenum',{},...
                        'Fullpath',{},...
                        'ResultExist',{},...
                        'Report',{},...
                        'ReportExist',{},...
                        'Summary',{},...
                        'AcqDatenum',{},...
                        'AcqDate',{},...
                        'AcqComments',{});  

if iscell(basepth)
    for i=1:numel(basepth)
       ds =  findAllDatasets(basepth{i},arg); 
       DatasetList = [DatasetList(:); ds(:)]; 
    end
    return
end
                    
%% get Userlist
Userlist = dir(basepth); 
Userlist = Userlist(3:end); 
Userlist(~[Userlist.isdir])=[]; 

%% for each user, create projects
for u=1:numel(Userlist)
    %%
    userpath = fullfile(basepth,Userlist(u).name); 
    ProjectListForUser = dir(userpath); 
    ProjectListForUser = ProjectListForUser(3:end); 
    ProjectListForUser(~[ProjectListForUser.isdir])=[]; 
    ProjectListForUser(cellfun(@(f) strcmp(f,'Delme'),{ProjectListForUser.name}))=[]; 
    
    %% for each project per user, find datasets
    for p = 1:numel(ProjectListForUser)
        projectpath = fullfile(userpath,ProjectListForUser(p).name); 
        DatasetsInProject = dir(projectpath); 
        DatasetsInProject = DatasetsInProject(3:end); 
        DatasetsInProject(~[DatasetsInProject.isdir])=[];
        
        for d = 1:numel(DatasetsInProject)
            
            
            %% find if a report exists
            DatasetPath = fullfile(projectpath,DatasetsInProject(d).name); 
            if exist(fullfile(DatasetPath,'Report','Report.html'),'file')
                Report = fullfile(DatasetPath,'Report','Report.html'); 
                ReportExist = true; 
            else
                Report = ''; 
                ReportExist = false;
            end
            
            if exist(fullfile(DatasetPath,'Results.mat'),'file'); 
                ResultExist = 1; 
            else
                ResultExist = 0; 
            end
            
            splt = regexp(DatasetsInProject(d).name,'_','split'); 
            AcqDate = splt{end}; 
            try
                AcqDatenum = datenum(AcqDate,'yyyymmmdd');
            catch  %#ok<CTCH>
                AcqDate = DatasetsInProject(d).date; 
                AcqDatenum = datenum('2012Mar01','yyyymmmdd');
            end
            % if empty - somethig didn't work out- change acq to last file
            % change by defaults
            if isempty(AcqDatenum)
                AcqDate = DatasetsInProject(d).date; 
                AcqDatenum = datenum('2012Mar01','yyyymmmdd');
            end
            
            %% if there are tif files in that dataset add it 
            % if there is a space in the path - skip it - not a legit
            % dataset
            if ~isempty(strfind(DatasetPath,' '))
                continue, 
            end
            DatasetPath = regexprep(DatasetPath,'(','\\('); 
            DatasetPath = regexprep(DatasetPath,')','\\)'); 
            [flag,str] = system(sprintf('find %s/ -name "*.tif"',DatasetPath));
            assert(~flag,'Couldn''t check for existance of tif files in dataset %s',DatasetPath); 
            if ~isempty(str) 
                %% combine all to a Struct
                DS = struct('User',Userlist(u).name,...
                    'Project',ProjectListForUser(p).name,...
                    'Dataset',DatasetsInProject(d).name,...
                    'Date',DatasetsInProject(d).date,...
                    'Datenum',DatasetsInProject(d).datenum,...
                    'Fullpath',DatasetPath,...
                    'ResultExist',ResultExist,...
                    'Report',Report,...
                    'ReportExist',ReportExist,...
                    'Summary','Summary Missing!',...
                    'AcqDatenum',AcqDatenum,...
                    'AcqDate',AcqDate,...
                    'AcqComments','None');
                
                
                DatasetList(end+1) = DS;  %#ok<AGROW>
            end
        end
    end
end

