classdef MultiRunComputationalResults < ComputationalResults
    properties
        id= {};
    end
    
    properties (Transient = true)
        j % a sibmitted job
                
    end
                           
    methods (Static)
        function R =loadobj(S)
            R = MultiRunComputationalResults;
            R = reload(R,S);
        end    
    end
    methods
        
        function S = saveobj(R)
            S = toStruct(R);
        end
        function R = reload(R,S)
            R = reload@ComputationalResults(R,S); 
            if isfield(S,'id')
                R.id = S.id;             
            end
        end
        function S = toStruct(R)
            % call the superclass method to start the transition to a
            % struct
            S = toStruct@ComputationalResults(R);
            % add all the new fields
            S.id = R.id;

        end
        % constructor
        function R = MultiRunComputationalResults(pth,reset,varargin)
            arg.pathtorepos = '';
            arg.changesetid = '';
            arg.runpth = '';
            arg.runinput = {};
            arg.status = '';
            arg.executeon = '';
            arg.id = {};
            arg.clustercomp = [];
            arg.attachedfiles = {};
            arg = parseVarargin(varargin,arg);
            if nargin==0
                pth = '';
                reset = false;
            end
            if nargin==1
                reset=false;
            end
            R@ComputationalResults(pth,reset,'pathtorepos',arg.pathtorepos,'changesetid',arg.changesetid,'runpth',arg.runpth,'runinput',arg.runinput,'status',arg.status,'executeon',arg.executeon,'clustercomp',arg.clustercomp,'attachedfiles',arg.attachedfiles);  
            if isempty(R.id)
                R.id = arg.id;
            end
        end
        
        % overall all the save & load stuff
        
        % Overload 
        function execute(R,varargin)
            % rewrite the execute stuff such that Rtmp is in an array. in
            % local we run that in a parfor loop in TSCC we submit multiple
            % tasks in a job and it calls merge in the end to merge then by
            % Runid
            arg.executeon = '';
            arg = parseVarargin(varargin,arg);

            %% execute
            if ~isempty(arg.executeon)
                R.executeon = arg.execute;
            end
            R.status = 'executing';
            switch R.executeon
                case 'local'
                    parallel.defaultClusterProfile('local');
                    RtmpArray(length(R.id)) = MultiRunComputationalResults; 
                    % Run simulations 
                    for i =1:length(RtmpArray)
                        Rtmp = RtmpArray(i);
                        Rtmp.pth = fullfile(R.runpth,R.id{i});
                        Rtmp.runpth = Rtmp.pth;
                        Rtmp.id = {R.id{i}};
                        rRunInput = R.runinput{i};
                        Rtmp.runinput = {rRunInput};
                        C = strsplit(R.analysisScript,'.');
                        analysisFunc = strcat('@',C{1});
                        Rtmp = feval(str2func(analysisFunc),Rtmp);
                        RtmpArray(i) = Rtmp;
                    end
                    % merge the RtmpArray to R
                    Rmerge = RtmpArray.merge;
                    R.id = Rmerge.id;
                    R.runinput = Rmerge.runinput;
                    R.Data = Rmerge.Data;
                    
                % TSCC independent job
                case 'TSCC_Independent' 
                    parallel.defaultClusterProfile('GenericProfile1');
                    RtmpStructArray = cell(length(R.id),1);
                    % Run simulations
                    for i =1:length(RtmpStructArray)
                        Rtmp = MultiRunComputationalResults; 
                        Rtmp.pth = fullfile(R.runpth,R.id{i});
                        Rtmp.runpth = Rtmp.pth;
                        Rtmp.id = {R.id{i}};
                        rRunInput = R.runinput{i};
                        Rtmp.runinput = {rRunInput};
                        % assign the temporary object
                        S = toStruct(Rtmp);
                        RtmpStructArray{i} = S; 
                    end
                    %  submit the job
                    % create a cluster
                    clustercomp = R.clustercomp;
                    cluster = getCluster(clustercomp.username,clustercomp.account,clustercomp.clusterhost,clustercomp.ppn,clustercomp.queue,clustercomp.time,clustercomp.datalocation, clustercomp.remotedatalocation,clustercomp.keyfile,clustercomp.matlabroot);
                    R.j = createJob(cluster);
                    % transfer the attached files
                    for i =1:length(R.attachedfiles)
                        attachedfile = fullfile(R.pth,'runenv',R.attachedfiles{i});
                        R.attachedfiles{i} = attachedfile; 
                    end
                    % attach the files related to the class and superclasses of the results
                    
                    objectFilesAttach = {};
                    objectFilesAttach{1} = which(class(R));
                    outSuperclassesCell = superclasses(R);
                    for i =1:length(outSuperclassesCell)-1
                        objectFilesAttach = {objectFilesAttach{:},which(outSuperclassesCell{i})};
                        
                    end
                    R.j.AttachedFiles = {R.attachedfiles{:},objectFilesAttach{:}};
                    
                    %j.AttachedFiles = R.attachedfiles;
                    C = strsplit(R.analysisScript,'.');
                    analysisFunc = strcat('@',C{1});   
                    t = createTask(R.j,str2func(analysisFunc),1,RtmpStructArray);
                    % submit the job 
                    submit(R.j);
                    % 
                    flag = 0;
                    while flag == 0
                        if ~isempty(R.j.StartTime)
                            flag =1;
                        end
                    end
                    pause(60);
                    %taskMap(R.j,'RefreshPeriod',20);
            end     
        end
        
        % merge the simulation results from individual runs into one big
        % result object
        function R = merge(Rvec,varargin)
                %instantiate a new results object
                R = MultiRunComputationalResults;            
                % merge the data under the same header
                if length(Rvec) == 1
                    R.Header = Rvec.Header;   
                else
                    R.Header = Rvec(1).Header;  
                end
                % initialize the cells of data to read
                idArray = {};
                runinputArray = {};
                RDataCell = cell(1,length(R.Header));
                for i =1:length(RDataCell)
                    RDataCell{i} = {};                   
                end
                
                if length(Rvec) == 1
                    Rtmp = Rvec;
                    % read the id 
                    idArray = {idArray{:}; Rtmp.id{1}};
                    % read the runinput
                    runinputArray = {runinputArray{:}; Rtmp.runinput{1}};
                    % read the data from the small results object
                    currData = Rtmp.Data;
 
                    % go through each of the array of results to assign the data 
                    for j =1:length(currData)
                        currDataCell = RDataCell{j};
                        RDataCell{j} = {currDataCell{:};currData{1}};
                    end
                    R.id = idArray;
                    R.runinput = runinputArray;
                    R.Data = RDataCell;               
                    
                else
                    for i =1:length(Rvec)
                        Rtmp = Rvec(i);
                        % read the id 
                        idArray = {idArray{:}; Rtmp.id{1}};
                        % read the runinput
                        runinputArray = {runinputArray{:}; Rtmp.runinput{1}};
                        % read the data from the small results object
                        currData = Rtmp.Data;

                        % go through each of the array of results to assign the data 
                        for j =1:length(currData)
                            currDataCell = RDataCell{j};
                            RDataCell{j} = {currDataCell{:};currData{1}};
                        end
                    end
                    % 
                    R.id = idArray;
                    R.runinput = runinputArray;
                    R.Data = RDataCell;                    
                end
  
        end
        
        % This function gather in the results 
        function  gather(R,varargin)
            arg.wait = false; 
            arg.gatherpartial = false; 
            
            if ~isfinished(R)
                if arg.wait
                    R.waitForExecution; 
                elseif ~arg.gatherpartial
                    error('Cannot gather, execution is not finished')
                end
            end
            
            % create a job for gathering all task based Results.mat
            % the gatherfunc (always the same, saved in repo outside
            % Resutls classes) will get two input arguments, a pth and a
            % cell array of ids and will return a cell array of structs
            % that we will get by calling fetchOutput
            clustercomp = R.clustercomp;
            cluster = getCluster(clustercomp.username,clustercomp.account,clustercomp.clusterhost,clustercomp.ppn,clustercomp.queue,clustercomp.time,clustercomp.datalocation, clustercomp.remotedatalocation,clustercomp.keyfile,clustercomp.matlabroot);
            gatherJob = createJob(cluster);
            gatherJob.AttachedFiles = R.attachedfiles;
            inputCell = {R.runpth,R.id};
            % add gatherfunc with R.runpth and R.ids
            t = createTask(gatherJob,@GatherSMCMCDataFromTSCC,1,inputCell);
            submit(gatherJob);
            wait(gatherJob);
            CofS = gatherJob.fetchOutputs; 
            
            assert(numel(CofS)==numel(R.id),'error in retrival of output')
            CofS(cellfun(@isempty,CofS))=[]; 
            RtmpArray(numel(CofS)) = MultiRunComputationalResults ;
            for i=1:numel(CofS)
                S = CofS{i};
                RtmpArray(i).reload(S{1});
            end
            
            Rmerge = RtmpArray.merge;
            % re-assign the id, runinput and Data
            R.Header = Rmerge.Header;
            R.id = Rmerge.id;
            R.runinput = Rmerge.runinput;
            R.Data = Rmerge.Data; 
        end
        
        
        function waitForExecution(R,varargin)
            arg.verbose = false; 
            while ~R.isfinished
                pause(10)
            end
        end
        
        % This function finishes 
        function flag = isfinished(R)
            if isempty(R.j.FinishTime)
                flag = 0;
            else
                flag = 1; % job is finished more than 1 min ago
            end
        end
        

        
    end
        
end