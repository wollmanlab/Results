classdef ComputationalResults < Results
    
    properties
        pathtorepos = ''; % The directory of the local repository of codes 
        runpth = ''; % the directory of the results object in the TSCC
        runinput = {}; % cell of cells of inputs
        status = ''; 
        executeon = '';
        % processor, clusterhost, username, account, queue, time,
        % datalocation, remotedatalocation, keyfile, matlabroot, ppn,
        % attachedfiles
        % clustercomp contains the data for the   
        clustercomp = [];
        attachedfiles = {};

    end
    
    properties (SetAccess = protected)
        changesetid = [];
    end
    
    
    methods (Static)
        function R = loadobj(S)
            R = ComputationalResults; 
            R = reload(R,S); 
        end
    end
    
    methods
        
        
        function S = saveobj(R)
            S = toStruct(R); 
        end
        
        function R = reload(R,S)
            R = reload@Results(R,S); 

            if isfield(S,'pathtorepos') % is the fields exists load them, if not, 
                                     % they will be filled with default values
                                     % effectivly upcasting an object. 
                R.pathtorepos = S.pathtorepos;
            end
            if isfield(S,'changesetid')
                R.changesetid = S.changesetid;
            end

            if isfield(S,'runpth')
                R.runpth = S.runpth;
            end
            if isfield(S,'runinput')
                R.runinput = S.runinput;
            end
            if isfield(S,'status')
                R.status = S.status;
            end
            if isfield(S,'executeon')
                R.executeon = S.executeon;
            end
            if isfield(S,'clustercomp')
                R.clustercomp = S.clustercomp;
            end
            if isfield(S,'attachedfiles')
                R.attachedfiles = S.attachedfiles;
            end
        end
        
        function S = toStruct(R)
            % call the superclass method to start the transition to a
            % struct
            S = toStruct@Results(R);
            % add all the new fields
            S.pathtorepos = R.pathtorepos;
            S.changesetid = R.changesetid;
            S.runpth = R.runpth;
            S.runinput = R.runinput;
            S.status = R.status;
            S.executeon = R.executeon;
            S.clustercomp= R.clustercomp;   
            S.attachedfiles = R.attachedfiles;
        end
        
        function R = ComputationalResults(pth,reset,varargin)
            arg.pathtorepos = '';
            arg.changesetid = '';
            arg.runpth = '';
            arg.runinput = {};
            arg.status = '';
            arg.executeon = '';
            arg.clustercomp = [];
            arg.attachedfiles = {};
            arg = parseVarargin(varargin,arg); 
            if nargin==0
                pth ='';
                reset = false;
            end
            if nargin==1
                reset=false;
            end
            R@Results(pth,reset);
            % if nargin ==1 then just read the results object
            if isempty(R.pathtorepos) 
                R.pathtorepos = arg.pathtorepos;
            end
            if  isempty(R.changesetid)
                R.changesetid = arg.changesetid;
            end
            if isempty(R.runpth) 
                R.runpth = arg.runpth;      
            end
            if isempty(R.runinput)
                R.runinput = arg.runinput;     
            end
            if isempty(R.status)
                R.status = arg.status;
            end
            if isempty(R.executeon)
                R.executeon = arg.executeon;
            end
            if isempty(R.clustercomp)
                R.clustercomp = arg.clustercomp;
            end
            if isempty(R.attachedfiles)
                R.attachedfiles = arg.attachedfiles;
            end
        end
        
        function pth = get.pathtorepos(R)
            pth = R.pathtorepos;
        end
       
        function set.pathtorepos(R,pth)
            % verify that a .hg file exist in that path
            if ~exist(fullfile(pth,'.hg'),'dir')
                warning('Repo not found')
            end
            R.pathtorepos = pth;
        end
        
        function setID = get.changesetid(R)
            setID = R.changesetid;
        end        
        
        function runAnalysisFunction(R,varargin)
            arg.changesetid = '';
            arg.executeon = ''; 
            arg.clustercomp = [] ;
            arg.attachedfiles = {};
            arg = parseVarargin(varargin,arg); 
            if ~isempty(arg.executeon)
                R.executeon =arg.execute; 
            end
            % verify that analysis script was provided and it exist.
            analysisScriptFileName = fullfile(R.pathtorepos,R.analysisScript);
            assert(~isempty(analysisScriptFileName),'analysis script not provided');
            assert(exist(analysisScriptFileName,'file')>0,'analysis file doesn''t exist'); 
            
            %% If changesetid is empty, fill it up, either by creating a new one or by loading user defined on
            if isempty(R.changesetid)
                % If no changesetid provided, commit to create a new one,
                % assign it                 
                if isempty(arg.changesetid)
                    %TODO either change folder or figure out how to hg
                    %another place. 
                    cmd = sprintf(strcat('/opt/bin/hg commit -R ',R.pathtorepos,' -m "Computational Results Commit Before Run"')); 
                    [ok, cmdout]= system(cmd); 
                    cmd = sprintf(strcat('/opt/bin/hg id -R ',R.pathtorepos,' -i'));
                    [out,cmdout] = system(cmd);
                    R.changesetid = cmdout;
                    
                else
                    R.changesetid = arg.changesetid; 
                end
            end
                        
            
            %% If the run is on local server, then create a folder under the directory of runpth. If the run is on TSCC, then create a folder under the directory of runpth 
            if strcmp(R.executeon,'local')               
                if exist(fullfile(R.runpth,'runenv'))>0
                    cmd = ['rm -r ',fullfile(R.runpth,'runenv')];
                    system(cmd);
                end
                mkdir(fullfile(R.runpth,'runenv'));
                addpath(fullfile(R.runpth,'runenv'));
                % check the code from the repo with ChageSetID into runenv
                cmd = sprintf('/opt/bin/hg clone -r %s %s %s',strtrim(R.changesetid),strtrim(R.pathtorepos),strtrim(fullfile(R.runpth,'runenv')));
                ok = system(cmd); 
                assert(ok==0,'There was a problem with checking out files')
                
            elseif strcmp(R.executeon,'TSCC_Independent')
                if exist(fullfile(R.pth,'runenv'))>0
                    cmd = ['rm -r ',fullfile(R.pth,'runenv')];
                    system(cmd);
                end
                mkdir(fullfile(R.pth,'runenv'));   
                addpath(fullfile(R.pth,'runenv'));
                % check the code from the repo with ChageSetID into runenv
                cmd = sprintf('/opt/bin/hg clone -r %s %s %s',strtrim(R.changesetid),strtrim(R.pathtorepos),strtrim(fullfile(R.pth,'runenv')));
                ok = system(cmd); 
                assert(ok==0,'There was a problem with checking out files')
            
            else
                error('Error, executeon invalid. ');
            end
            
            if isempty(R.clustercomp)
                if isempty(arg.clustercomp)
                    error('no clustercomp provided');
                else
                    R.clustercomp = arg.clustercomp;
                end
            end
            
            if isempty(R.attachedfiles)
                if isempty(arg.attachedfiles)
                    error('no attachedfiles provided');
                else
                    R.attachedfiles = arg.attachedfiles;
                end
            end
            %% 
            R.execute; 
            
 
           
        end
        
        function execute(R,varargin)
            arg.executeon = '';
            arg = parseVarargin(varargin,arg); 
            
            %% execute
            if ~isempty(arg.executeon)
                R.executeon =arg.execute; 
            end
            R.status = 'executing';
            switch R.executeon
                case 'local'
                    Rtmp = ComputationalResults; 
                    Rtmp.Input=R.runinput; 
                    Rtmp = feval(R.analysisScript,Rtmp); 
                case 'TSCC'
                    error('TSCC execution is not supported for a ComptationalResults class, only for MultiRunComputaitonalResults')
            end
            
            %% gather and save results
            switch R.executeon
                case 'local'
                    R.Data=Rtmp.Data; 
                    R.Header=Rtmp.Header; 
                case 'TSCC'
                    error('TSCC execution is not supported for a ComptationalResults class, only for MultiRunComputaitonalResults')
            end
        end
        


    end
    
    
    
end