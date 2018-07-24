
classdef Results < handle
    
    properties
        % Type / Data pairs
        Header = {}; % Headers for all different data
        Data = {}; % actual data
        
        Movies = struct('name',{},'youtube',{},'tags',{}); % a list of movies create for this 
        
        Conclusions = 'missing'; % some final conclustion text.
        
        analysisScript = ''; % name of script that analyzes the raw data to create the Data & Headers 
        reportScript = '';
        
        reportPth = 'Report'; 
        
    end
    
    properties (Transient = true)
        pth
        verbose = true; 
    end
    
    methods (Static)
        function R = loadobj(S)
            R = Results; 
            R = reload(R,S); 
        end
        
    end
    
    methods 
        
        function S = saveobj(R)
            S = toStruct(R); 
        end
        
        function R = reload(R,S)
            R.Header = S.Header;
            R.Data = S.Data; 
            R.Conclusions = S.Conclusions; 
            R.analysisScript = S.analysisScript; 
            R.reportScript = S.reportScript;
            R.reportPth = S.reportPth; 
            % to maintain backward comp, only update if field exist
            if isfield(S,'Movies')
                R.Movies = S.Movies;
            end
        end
        
        function S = toStruct(R)
            S.pth = R.pth;
            S.Header = R.Header; 
            S.Data = R.Data; 
            S.Conclusions = R.Conclusions; 
            S.analysisScript = R.analysisScript; 
            S.reportScript = R.reportScript;
            S.reportPth=R.reportPth; 
            S.Movies = R.Movies; 
        end
        
        % constructor -> read from file
		function R = Results(pth,reset)
            if nargin == 0 || isempty(pth)
                return
            end
            if nargin==1
                reset=false;
            end
            R.pth = pth;
            if exist(fullfile(pth,'Results.mat'),'file') && ~reset
                s = load(fullfile(pth,'Results.mat'));
                S = toStruct(s.R); % just in case object was saved as Results / or a subclass of Results and not a struct
                R = reload(R,S);  % recreate the Results object (or part of object if its a subclass)
                
                % update Movie info
                try
                    R.updateMoviesInfoFromFSandYouTube; 
                catch e
                    warning('Cann''ot update YouTube movies, message: %s, moving on', e.message); 
                end
            end
		end
        
        function saveResults(R,pth)
           % save results to path pth
           if nargin==1
               pth = R.pth;
           end
           % make sure this is a pth based results (and not a merged
           % one...)
           assert(~isempty(pth),'can''t save Resutls with empty pth'); 
           
           % saving Conclusion as a seperate file to make it easier to read
           Conclusions = R.Conclusions;  %#ok<PROP,NASGU>

           save(fullfile(pth,'Results.mat'),'R','Conclusions', '-v7.3')
           fid=fopen(fullfile(pth,'Conclusions.txt'),'w'); 
           fprintf(fid,'%s',R.Conclusions);  
           fclose(fid); 
           save(fullfile(pth,'Conclusions.mat'),'Conclusions', '-v7.3')
           fid=fopen(fullfile(pth,'Classtype.txt'),'w'); 
           fprintf(fid,'%s',class(R)); 
           fclose(fid); 
           
        end
        
        function HTML = showMovie(R,name,varargin)
            arg.width=600; 
            arg.height=400; 
            arg.tags={}; 
            arg.htmltag=true; 
            arg = parseVarargin(varargin,arg); 
            
            HTML=''; 
            
            %% deal with a call by tags and not name
            if isempty(name) % then use tags
                if isempty(arg.tags)
                    warning('No name or tags, not showing a movie!')
                    return
                end
                tf = false(size(R.Movies)); 
                for i=1:numel(tf)
                    if ~isempty(R.Movies(i).tags)
                        tf(i) = all(ismember(arg.tags,R.Movies(i).tags));
                    end
                end
                NamesToShow = {R.Movies(tf).name}; 
                arg.tags={};
                arg.htmltag=false; 
                for i=1:numel(NamesToShow)
                    HTML=[HTML R.showMovie(NamesToShow{i},arg)];  %#ok<AGROW>
                end
                HTML=['<html>' HTML '</html>'];
                disp(HTML);
            else
                %% from here down we are working only by name
                
                if isempty(getappdata(0,'publishing')) || ~getappdata(0,'publishing')
                    %% figure out movie size and load it using Fiji (as avi)
                    filename = fullfile(R.pth,'Movies',name);
                    Mov = VideoReader(filename);
                    nFrames = Mov.NumberOfFrames;
                    cmd = sprintf('select=%s first=1 last=%g',filename,nFrames);
                    if ~exist('MIJ','class')
                        Miji;
                    end
                    MIJ.run('AVI...', cmd);
                    HTML=''; 
                else
                    %% create the html for a youtube movie
                    Names={R.Movies.name};
                    if nnz(ismember(Names,name))==0
                        error('Movie %s not found, please use stack2movie first',name)
                    end
                    if nnz(ismember(Names,name)>1)
                        error('Duplicate Movie names (%s) please check')
                    end
                    Codes={R.Movies.youtube};
                    if isempty(Codes{ismember(Names,name)})
                        % missing a youtube code. Try to add it:
                        addMissingYouTubeCodes(R)
                        Codes={R.Movies.youtube};
                        if isempty(Codes{ismember(Names,name)})
                            error('Could not get the YouTube code for this movie, was it uploaded?')
                        end
                    end
                    YouTubeID = Codes{ismember(Names,name)};
                    
                    HTML = sprintf('<iframe width="%g" height="%g" src="http://www.youtube.com/embed/%s?loop=1&amp;playlist=&amp;vq=hd720&amp;rel=0&amp;start=&amp;end=&amp;fs=0&amp;iv_load_policy=3&amp;modestbranding=1&amp;color=white" frameborder="0" allowfullscreen=""></iframe>',arg.width,arg.height,YouTubeID);
                    if arg.htmltag
                        HTML=['<html>' HTML '</html>'];
                        disp(HTML);
                    end
                end
            end
            
        end

        
        function stack2movie(R,Stk,name,varargin)
            
            %% input arguments
            clr=gray(256); 
            clr(1,:)=[0 0 0];
            arg.colormap=clr; 
            arg.framerate=7; 
            arg.extension = 'avi';
            arg.redo = false; 
            arg.tags = {}; 
            arg = parseVarargin(varargin,arg);
            
            assert(size(arg.colormap,1)<=256,'Cannot make a movie with more than 256 colormap')
            
            [~,~,ext] = fileparts(name); 
            if isempty(ext)
                name = [name '.' arg.extension]; 
            end
            % make sure there are no spaces
            name = regexprep(name,' ','_');
            
            if ~arg.redo && ismember(name,{R.Movies.name})
                warning('Movie %s already exist. To replace, change the redo flag to true\n',name)
                return
            end
            
            %% create folder if needed
            if ~exist(fullfile(R.pth,'Movies'),'dir')
                mkdir(fullfile(R.pth,'Movies'));
            end
            movfilename = fullfile(R.pth,'Movies',name);
            
           
            
            %% deal with duplicate
            if ismember(name,{R.Movies.name})
                warning('Movie %s already exist in names, will overwrite it',name) 
                delete(movfilename); 
                ix=find(ismember({R.Movies.name},name));
            else
                R.Movies(end+1).name=name; 
                ix=numel(R.Movies); 
            end
            
            %% transform Stack to a [N M 3 K] RGB if needed. 
            if ndims(Stk)==3 % assume a grayscale that need to be converted.  
                Stk = gray2ind(Stk,size(arg.colormap,1));
                Stk = permute(Stk,[1 2 4 3]);
                Stk = immovie(Stk,arg.colormap);
                Stk = cat(4,Stk.cdata);
            end
            
            % verify input
            assert(ndims(Stk)==4 && size(Stk,3)==3,'error, must supply a 3D stack of a 4D with NxMx3xK stacks'); 
            
            %% save to drive
            mov = VideoWriter(movfilename);
            mov.FrameRate=arg.framerate;
            mov.open;
            mov.writeVideo(Stk);
            mov.close;
            
            %% deal with tags
            
            R.Movies(ix).tags=arg.tags; 
            
            % create a tag string
            tag='"'; 
            for i=1:numel(arg.tags)
                tag=[tag arg.tags{i} ','];  %#ok<AGROW>
            end
            tag(end)='"';
            if numel(tag)==1
                tag(2)='"'; 
            end
            
            %% upload to youtube and get the code! 
            cmd=sprintf(['youtube-upload --title="%s" --description="%s" --tags=%s ' ... 
                         '--client-secrets=/home/rwollman/Google/my_client_secret.json --credentials-file=/home/rwollman/.youtube-upload-credentials.json "%s"'],...
                         name,R.pth,tag,movfilename);
            [err,yout]=system(cmd);        
            if err
                error('Could not upload to YouTube, message %s, filename %s',yout,movfilename)
            else
                splt = strsplit(yout,'\n'); 
                splt(cellfun(@isempty,splt))=[]; 
                R.Movies(ix).youtube=splt{end}; 
            end
        end
        
        function updateMoviesInfoFromFSandYouTube(R)
            if ~exist(fullfile(R.pth,'Movies'),'dir')
                return
            end
            lst = dir(fullfile(R.pth,'Movies/*.avi')); 
            onFS = {lst.name};  
            inR = {R.Movies.name}; 
            missing = setdiff(onFS,inR); 
            for i=1:numel(missing)
                R.Movies(end+1).name=missing{i}; 
            end
            R.addMissingYouTubeCodes; 
        end
        
        function addMissingYouTubeCodes(R)
           % function to look

           % see if there are any missing codes
           Codes = {R.Movies.youtube}; 
           if sum(~cellfun(@isempty,Codes))==numel(R.Movies)
               return
           end
           
           missing_ix = find(cellfun(@isempty,Codes)); 
           Names = {R.Movies.name}; 
           Names = Names(missing_ix); 
           
           %% start by looking up the google spreadsheet for this pth & Movie names
           system('wget --output-document /home/rwollman/bin/VideoList.csv https://docs.google.com/spreadsheets/d/1eFkyX68lmi0FPMqODp4YmOrqQxm4x9WCBHJruCmS_F0/export?gid=0\&format=csv')
           fid=fopen('/home/rwollman/bin/VideoList.csv'); 
           VideoList=textscan(fid,'%s %s %s','delimiter',',');
           fclose(fid);
           
           % keep only entires from this dataset
           VideoList{1}(~ismember(VideoList{2},R.pth))=[]; 
           VideoList{3}(~ismember(VideoList{2},R.pth))=[]; 
           
           if isempty(VideoList{1})
               warning('Did not find any youtube videos for this Results object')
           end
                      
           for i=1:numel(VideoList{1})
               ix = ismember(Names,VideoList{1}{i});
               if any(ix)
                   R.Movies(missing_ix(ix)).youtube=VideoList{3}{i};
               end
           end

        end
        
        function publish(R,varargin)
            %% supporss all output
            setappdata(0,'publishing',true)
            arg.redo = true; 
            arg = parseVarargin(varargin,arg);
            warning off %#ok<*WNOFF>
            % use the matlab publish command to create an html document
            % that summarises the resutls from that experiment.
            if arg.redo 
                mydoc = publish(R.reportScript,'outputDir',fullfile(R.reportPth));
                pathstr = fileparts(mydoc); 
                newfilename = fullfile(pathstr,'Report.html'); 
                copyfile(mydoc,newfilename)
            end
            warning on %#ok<*WNON>
            setappdata(0,'publishing',false)
        end
        
        function add(R,H,D,varargin)
            arg.redo = true; 
            arg = parseVarargin(varargin,arg); 
            if ~arg.redo &&  ismember(H,R.Header)
                error('Can''t add data with header %s - it already exist!',H); 
            elseif arg.redo && ismember(H,R.Header)
                setData(R,H,D); 
            else
                R.Header{end+1}=H;
                R.Data{end+1}=D;
            end
        end
        
        function deleteData(R,Headers) %added by Yanfei, use with caution
            % Headers should be in a cell array of strings
            % example: {'JNK_cyto', 'JNK_nuc', 'JNK_ring'}
            ix = ismember(R.Header, Headers);
            R.Data(ix) = [];
            R.Header(ix) = [];
            R.TimeVecs = rmfield(R.TimeVecs,Headers);
        end
        
        function D = getData(R,H)
            ix = find(ismember(R.Header,H)); 
            if isempty(ix)
                warning('Header %s does not exist',H); 
            end
            D=R.Data{ix};
        end
        
        function setData(R,H,D)
            if ~ischar(H)
                error('can only set a SINGLE header and it must be string!'); 
            end
            ix = find(ismember(R.Header,H), 1); 
            if isempty(ix)
                warning('Header %s does not exist, preforming an ADD operation!',H); 
                add(R,H,D); 
            else
                R.Data{ix}=D; 
            end
        end
    end
    
end