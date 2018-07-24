classdef MultiPositionWoundResults < MultiPositionResults
   properties
       PIVlbl
   end
   
   properties (Dependent = true) 
       TimeVecs
   end
   
    methods (Static)
        function R = loadobj(S)
            R = MultiPositionWoundResults; 
            R = reload(R,S); 
        end
    end
   
    
   methods
        
        function S = saveobj(R)
            S = toStruct(R); 
        end
        
        function R = reload(R,S)
            R = reload@MultiPositionResults(R,S); 
            if isfield(S,'PIVlbl') % is the fields exists load them, if not, 
                                % they will be filled with default values
                                % effectivly upcasting an object. 
                %R.WoundLbl = S.WoundLbl;
                R.PIVlbl = S.PIVlbl;
%                 R.TimeVecs = S.TimeVecs;
            end
        end
        
        function S = toStruct(R)
            % call the superclass method to start the transition to a
            % struct
            S = toStruct@MultiPositionResults(R);
            
            %S.WoundLbl = R.WoundLbl;
            S.PIVlbl = R.PIVlbl;
%             S.TimeVecs = R.TimeVecs;
        end
        
        function R = MultiPositionWoundResults(pth,reset) %constructor
            if nargin==0
                pth=''; 
                reset=false; 
            end
            if nargin==1
                reset=false; 
            end
            R@MultiPositionResults(pth,reset);
        end
        
        %function TimeVecs = get.TimeVecs(R)
        %    TimeVecs = getTimeVecs(R);
        %end
        
        function TimeVecs = get.TimeVecs(R)
             for i=1:R.Np             
                     pos = R.PosNames{i};
                     ix = ismember(R.PosNames,pos);
                 if ~isempty(R.PIVlbl{ix})
                    TimeVecs{ix}=cell2mat(R.PIVlbl{ix}.Tvec);
                 else
                    TimeVecs{ix} = [];
                 end
             end;
        end
        
         function T = getTimeVecByPosition(R,pos)
             assertPositionExist(R,pos); 
             ix = ismember(R.PosNames,pos);  
             T =  R.TimeVecs{ix};
         end
        
        
        
        function setPIVLbl(R,PLbl,pos,varargin)
             if numel(R.PIVlbl)<R.Np
                R.PIVlbl{R.Np}=[];
             end
            
            arg.redo = false; 
            arg = parseVarargin(varargin,arg); 
            
            % do some checks
            assert(isa(PLbl,'PIVLbl'),'second variable PLbl must be of Class PIVLbl'); 
            assertPositionExist(R,pos); 
            if iscell(pos)
                assert(numel(PLbl)==numel(pos),'if position list is a cell array, Lbl must be of same size')
            else
                pos={pos}; 
            end
            psname = R.PosNames; 
            ix = ismember(psname,pos);
            
            R.PIVlbl{ix}=PLbl;      
        end
        
        function PL = getPIVLbl(R,pos)
            assertPositionExist(R,pos);
            ix = ismember(R.PosNames,pos); 
            PL = R.PIVlbl{ix}; 
        end
        
        function R = merge(Rvec,varargin)
            % TODO: extend support for merging object with different
            % TimeVec structures, for now I'm assuming they are the same!
            arg.prefix='%g_'; 
            arg = parseVarargin(varargin,arg); 
            arg.Results = MultiPositionWoundResults; 
            R = merge@MultiPositionResults(Rvec,arg); 
            
            R.PIVlbl=PIVLbl.empty(1,0);
            for i=1:numel(Rvec)
                R.PIVlbl = [R.PIVlbl Rvec(i).PIVlbl]; 
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
      
        function HeatMapData(R, dataname,pos, frame, varargin)
                arg.clims = [-10, 10];
                arg = parseVarargin(varargin,arg); 

               a = R.getData(dataname,pos);
               a = a(:,frame);
               if min(a(:))>=0
                colormap(plasma())
               else
                colormap(makeColorMap([0.6 0 0.6],[1 1 1],[0.6 0.6 0]))
               end;
               imagesc(unique(R.PIVlbl{1}.X), unique(R.PIVlbl{1}.Y), reshape(a,numel(unique(R.PIVlbl{1}.X)),numel(unique(R.PIVlbl{1}.Y)))',arg.clims);
               set(gcf,'color','w');
               axis equal
               title(dataname)
               colorbar
               set(gca,'ydir','normal')
               shg;
        end
   end
   
   
    
end