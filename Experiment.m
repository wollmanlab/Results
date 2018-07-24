classdef Experiment < handle
   % An Experiment object is a container of succesful Results object from
   % multiple Datasets. The goal is to be able to define an experiment 
   
   
    properties
                
        DatasetsNames
         
        projectpth
        Name
        
        Figures % a collection of Figure_ objects that are associated with this Experiment object
        Movies
        
    end
    
    properties (Dependent = true)
        User
        Project
        pth
    end
    
    methods
        
        function E = Experiment(projectpth,Name,varargin)
            E.Name = Name; 
            E.projectpth = projectpth; 
            if exist(fullfile(E.pth,'Experiment.mat'),'file')
                s = load(fullfile(E.pth,'Experiment.mat'));
                E = s.E;
            end
        end
        
        function saveExperiment(E)
            if ~exist(E.pth,'dir')
                mkdir(E.pth); 
            end
            save(fullfile(E.pth,'Experiment.mat'),E); 
        end
               
        function publish(E,varargin)
        end
        
        function addDataset(E,R)
            % what datasets are associated with this expriment
            E.DatasetsNames = unique([E.DatasetsNames; {R.pth}]); 
        end
        
        function addFigure(E,F)
            % add a figure to Experiment object 
             existingNames = cell(size(E.Figures)); 
            for i=1:numel(R.Figures)
                existingNames{i}=E.Figures(i).name; 
            end
            ix = ismember(existingNames,F.name); 
            if any(ix)
                E.Figures(ix)=F; 
            else
                E.Figures = [E.Figures; F];
            end
        end
                
        function showAllFigures(E)
            % DFigs is a 
        end
        
        function pth = get.pth(E)
            pth = fullfile(E.projectpth,'Experiments',E.Name);
        end
            
                
    end
    
    
end
    
    