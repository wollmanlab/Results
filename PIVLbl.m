classdef PIVLbl < handle
    
    properties
        flist
        PosName
        pth
        WoundLbl
        Tvec
        ImageDims
    end
    
    properties (Transient = true)
        verbose = true;
    end
    
    properties (Dependent = true)
        X
        Y
    end
    
    
    methods
        
        function X = get.X(W)
            for i=1%:length(W.flist)
                fpthpos = [W.pth filesep W.PosName];
                data = load([fpthpos filesep W.flist{i}]);
                X = data(:,1);
            end
        end
        
        function Y = get.Y(W)
            for i=1%:length(W.flist)
                fpthpos = [W.pth filesep W.PosName];
                data = load([fpthpos filesep W.flist{i}]);
                Y = data(:,2);
            end
        end
        
        function U = U(W,frame)
            if nargin<2
                frame = 1:size(W.flist,2);
            end
            U = [];
            for i=frame%1%:length(W.flist)
                fpthpos = [W.pth filesep W.PosName];
                data = load([fpthpos filesep W.flist{i}]);
                U = [U ,data(:,3)];
            end
        end
        
        function UMasked = UMasked(W,frame)
            if nargin<2
                frame = 1:size(W.flist,2);
            end
            assert(strcmp(W.WoundLbl.PosName,W.PosName),'PIV and Wound must be from the same well!')
            for i=1:numel(frame)
                indInWound = W.WoundLbl.InWound(frame(i), [W.X W.Y]);
                UMasked(:,i) = W.U(frame(i)).*~indInWound;
            end
        end
        
        function UMaskedCorrected = UMaskedCorrected(W,frame)
            if nargin<2
                frame = 1:size(W.flist,2);
            end
            assert(strcmp(W.WoundLbl.PosName,W.PosName),'PIV and Wound must be from the same well!')
            for i=1:numel(frame)
                indInWound = W.WoundLbl.InWound(frame(i), [W.X W.Y]);
                UMaskedCorrected(:,i) = (W.U(frame(i))-mean(W.U(frame(i)))).*~indInWound;
            end
        end
        
        function V = V(W,frame)
            if nargin<2
                frame = 1:size(W.flist,2);
            end
            V = [];
            for i=frame%1%:length(W.flist)
                fpthpos = [W.pth filesep W.PosName];
                data = load([fpthpos filesep W.flist{i}]);
                V = [V ,data(:,4)];
            end
        end
        
        function VMasked = VMasked(W,frame)
            if nargin<2
                frame = 1:size(W.flist,2);
            end
            assert(strcmp(W.WoundLbl.PosName,W.PosName),'PIV and Wound must be from the same well!')
            for i=1:numel(frame)
                indInWound = W.WoundLbl.InWound(frame(i), [W.X W.Y]);
                VMasked(:,i) = W.V(frame(i)).*~indInWound;
            end
        end
        
        
        function VMaskedCorrected = VMaskedCorrected(W,frame)
            if nargin<2
                frame = 1:size(W.flist,2);
            end
            assert(strcmp(W.WoundLbl.PosName,W.PosName),'PIV and Wound must be from the same well!')
            for i=1:numel(frame)
                indInWound = W.WoundLbl.InWound(frame(i), [W.X W.Y]);
                VMaskedCorrected(:,i) = (W.V(frame(i))-mean(W.V(frame(i)))).*~indInWound;
            end
        end
        
        function indInWound = indInWound(W,frame)
            assert(strcmp(W.WoundLbl.PosName,W.PosName),'PIV and Wound must be from the same well!')
            for i=1:numel(frame)
                indInWound(:,i) = W.WoundLbl.InWound(frame(i), [W.X W.Y]);
            end
        end
        
        
        function setWoundLbl(W,WLbl,varargin)
            % do some checks
            assert(isa(WLbl,'WoundLbl'),'second variable WoundLbl must be of Class WoundLbl');
            assert(strcmp(WLbl.PosName,W.PosName),'PIV and Wound must be from the same well!')
            
            W.WoundLbl=WLbl;
        end
        
        function WL = getWoundLbl(W)
            WL = W.WoundLbl;
        end
        
        function plotPIVframe(W,frame, varargin)
            if nargin==2
                disp('Using automatic range')
            end
            X = W.X;
            Y = W.Y;
            U = W.UMasked(frame);
            V = W.VMasked(frame);
            
            arg.range=0.99*max(sqrt(U.^2+V.^2));
            arg = parseVarargin(varargin,arg);
            
            close all
            figure('position',[10, 10, 612/0.9, 512],'color','k')
            %axes('position',[0 0 0.9 1])
            ncquiverref_PIVMod(X,Y,U,V,[],[],[],'col',linspace(1,arg.range,20));
            set(gca,'color','k','visible','off','Ydir','reverse','xlim',[0 2448],'ylim',[0 2048],'position',[0 0 0.9 1]);
        end
        
    end
end
