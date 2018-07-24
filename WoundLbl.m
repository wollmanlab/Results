classdef WoundLbl < handle
    
    properties 
        PolyXY
        IsThereAWound
        PosName
        pth
        ImageDims
        %Tvec
    end
    
    properties (Transient = true)
       verbose = true; 
    end
    
    properties (Dependent = true) 
       %Mask
       %Centroid
       %Area
       filename
    end
    

    
    methods



%         function filename = get.filename(WoundLbl)
%             filename = fullfile(WoundLbl.pth,sprintf('WoundLabel_%s.mat',WoundLbl.pos));
%         end

        function Mask = Mask(W,frame)
                %Mask = ones(W.ImageDims(1),W.ImageDims(2),size(W.PolyXY,3));
                Mask = ~poly2mask(double(W.PolyXY{frame}(:,1)),double(W.PolyXY{frame}(:,2)),W.ImageDims(1),W.ImageDims(2));
        end
        
        function Centroid = Centroid(W,frame)
            if nargin<2
                frame =1:numel(W.PolyXY);
            end
            Centroid = zeros(numel(frame),2);
            for i=1:numel(frame)
                if W.IsThereAWound(frame(i))
                    Mask = ~W.Mask(frame(i));
                    A = regionprops(Mask,'Centroid');
                    try
                        Centroid(i,:) = A.Centroid;
                    catch
                        Centroid(i,:) = Centroid(i-1,:);
                    end
                elseif W.IsThereAWound(1)
                    lastWframe = find(diff(W.IsThereAWound));
                    Centroid(i,:) = Centroid(max(lastWframe-5,1));
                else
                    Centroid(i,:) = ceil(W.ImageDims./2);
                    %warning('no wound, returning center of frame.')
                end
            end
        end
        
        function Cent = SmoothCentroid(W,frame)
            if nargin<2
                frame =1:numel(W.PolyXY);
            end
            flt = GaussianFit([1, 0, 1], -5:5);
            Cent = zeros(numel(frame),2);
            PadBefore = -(frame(1)-5-1)*((frame(1)-5)<1);
            PadAfter =(frame(end)+5-numel(W.IsThereAWound))*((frame(end)+5)>numel(W.IsThereAWound));
            CentRough = W.Centroid(max(frame(1)-5,1):min(frame(end)+5,length(W.IsThereAWound)));
            
            CentroidEnv = [repmat(CentRough(1,:),PadBefore,1);...
            CentRough;...
            repmat(CentRough(end,:),PadAfter,1)];

            for i=1:numel(frame)
                    Cent(i,:) = sum(repmat(flt',1,2).*CentroidEnv(i:i+10,:));
            end
        end
        
        function Area = Area(W,frame)
            if nargin<2
                frame =1:numel(W.PolyXY);
            end
            Area = zeros(numel(frame),1);
            %for frame = 1:size(W.PolyXY,3);
            for i=1:numel(frame)
                if W.IsThereAWound(frame(i))
                    Mask = ~W.Mask(frame(i));
                    A = regionprops(Mask,'Area');
                    try
                        Area(i) = A.Area;
                    catch
                        Area(i)=0;
                    end;
                else
                    Area(i) = 0;
                    %warning('no wound, Area is 0.')
                end
            end
            %end
        end

        function ix = InWound(W, frame,XY)
            Wpoly = W.PolyXY{frame};
            if numel(unique(Wpoly,'rows'))>2
                ix = inpolygon(XY(:,1), XY(:,2), Wpoly(:,1), Wpoly(:,2));
            else
                ix = zeros(size(XY,1),1);
            end;
        end

    end
    
end
        