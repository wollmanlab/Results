function [Lbl,NucLabels,T,msk,Reg] = segmentNucleiOnly(MD,well,varargin)

%% define analysis parameters
% parameters for nuclei detection
arg.nuc_erode = strel('disk',3); % initial erosion to enhance nuclei centers
arg.nuc_smooth = fspecial('gauss',7,5); % filtering to smooth it out
arg.nuc_suppress = 0.05; % supression of small peaks - units are in a [0 1] space
arg.nuc_minarea = 30; % smaller then this its not a nuclei
arg.nuc_maxarea = Inf; 
arg.nuc_stretch = [1 99]; 
arg.nuc_channel = 'DeepBlue';
arg.mindistancefromedge =150;
arg.shrinkmsk = strel('disk',50);
arg.removetopprcentile = 0; 
arg.positiontype = 'Position'; 
arg.register = []; % optional registration object
arg.registerreference = []; 
arg.timefunc = @(t) true(size(t));
arg.project = false; 
arg.specificframeonly = []; % will duplicate that frame for all timepoints
arg.singleframe = false; % must be used in conjunction with specificframeonly. IF true will only return single timepoint
arg.track_method = 'none';
arg.threshold_method = 'otsu';

arg = parseVarargin(varargin,arg); 

if isempty(arg.project)
    arg.project=false; 
end

%% read the Hoescht stack
if iscell(arg.nuc_channel)
    nuc=cell(numel(arg.nuc_channel),1); 
    for j=1:numel(arg.nuc_channel)
        nuc{j} = stkread(MD,arg.positiontype,well,'Channel',arg.nuc_channel{j} ,'timefunc',arg.timefunc);
    end
    nuc=mean(cat(4,nuc{:}),4);
    arg.nuc_channel = arg.nuc_channel{1};
else
    nuc = stkread(MD,arg.positiontype,well,'Channel',arg.nuc_channel ,'timefunc',arg.timefunc);
end

%% Create the CellLabel object
Lbl = CellLabel;

%% get timepoints for the Label matrices
T = MD.getSpecificMetadata('TimestampFrame','Channel',arg.nuc_channel,arg.positiontype,well,'timefunc',arg.timefunc);
T = cat(1,T{:});
[T,ordr]= sort(T);
nuc=nuc(:,:,ordr); 

%% register if requested or if Registration object was passed as an input arguemtn
if islogical(arg.register) && arg.register
    if isempty(arg.registerreference)
        %% find the frames where potential shifts occured (between acq)
        if numel(unique(MD,'acq'))==1
            possibleShiftingFrames=[];
        else
            tbl=MD.tabulate('acq','Position',well,'Channel',arg.nuc_channel,'timefunc',arg.timefunc);
            tbl=cat(1,tbl{:,2});
            possibleShiftingFrames=cumsum(tbl(1:end-1))+1;
        end
    [nuc,Tforms] = registerStack(nuc,'maxdisp',100,'onlyspecificframes',possibleShiftingFrames);
    else
        [nuc,Tforms] = registerStack(nuc,'reference',arg.registerreference);
    end
    Reg = Registration(T,Tforms);
    arg.register = Reg; 
elseif ~isempty(arg.register) && isa(arg.register,'Registration')
    Reg = arg.register; 
    nuc = Reg.register(nuc,T); 
else
    Reg = []; 
end

%% projection
% define NucLabel before projection to maintain size
NucLabels = cell(size(nuc,3),1);
if arg.project  
    nuc = mean(nuc,3);
end

if ~isempty(arg.specificframeonly)
    nuc=nuc(:,:,arg.specificframeonly); 
end

if arg.singleframe && size(nuc,3)==1
    NucLabels(2:end)=[]; 
    T=T(arg.specificframeonly); 
end

%% first subtrack background for entire stack
% subtract bacgkround using a mask to avoid corner issues
nucprj = mean(nuc,3);
msk = nucprj>prctile(nucprj(:),5) & nucprj<=prctile(nucprj(:),100-arg.removetopprcentile);
  
msk = imfill(msk,'holes');
msk  = bwareaopen(msk,10000);
if ~isempty(arg.shrinkmsk)
    msk(1,:)=0; 
    msk(:,1)=0; 
    msk(end,:)=0; 
    msk(:,end)=0; 
    msk = imerode(msk,arg.shrinkmsk); 
end

bnd=bwboundaries(msk);
bnd=bnd{1}(:,[2 1]);

% nuc = backgroundSubtraction(nuc,'msk',logical(msk));


%% for each well, segment all nucleri
for i=1:size(nuc,3)
    %%
    %% segment nucleus
    % overall strategy -
    % perform a few operations to transform the nuclei into small hills
    % then find the local maxima to create a seed per each nucleur. Then
    % use watershed to segment the nuclei based on these seeds.
    nucbw = optThreshold(nuc(:,:,i),'method',arg.threshold_method,'msk',logical(msk),'transform','log');
    
    
    % identify peaks in the nuc

    % erode to emphasize the peaks
    nucpeaks = imerode(nuc(:,:,i),arg.nuc_erode);

    % smooths will get read of all the high frequency information and will
    % leabe a single peak per nuclei. We smooth with an input parametner
    % arg.nuc_smooth
    nucpeaks = imfilter(nucpeaks,arg.nuc_smooth);
    
    % next step is supression by intensity, to move to relative units, I am
    % rescaleing the image such that [0 1] are the 1% and 99% grayscale
    % values
    nucpeaks = mat2gray(nucpeaks,double(prctile(nucpeaks(nucbw),arg.nuc_stretch)));
    
    % imhmax supresses peaks below a specific height
    nucpeaks = imhmax(nucpeaks,arg.nuc_suppress);
    
    % at this point in the code each nuclei should have a unique sinlgle
    % peak, so we will use a regional max operation to find it. The &nucbw
    % is just to make sure that the identified peaks are within area of
    % hoecht and not some random peak in the background.
    nucpeaks = imregionalmax(nucpeaks) & nucbw;
    
    % next function is a custom segmentation function Roy wrote. It uses
    % combination of bwdistgeodesic and watershed to segment and then
    % removes objects that are too small and deals with the numbering
    nuclbl = segmentUsingSeeds(nucbw,bwlabel(nucpeaks),...
        'method','watershed','mincellarea',arg.nuc_minarea,'maxcellarea',arg.nuc_maxarea);
    
    % remove cells too close to boundaries
    prps = regionprops(nuclbl,'Centroid');
    cntr = cat(1,prps.Centroid)';
    DistFromEdge = min(distance(cntr,bnd'),[],2);
    nucbw = ismember(nuclbl,find(DistFromEdge>arg.mindistancefromedge));
    nucbw = bwareaopen(nucbw,arg.nuc_minarea);
    nuclbl = bwlabel(nucbw);
   
    %% add labels to cell array
    NucLabels{i} = nuclbl; 
    
end

if arg.project || ~isempty(arg.specificframeonly)
    addLbl(Lbl,repmat(NucLabels{1},[1 1 numel(NucLabels)]),'base',T,'relabel','none');
else
    %% add label to CellLabel object
    for i=1:numel(NucLabels)
        addLbl(Lbl,NucLabels{i} ,'base',T(i),'relabel',arg.track_method);
        addLbl(Lbl,NucLabels{i} ,'nuc',T(i),'relabel',arg.track_method);
    end
end

if isa(arg.register,'Registration')
    Lbl.Reg = Reg;
end
