function [Lbl,NucLabels,CellLabels,CytoLabels,Timing] = segmentNucleiAndRingsWithoutTracking(MD,well,varargin)


%%all args have to be lower case
arg.nuc_channel = 'DeepBlue'; % nuc is usually Deepblue, using Hoechst
arg.nuc_smooth1 = 5; % sigma of filtering done to improve thresholding. The size of the nuclei will be eroded by this sigma as well
arg.nuc_smooth2 = fspecial('gauss',15,9);
arg.nuc_smooth = strel('disk',5); 
arg.nuc_minsize=0; %exclude small nuc
arg.nuc_maxsize=100000; % exclude big nuc, 
arg.cyto_channels = {'Cyan','Yellow','Red'}; 
arg.positiontype = 'Position'; 
arg.register = []; % optional registration object
arg.timefunc = @(t) true(size(t));
arg.cyto_ringstrel = strel('disk',15); 
arg.cyto_transform='log'; % other alternative: non, adjust, positive, no linear as Roy has indicated
arg.cyto_thresholdmethod = 'otsu'; % other alternatives: gm,minerr,robust,localotsu,kmeans
arg.shrinkmsk = strel('disk',50);
arg.ring_spacer = strel('disk',1); 
arg.ring_width = strel('disk',4);
arg.sz =  [2048  2064]; 
arg.track = 'none'; 
arg.verbose = false; 
arg.bkgperc = 5; %Background percentage, for step3, added by Yanfei
arg.subsamplefactor= 1; % subsample for gm method used for cyto segmenation, added by Yanfei, better not change it when use other methods
arg.cellerode = 0; %one more cell erode step for step 5
arg.posteriorthresh = 0.5; % depend on cell coverage, used in cyto mask calculation, added by Yanfei

arg = parseVarargin(varargin,arg); 

if ~iscell(arg.cyto_channels)%check whether the argument is cell type, usually not if there is only one channel
    arg.cyto_channels={arg.cyto_channels}; 
end

t0=now; 
t00=now; 

%step(1) get timepoints for the Label matrices
%step(2) register if requested or if Registration object was passed as an input arguemtn
%step(3) first subtrack background for entire stack
%step(4) for each well, segment all nuclei
%step(5) create a whole cell binary mask
%step(6) segment cell using nuc and finalize three labels matriceis
%step(7) create the Lbl and populate it

%%step(1) get timepoints for the Label matrices

%defaultTypes = {   'Z' 'Channel'    'Exposure'    'Fluorophore'    'Marker'    'PixelSize'    'PlateType'    'Position'    'Skip'    'TimestampFrame'    'TimestampImage'    'XY'    'acq'    'frame'    'group' 'AllInputs'}' 
Time = MD.getSpecificMetadata('TimestampFrame','Channel',arg.nuc_channel,arg.positiontype,well,'timefunc',arg.timefunc); %better to use Time here, T(Roy uses it) sometimes means types
Time = cat(1,Time{:});%convert cell to array

%the following three lines were original by Roy, change to the following two lines by Yanfei
%nuc = stkread(MD,'TimestampFrame',Time,'Channel',arg.nuc_channel); 
%[Time,ordr]= sort(Time);
%nuc=nuc(:,:,ordr); 

Time= sort(Time);
nuc = MD.stkread('TimestampFrame',Time,'Channel',arg.nuc_channel); %stkread has been overloaded in the Metadata class stkread(MD,prop1,value1)

Timing.readnucimages=now-t0;
t0=now; 

%why convert 1 to true and 0 to false, it looks unnecessary
%because there is a islogical check point
if arg.register==1
    arg.register=true; 
end
if arg.register==0
    arg.register=false; 
end

%%step(2) register if requested or if Registration object was passed as an input arguemtn
if islogical(arg.register) && arg.register
    
    %% find the frames where potential shifts occured (between acq)
    if numel(unique(MD,'acq'))==1
        possibleShiftingFrames=[]; 
    else
        tbl=MD.tabulate('acq','Position',well,'Channel',arg.nuc_channel,'timefunc',arg.timefunc);
        tbl=cat(1,tbl{:,2});
        possibleShiftingFrames=cumsum(tbl(1:end-1))+1;
    end
    [nuc,Tforms] = registerStack(nuc,'maxdisp',100,'onlyspecificframes',possibleShiftingFrames);
    Reg = Registration(Time,Tforms);
    arg.register = Reg; 
elseif ~isempty(arg.register) && isa(arg.register,'Registration')
    Reg = arg.register; 
    nuc = Reg.register(nuc,Time); 
else
    Reg = []; 
end

Timing.registration=now-t0;
arg.verbose && fprintf('Finished registration %s\n',datestr(now-t00,13)); 

%%step(3) first subtrack background for entire stack
% subtract bacgkround using a mask to avoid corner issues
nucprj = mean(nuc,3);
bkgperc = arg.bkgperc;
msk = nucprj>prctile(nucprj(:),bkgperc); % bkgperc was 5%, not tunable. Yanfei changed it to an tunable option
msk = imfill(msk,'holes');
msk  = bwareaopen(msk,10000);
if ~isempty(arg.shrinkmsk)
    msk = imerode(msk,arg.shrinkmsk); 
end
nuc = backgroundSubtraction(nuc,'msk',logical(msk),'smoothstk',false);


NucLabels = zeros([arg.sz size(nuc,3)],'uint16');

Timing.backgroundsubtraction=now-t0;
t0=now; 
arg.verbose && fprintf('Finished background subtraction %s\n',datestr(now-t00,13)); 

%% step(4) for each well, segment all nuclei
nuc_smooth1 = arg.nuc_smooth1; 
nuc_smooth2=arg.nuc_smooth2; 
nuc_maxsize = arg.nuc_maxsize; 
nuc_minsize = arg.nuc_minsize; 
parfor i=1:numel(Time)
    %% segment all nuclei
    nucfilt = imfilter(nuc(:,:,i),fspecial('gauss',2*nuc_smooth1,nuc_smooth1)) ;
    bw = optThreshold(nucfilt,'msk',msk,'method','localotsu','transform','none');
    pks = imregionalmax(imfilter(nucfilt,nuc_smooth2));
    bw  = bw &~bwareaopen(bw,nuc_maxsize); %discard large connected cells, default 1000, change it you your image is 1.5X
    bw = bwareaopen(bw,nuc_minsize); %discard small pieces, default is 0; 
    NucLabels(:,:,i) = segmentUsingSeeds(bw,bwlabel(pks & bw));
end

Timing.nucleisegmentation=now-t0;
t0=now; 
arg.verbose && fprintf('Finished nuclei segmentation %s\n',datestr(now-t00,13)); 

%% step(5) create a whole cell binary mask
CellLabels = zeros([arg.sz size(nuc,3)],'uint16');
CytoLabels = zeros([arg.sz size(nuc,3)],'uint16');
RingLabels = zeros([arg.sz size(nuc,3)],'uint16'); %these three lines moved from step 3 to here

cyto_thresholdmethod=arg.cyto_thresholdmethod;
cyto_transform=arg.cyto_transform; %there two line moved from inside of the following loop 
posteriorthresh = arg.posteriorthresh;
BW = false(size(nuc)); %was true , a bug
BWtmp=false(size(nuc)); %moved from inside the loop
subsample = numel(BW)/arg.subsamplefactor;

for j=1:numel(arg.cyto_channels)
    MAPK = stkread(MD,'TimestampFrame',Time,'Channel',arg.cyto_channels{j});
    if ~isempty(Reg)
        MAPK = Reg.register(MAPK,Time); 
    end
    MAPK = imfilter(MAPK,fspecial('gauss',5,3)); 

    parfor i=1:numel(Time)
        %% segment cytoplasm
        % next line does an AND and all cyto channels whereas the one
        % afterwards does an OR. 
%         BW(:,:,i) = BW(:,:,i) & optThreshold(MAPK(:,:,i),'msk',msk,'method',cyto_thresholdmethod,'transform',cyto_transform);
        BWtmp(:,:,i) = optThreshold(MAPK(:,:,i),'msk',msk,'method',cyto_thresholdmethod,'transform',cyto_transform, 'subsample', subsample, 'posteriorthresh',posteriorthresh);
    end
    BW=BW | BWtmp; % why use BWtmp here? because there might be two or more cyto channels
end
BW = imerode(BW, strel('disk',arg.cellerode)); %an extra step that roy does have it. Yanfei, cellerode default is 0;


Timing.cytothreshold=now-t0;
t0=now; 
arg.verbose && fprintf('Finished cyto thresholding %s\n',datestr(now-t00,13)); 

%%step(6) segment cell using nuc and finalize three labels matriceis
clear MAPK
cyto_ringstrel=arg.cyto_ringstrel; 
ring_spacer = arg.ring_spacer; 
ring_width=arg.ring_width; 

parfor i=1:numel(Time)
    
    %% limit cell signal to some positive distance from cyto
    nucmskdil =  imdilate(NucLabels(:,:,i)>0,cyto_ringstrel); 
    bw = BW(:,:,i) & nucmskdil;
        
    %% nuc must be where there is a cell signal, i.e. bw is true
    nuclbl = NucLabels(:,:,i); 
    nuclbl(~bw)=0; 
    nuclbl = bwlabel(nuclbl>0); 
    
    celllbl = segmentUsingSeeds(bw,nuclbl); 
    cytlbl = celllbl; 
    cytlbl(nuclbl>0)=0; 
    
      
    %% some bookeeping to make sure the numbers of all three match
    % get rid of nuclei without cytoplasm
    nuclbl(ismember(nuclbl,setdiff(nuclbl,cytlbl)))=0; 
    assert(isempty(setdiff(nuclbl,cytlbl)), 'differnt number of nuc and cyto') 
    nuclbl = bwlabel(nuclbl>0); 

    % redo segmentation but now with only "correct" nuclei, but this is an
    % inefficient way to do 
    celllbl = segmentUsingSeeds(bw,nuclbl); 
    cytlbl = celllbl; 
    cytlbl(nuclbl>0)=0; 
    
    % create a ring, overlap of rings is fixed. By Yanfei
    rnglbl=celllbl; 
    extnucbw = imdilate(nuclbl>0,ring_spacer);
    rngbw = imdilate(extnucbw,ring_width);
    rngbw = rngbw & ~extnucbw & cytlbl>0; 
    rnglbl(~rngbw) =0; % remove nuclei
           
     % get rid of nuclei without ring added by Yanfei
    nuclbl(ismember(nuclbl,setdiff(nuclbl,rnglbl)))=0; 
    nuclbl = bwlabel(nuclbl>0); 
    
    % redo segmentation but now with only "correct" nuclei, but this is an
    % inefficient way to do 
    celllbl = segmentUsingSeeds(bw,nuclbl); 
    cytlbl = celllbl; 
    cytlbl(nuclbl>0)=0; 
    
    % create a ring
    rnglbl=celllbl; 
    extnucbw = imdilate(nuclbl>0,ring_spacer);
    rngbw = imdilate(extnucbw,ring_width);
    rngbw = rngbw & ~extnucbw & cytlbl>0; 
    rnglbl(~rngbw) =0; % remove nuclei           
    
    % assign everything
    CellLabels(:,:,i)=celllbl;
    CytoLabels(:,:,i)=cytlbl; 
    NucLabels(:,:,i) = nuclbl; 
    RingLabels(:,:,i) = rnglbl; 
end

Timing.cytosegmentation=now-t0;
t0=now; 
arg.verbose && fprintf('Finished cyto segmentation %s\n',datestr(now-t00,13)); 

%%step(7) create the Lbl and populate it
Lbl = CellLabel;
Lbl.saveToFile=true; 
Lbl.pth=MD.pth; 
Lbl.posname = well; 
Lbl.Reg = Reg; 
addLbl(Lbl,CellLabels,'base',Time,'relabel',arg.track);
for i=1:numel(Time)
    addLbl(Lbl,CytoLabels(:,:,i),'cyto',Time(i),'relabel',arg.track);
    addLbl(Lbl,NucLabels(:,:,i),'nuc',Time(i),'relabel',arg.track);
    addLbl(Lbl,RingLabels(:,:,i),'ring',Time(i),'relabel',arg.track);
end

arg.verbose && fprintf('Finished creating cell label %s\n',datestr(now-t00,13));  %#ok<*VUNUS>
Timing.createlabelobject=now-t0;



