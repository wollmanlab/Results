classdef MultiPositionSingleCellWoundResults < MultiPositionResults
    properties
        Frames = {};
        PIVlbl
        CorneaCellLbls
        Tracks
    end
    
    properties (Dependent = true)
        TimeVecs
    end
    
    methods (Static)
        function R = loadobj(S)
            R = MultiPositionSingleCellWoundResults;
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
            end
            if isfield(S,'CorneaCellLbls') % is the fields exists load them, if not,
                % they will be filled with default values
                % effectivly upcasting an object.
                %R.WoundLbl = S.WoundLbl;
                R.CorneaCellLbls = S.CorneaCellLbls;
            end
            if isfield(S,'Frames') % is the fields exists load them, if not,
                % they will be filled with default values
                % effectivly upcasting an object.
                R.Frames = S.Frames;
            end
            if isfield(S,'Tracks') % is the fields exists load them, if not,
                % they will be filled with default values
                % effectivly upcasting an object.
                R.Tracks = S.Tracks;
            end
        end
        
        function S = toStruct(R)
            % call the superclass method to start the transition to a
            % struct
            S = toStruct@MultiPositionResults(R);
            
            %S.WoundLbl = R.WoundLbl;
            S.PIVlbl = R.PIVlbl;
            S.CorneaCellLbls = R.CorneaCellLbls;
            S.Tracks = R.Tracks;
            S.Frames = R.Frames;
            %             S.TimeVecs = R.TimeVecs;
        end
        
        function R = MultiPositionSingleCellWoundResults(pth,reset) %constructor
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
        
        
        
        
        
        
        function setTracks(R,trcks,pos,varargin)
            if numel(R.Tracks)<R.Np
                R.Tracks{R.Np}=[];
            end
            
            arg.redo = false;
            arg = parseVarargin(varargin,arg);
            
            % do some checks
            %assert(isa(trcks,'struct'),'second variable trcks must be of Class struct');
            assertPositionExist(R,pos);
            if iscell(pos)
                if numel(pos)>1
                    error('One position at a time, please.')
                end
            else
                pos={pos};
            end
            psname = R.PosNames;
            ix = ismember(psname,pos);
            
            R.Tracks{ix}=trcks;
        end
        
        function PL = getTracks(R,pos)
            assertPositionExist(R,pos);
            ix = ismember(R.PosNames,pos);
            PL = R.Tracks{ix};
        end
        
        
        
        
        
        function setCorneaCellsLbl(R,CLbl,pos,varargin)
            if numel(R.CorneaCellLbls)<R.Np
                R.CorneaCellLbls{R.Np}=[];
            end
            
            arg.redo = false;
            arg = parseVarargin(varargin,arg);
            
            % do some checks
            assert(isa(CLbl{1},'CorneaCellsLbl'),'second variable PLbl must be of Class CorneaCellsLbl');
            assertPositionExist(R,pos);
            if iscell(pos)
                if numel(pos)>1
                    error('One position at a time, please.')
                end
            else
                pos={pos};
            end
            psname = R.PosNames;
            ix = ismember(psname,pos);
            
            R.CorneaCellLbls{ix}=CLbl;
        end
        
        function PL = getCorneaCellsLbl(R,pos)
            assertPositionExist(R,pos);
            ix = ismember(R.PosNames,pos);
            PL = R.CorneaCellLbls{ix};
        end
        
        
        
        
        
        
        
        function Link(R,pos)
            %% Now, we'll see how well we can track between adjacent frames (lap)
            
            %init assignment matrices
            Link12MatCell = {};
            Link21MatCell = {};
            
            searchRadius = 35;
            maxAmpRatio = 6;
            epsilon = 10^-4;%prevent cost==0
            CorneaCells = R.getCorneaCellsLbl(pos);
            inds = find(cellfun(@(x) ~isempty(x), CorneaCells));
            
            parfor i=inds(1:end-1)';
                i
                %build cost function
                Dists  = createDistanceMatrix(CorneaCells{i}.Centroids,CorneaCells{i+1}.Centroids);
                
                a1 = repmat(CorneaCells{i}.Intensities,1,numel(CorneaCells{i+1}.Intensities));
                a2 = repmat(CorneaCells{i+1}.Intensities',numel(CorneaCells{i}.Intensities),1);
                
                %divide the larger of the two amplitudes by the smaller value
                ampRatio = a1./a2;
                J = ampRatio < 1;
                ampRatio(J) = 1./ampRatio(J);
                %Keep only epithelium for tracking
                epiMatches = single(~isnan(CorneaCells{i}.epiScore))*single(~isnan(CorneaCells{i+1}.epiScore))';
                
                costMat = Dists.*(log2(ampRatio)+epsilon).*epiMatches;
                
                costMat(Dists>searchRadius) = 0;
                %costMat(ampRatio>maxAmpRatio) = 0;
                costMat(isnan(costMat)) = 0;
                costMat = sparse(double(costMat));
                %
                [Links12, Links21] = lap(costMat,[],[],1);
                %Make matrix of connections found
                Link12Mat =repmat(Links12(1:length(CorneaCells{i}.Centroids(:,1))),1,length(CorneaCells{i+1}.Centroids(:,1)));
                Link21Mat =repmat(Links21(1:length(CorneaCells{i+1}.Centroids(:,1))),1,length(CorneaCells{i}.Centroids(:,1)));
                
                LinkMat1 = meshgrid(1:length(CorneaCells{i+1}.Centroids(:,1)), 1:length(CorneaCells{i}.Centroids(:,1)));
                LinkMat2 = meshgrid(1:length(CorneaCells{i}.Centroids(:,1)), 1:length(CorneaCells{i+1}.Centroids(:,1)));
                
                Link12Mat = LinkMat1==Link12Mat;
                Link21Mat = LinkMat2==Link21Mat;
                
                Link12MatCell{i}=sparse(Link12Mat);
                Link21MatCell{i}=sparse(Link21Mat);
                %CorneaCells{i}.Link12Mat=sparse(Link12Mat);
                %CorneaCells{i}.Link21Mat=sparse(Link21Mat);
                
            end
            
            for i=inds(1:end-1)';
                CorneaCells{i}.Link12Mat = Link12MatCell{i};
                CorneaCells{i}.Link21Mat = Link21MatCell{i};
            end
            R.setCorneaCellsLbl(CorneaCells,pos)
        end
        
        function closeGaps(R,pos)
            CorneaCells = R.getCorneaCellsLbl(pos);
            inds = find(cellfun(@(x) ~isempty(x), CorneaCells));
            %% Book keeping: Make all tracks fragmants
            numFeatures = cellfun(@(x) x.num, CorneaCells(inds))';
            trackedFeatureIndx = (1:numFeatures(1))';
            numFrames = numel(inds);
            
            %initialize auxiliary matrices for storing information related to tracks
            %fragments
            numTracksWorstCase = round(sum(numFeatures)/10); %arbitrary large number
            trackedFeatureIndxAux = zeros(numTracksWorstCase,numFrames);
            rowEnd = numTracksWorstCase; %We'll fill this from the bottom up
            
            for i=1:numFrames-1
                %get indices of features in 2nd frame that are connected to features in 1st frame
                %indx1C - indexes in frame 1 that are linked to frame 2
                %indx2C - indexes in frame 2 that are linked to indx1C in frame 1
                numFeaturesFrame1 = numFeatures(i);
                numFeaturesFrame2 = numFeatures(i+1);
                [indx2C,indx1C] = find(CorneaCells{i}.Link21Mat);
                
                %%
                %find existing tracks that are not connected to features in 2nd frame
                numExistTracks = size(trackedFeatureIndx,1);
                indx1U = setdiff(1:numExistTracks,indx1C); %features in 1 not connected to 2
                numRows = length(indx1U);
                %%
                %determine where to store these tracks in auxiliary matrix
                %extend auxiliary matrices if necessary
                rowStart = rowEnd - numRows + 1;
                if rowStart <= 1
                    trackedFeatureIndxAux = [zeros(numTracksWorstCase,numFrames); ...
                        trackedFeatureIndxAux];
                    rowEnd = rowEnd + numTracksWorstCase;
                    rowStart = rowStart + numTracksWorstCase;
                end
                
                %% move rows of tracks that are not connected to points in
                %2nd frame to auxilary matrix
                trackedFeatureIndxAux(rowStart:rowEnd,1:i) = trackedFeatureIndx(indx1U,:);
                
                %%
                %assign space for new connectivity matrix
                tmp = zeros(numFeaturesFrame2,i+1);
                %fill in the feature numbers in 2nd frame
                tmp(1:numFeaturesFrame2,i+1) = (1:numFeaturesFrame2)';
                %shuffle existing tracks to get the correct connectivity with 2nd frame
                tmp(indx2C,1:i) = trackedFeatureIndx(indx1C,:);
                %update the connectivity matrix "trackedFeatureIndx"
                trackedFeatureIndx = tmp;
                
                %update rowEnd to indicate until which row the auxiliary
                %matrices are ampty
                rowEnd = rowStart - 1;
            end
            %
            %add information from last frame to auxiliary matrices
            numRows = size(trackedFeatureIndx,1);
            rowStart = rowEnd - numRows + 1;
            if rowStart <= 1
                trackedFeatureIndxAux = [zeros(numRows,numFrames); ...
                    trackedFeatureIndxAux];
                rowEnd = rowEnd + numRows;
                rowStart = rowStart + numRows;
            end
            trackedFeatureIndxAux(rowStart:rowEnd,:) = trackedFeatureIndx;
            
            %remove all empty rows
            trackedFeatureIndx = trackedFeatureIndxAux(rowStart:end,:);
            clear trackedFeatureIndxAux
            
            % get total number of tracks
            numTracks = size(trackedFeatureIndx,1);
            
            %find the frame where each track begins and then sort the vector
            frameStart = zeros(numTracks,1);
            for i=1:numTracks
                frameStart(i) = find((trackedFeatureIndx(i,:)~=0),1,'first');
            end
            [frameStart,indx] = sort(frameStart);
            
            %rearrange "trackedFeatureIndx" such that tracks are sorted in ascending order by their
            %starting point. Note that this ends up also arranging tracks starting at the
            %same time in descending order from longest to shortest.
            trackedFeatureIndx = trackedFeatureIndx(indx,:);
            
            
            %% Filter short tracks
            movieInfo = [];
            for i=inds'
                n = CorneaCells{i}.num;
                movieInfo(i).xCoord = [CorneaCells{i}.Centroids(:,1) zeros(n,1)];
                movieInfo(i).yCoord = [CorneaCells{i}.Centroids(:,2) zeros(n,1)];
                movieInfo(i).zCoord = [CorneaCells{i}.Centroids(:,3) zeros(n,1)];
                movieInfo(i).amp = [CorneaCells{i}.Intensities  zeros(n,1)];
                movieInfo(i).num = CorneaCells{i}.num;
            end
            probDim = 3;
            trackedFeatureInfo = coordAmpMatFromIndicesSparse(trackedFeatureIndx,movieInfo,...
                numFrames,probDim);
            trackSEL = getTrackSEL(trackedFeatureInfo);
            
            minTrackLen=3;
            
            %remove tracks whose length is less than minTrackLen
            indxKeep = find(trackSEL(:,3) >= minTrackLen);
            trackSEL = trackSEL(indxKeep,:);
            trackedFeatureIndx = trackedFeatureIndx(indxKeep,:);
            trackedFeatureInfo = trackedFeatureInfo(indxKeep,:);
            numTracks = size(trackSEL,1)
            clear movieInfo;
            
            
            
            
            %% Find all possible links based on thresholds
            maxTimeJump = 5;
            maxStep = 25;
            trackStartTime = trackSEL(:,1);
            trackEndTime   = trackSEL(:,2);
            %CentroidsStarts = [];
            %CentroidsEnds = [];
            
            %             for ind=1:numTracks
            %                 CentroidsStarts = [CentroidsStarts; full(trackedFeatureInfo(ind,8*(trackStartTime(ind)-1)+1:8*(trackStartTime(ind)-1)+3))];
            %                 CentroidsEnds = [CentroidsEnds; full(trackedFeatureInfo(ind,8*(trackStartTime(ind)-1)+1:8*(trackStartTime(ind)-1)+3))];
            %             end
            CentroidsStarts = cell2mat(arrayfun(@(a) full(trackedFeatureInfo(a,8*(trackStartTime(a)-1)+1:8*(trackStartTime(a)-1)+3))',1:numTracks,'UniformOutput',false))';
            CentroidsEnds = cell2mat(arrayfun(@(a) full(trackedFeatureInfo(a,8*(trackEndTime(a)-1)+1:8*(trackEndTime(a)-1)+3))',1:numTracks,'UniformOutput',false))';
            
            indx1=[];
            indx2=[];
            Dists=[];
            Skips = [];
            
            f = waitbar(0,'Calculating ditances of possible links...');
            for ind1 = 1:numTracks
                if ~rem(ind1,100)
                    waitbar(ind1/numTracks,f,'Calculating ditances of possible links...');
                end
                for ind2=1: numTracks
                    dT = trackStartTime(ind2)-trackEndTime(ind1);
                    if dT>0 && dT<=maxTimeJump;%condition on time
                        dR = norm(CentroidsStarts(ind1,:)-CentroidsStarts(ind2,:));
                        if dR<=sqrt(dT)*maxStep;%condition on space
                            indx1=[indx1 ind1];
                            indx2=[indx2 ind2];
                            Dists=[Dists dR];
                            Skips=[Skips dT];
                        end
                    end
                end
            end
            close(f)
            %% calculate cost matrix
            f = waitbar(0,'Gap closing...');
            
            [trackStats,statsRelChange,errFlag] = getTrackStats(trackedFeatureInfo,0.70,maxTimeJump);
            dispSqTheta = trackStats.dispSqTheta;
            dispSqR = trackStats.dispSqR;
            addConst =   - dispSqR.*log(dispSqTheta) + log(gamma(dispSqR));%log(ampDiffStd)
            
            costs = [];
            
            for i=1 : numel(indx1)
                if ~rem(i,100)
                    waitbar(i/numel(indx1),f,'Calculating cost matrix for gap closing...');
                end
                dta = Skips(i);
                distA = Dists(i);
                dispSqT = dispSqTheta(dta);
                dispSqAR = dispSqR(dta);
                addCons = addConst(dta);
                costs = [costs, dispSqT.*distA.^2-(dispSqAR-1).*log(max(distA.^2,realmin))+addCons];
            end
            close(f)
            costMat = sparse(indx1,indx2,costs,numTracks,numTracks);
            
            
            %% Close gaps with merging/splitting
            
            %if there are gaps to close (i.e. if there are tracks that start after the
            %first frame and tracks that end before the last frame) ...
            numTracksLink = size(trackedFeatureIndx,1);
            mergeSplit = 0;
            %if any(trackStartTime > 1) && any(trackEndTime < numFramesEff)
            
            
            
            %calculate the cost matrix, which already includes the
            %costs of birth and death
            % -- USER DEFINED FUNCTION -- %
            
            %link tracks based on this cost matrix, allowing for birth and death
            [link12,link21] = lap(costMat,[],[], 1);
            link12 = double(link12);
            link21 = double(link21);
            
            %put the indices of all tracks from linking in one vector
            tracks2Link = (1:numTracksLink)';
            tracksRemaining = tracks2Link;
            
            %reserve memory space for matrix showing track connectivity
            compoundTrack = zeros(numTracksLink,600);
            
            %initialize compTrackIndx
            compTrackIndx = 0;
            %
            while ~isempty(tracksRemaining)
                
                %update compound track index by 1
                compTrackIndx = compTrackIndx + 1;
                
                %take first track as a seed to build a compound track with
                %closed gaps and merges/splits
                trackSeed = tracksRemaining(1);
                seedLength = 1;
                seedLengthOld = 0; %dummy just to get into the while loop
                
                %while current seed contains more tracks than previous seed, i.e.
                %whie new track segments are still being added to the compound
                %track
                while seedLength > seedLengthOld
                    
                    %store current seed for later comparison
                    seedLengthOld = seedLength;
                    
                    %find tracks connected to ends of seed tracks
                    tmpTracks = link12(trackSeed);
                    trackLink2End = tmpTracks(tmpTracks <= numTracksLink); %starts linked to ends
                    trackMerge = [];
                    if mergeSplit
                        trackMerge = indxMerge(tmpTracks(tmpTracks > numTracksLink & ...
                            tmpTracks <= numTracksLink+numMerge) - numTracksLink); %tracks that ends merge with
                    end
                    
                    %find tracks connected to starts of seed tracks
                    tmpTracks = link21(trackSeed);
                    trackLink2Start = tmpTracks(tmpTracks <= numTracksLink); %ends linked to starts
                    trackSplit = [];
                    if mergeSplit
                        trackSplit = indxSplit(tmpTracks(tmpTracks > numTracksLink & ...
                            tmpTracks <= numTracksLink+numSplit) - numTracksLink); %tracks that starts split from
                    end
                    
                    %put all tracks together as the new seed
                    trackSeed = [trackSeed; trackLink2End; trackLink2Start; ...
                        trackMerge; trackSplit];
                    
                    %remove repetitions and arrange tracks in ascending order
                    trackSeed = unique(trackSeed);
                    
                    %get number of tracks in new seed
                    seedLength = length(trackSeed);
                    
                    %expand new seed if merging/splitting are allowed
                    if mergeSplit
                        
                        %variables storing merge/split seed tracks
                        mergeSeed = [];
                        splitSeed = [];
                        
                        %go over all seed tracks
                        for iSeed = 1 : seedLength
                            
                            %get the location(s) of this track in indxMerge
                            mergeSeed = [mergeSeed; find(indxMerge == trackSeed(iSeed))];
                            
                            %get the location(s) of this track in indxSplit
                            splitSeed = [splitSeed; find(indxSplit == trackSeed(iSeed))];
                            
                        end
                        
                        %add numTracksLink to mergeSeed and splitSeed to determine
                        %their location in the cost matrix
                        mergeSeed = mergeSeed + numTracksLink;
                        splitSeed = splitSeed + numTracksLink;
                        
                        %find tracks merging with seed tracks
                        trackMerge = [];
                        for iSeed = 1 : length(mergeSeed)
                            trackMerge = [trackMerge; find(link12(1:numTracksLink)==mergeSeed(iSeed))];
                        end
                        
                        %find tracks splitting from seed tracks
                        trackSplit = [];
                        for iSeed = 1 : length(splitSeed)
                            trackSplit = [trackSplit; find(link21(1:numTracksLink)==splitSeed(iSeed))];
                        end
                        
                        %add these track to the seed
                        trackSeed = [trackSeed; trackMerge; trackSplit];
                        
                        %remove repetitions and arrange tracks in ascending order
                        trackSeed = unique(trackSeed);
                        
                        %get number of tracks in new seed
                        seedLength = length(trackSeed);
                        
                    end %(if mergeSplit)
                    
                end %(while length(trackSeed) > length(trackSeedOld))
                
                %expand trackSeed to reserve memory for connetivity information
                trackSeedConnect = [trackSeed zeros(seedLength,2)];
                
                %store the tracks that the ends of the seed tracks are linked to,
                %and indicate whether it's an end-to-start link (+ve) or a merge (-ve)
                tmpTracks = link12(trackSeed);
                if mergeSplit
                    tmpTracks(tmpTracks > numTracksLink & tmpTracks <= ...
                        numTracksLink+numMerge) = -indxMerge(tmpTracks(tmpTracks > ...
                        numTracksLink & tmpTracks <= numTracksLink+numMerge) - numTracksLink);
                end
                tmpTracks(tmpTracks > numTracksLink) = NaN;
                trackSeedConnect(:,2) = tmpTracks;
                
                %store the tracks that the starts of the seed tracks are linked to,
                %and indicate whether it's a start-to-end link (+ve) or a split (-ve)
                tmpTracks = link21(trackSeed);
                if mergeSplit
                    tmpTracks(tmpTracks > numTracksLink & tmpTracks <= ...
                        numTracksLink+numSplit) = -indxSplit(tmpTracks(tmpTracks > ...
                        numTracksLink & tmpTracks <= numTracksLink+numSplit) - numTracksLink);
                end
                tmpTracks(tmpTracks > numTracksLink) = NaN;
                trackSeedConnect(:,3) = tmpTracks;
                
                %store tracks making up this compound track and their connectivity
                compoundTrack(compTrackIndx,1:3*seedLength) = reshape(...
                    trackSeedConnect,3*seedLength,1)';
                
                %in the list of all tracks, indicate that these tracks have
                %been taken care of by placing NaN instead of their number
                tracks2Link(trackSeed) = NaN;
                
                %retain only tracks that have not been linked to anything yet
                tracksRemaining = tracks2Link(~isnan(tracks2Link));
                
            end %(while ~isempty(tracksRemaining))
            
            %remove empty rows
            maxValue = max(compoundTrack,[],2);
            compoundTrack = compoundTrack(maxValue > 0,:);
            
            %determine number of tracks after gap closing (including merge/split if
            %specified)
            numTracksCG = size(compoundTrack,1);
            
            %reserve memory for structure storing tracks after gap closing
            tracksFinal = repmat(struct('tracksFeatIndxCG',[],...
                'tracksCoordAmpCG',[],'seqOfEvents',[]),numTracksCG,1);
            
            f = waitbar(0,'Closing gaps...');
            
            %go over all compound tracks
            for iTrack = 1 : numTracksCG
                if ~rem(iTrack,100)
                    waitbar(iTrack/numTracksCG,f,'Closing gaps...');
                end
                %get indices of tracks from linking making up current compound track
                %determine their number and connectivity
                trackSeedConnect = compoundTrack(iTrack,:)';
                trackSeedConnect = trackSeedConnect(trackSeedConnect ~= 0);
                seedLength = length(trackSeedConnect)/3; %number of segments making current track
                trackSeedConnect = reshape(trackSeedConnect,seedLength,3);
                
                %get their start times
                segmentStartTime = trackStartTime(trackSeedConnect(:,1));
                
                %arrange segments in ascending order of their start times
                [segmentStartTime,indxOrder] = sort(segmentStartTime);
                trackSeedConnect = trackSeedConnect(indxOrder,:);
                
                %get the segments' end times
                segmentEndTime = trackEndTime(trackSeedConnect(:,1));
                
                %calculate the segments' positions in the matrix of coordinates and
                %amplitudes
                segmentStartTime8 = 8 * (segmentStartTime - 1) + 1;
                segmentEndTime8   = 8 * segmentEndTime;
                
                %instead of having the connectivity in terms of the original track
                %indices, have it in terms of the indices of this subset of tracks
                %(which are arranged in ascending order of their start times)
                for iSeed = 1 : seedLength
                    value = trackSeedConnect(iSeed,2);
                    if value > 0
                        trackSeedConnect(iSeed,2) = find(trackSeedConnect(:,1) == ...
                            value);
                    elseif value < 0
                        trackSeedConnect(iSeed,2) = -find(trackSeedConnect(:,1) == ...
                            -value);
                    end
                    value = trackSeedConnect(iSeed,3);
                    if value > 0
                        trackSeedConnect(iSeed,3) = find(trackSeedConnect(:,1) == ...
                            value);
                    elseif value < 0
                        trackSeedConnect(iSeed,3) = -find(trackSeedConnect(:,1) == ...
                            -value);
                    end
                end
                
                %get track information from the matrices storing linking information
                tracksFeatIndxCG = trackedFeatureIndx(trackSeedConnect(:,1),:);
                tracksCoordAmpCG = trackedFeatureInfo(trackSeedConnect(:,1),:);
                
                %convert zeros to NaNs where approriate for the case of sparse
                %matrices
                if issparse(tracksCoordAmpCG)
                    
                    %convert sparse to full
                    tracksCoordAmpCG = full(tracksCoordAmpCG);
                    
                    %go over all the rows in this compound track
                    for iRow = 1 : size(tracksCoordAmpCG,1)
                        
                        %find all the zero entries
                        colZero = find(tracksCoordAmpCG(iRow,:)==0);
                        colZero = colZero(:)';
                        
                        %find the columns of the x-coordinates corresponding to
                        %the zero columns
                        xCoordCol = colZero - mod(colZero-1,8*ones(size(colZero)));
                        
                        %keep only the columns whose x-coordinate is zero as
                        %well
                        colZero = colZero(tracksCoordAmpCG(iRow,xCoordCol)==0);
                        
                        %replace zero with NaN in the surviving columns
                        tracksCoordAmpCG(iRow,colZero) = NaN;
                        
                    end
                    
                end
                
                %perform all gap closing links and modify connectivity accordingly
                %go over all starts in reverse order
                for iSeed = seedLength : -1 : 2
                    
                    %find the track this track might be connected to
                    track2Append = trackSeedConnect(iSeed,3);
                    
                    %if there is a track (which is not a split)
                    if track2Append > 0
                        
                        %put track information in the relevant row
                        tracksFeatIndxCG(track2Append,segmentStartTime(iSeed):...
                            segmentEndTime(iSeed)) = tracksFeatIndxCG(iSeed,...
                            segmentStartTime(iSeed):segmentEndTime(iSeed));
                        tracksFeatIndxCG(iSeed,:) = 0;
                        tracksCoordAmpCG(track2Append,segmentStartTime8(iSeed):...
                            segmentEndTime8(iSeed)) = tracksCoordAmpCG(iSeed,...
                            segmentStartTime8(iSeed):segmentEndTime8(iSeed));
                        tracksCoordAmpCG(iSeed,:) = NaN;
                        
                        %update segment information
                        segmentEndTime(track2Append) = segmentEndTime(iSeed);
                        segmentEndTime8(track2Append) = segmentEndTime8(iSeed);
                        segmentEndTime(iSeed) = NaN;
                        segmentEndTime8(iSeed) = NaN;
                        segmentStartTime(iSeed) = NaN;
                        segmentStartTime8(iSeed) = NaN;
                        
                        %update connectivity
                        trackSeedConnect(track2Append,2) = trackSeedConnect(iSeed,2);
                        trackSeedConnect(trackSeedConnect(:,2) == iSeed,2) = track2Append;
                        trackSeedConnect(trackSeedConnect(:,3) == iSeed,3) = track2Append;
                        trackSeedConnect(trackSeedConnect(:,2) == -iSeed,2) = -track2Append;
                        trackSeedConnect(trackSeedConnect(:,3) == -iSeed,3) = -track2Append;
                        
                    end %(if track2Append > 0)
                    
                end %(for iSeed = seedLength : -1 : 2)
                
                %find rows that are not empty
                maxValue = max(tracksFeatIndxCG,[],2);
                rowsNotEmpty = find(maxValue > 0);
                
                %remove empty rows
                tracksFeatIndxCG = tracksFeatIndxCG(rowsNotEmpty,:);
                tracksCoordAmpCG = tracksCoordAmpCG(rowsNotEmpty,:);
                segmentEndTime   = segmentEndTime(rowsNotEmpty);
                segmentStartTime = segmentStartTime(rowsNotEmpty);
                trackSeedConnect = trackSeedConnect(rowsNotEmpty,:);
                
                %update connectivity accordingly
                %by now, only merges and splits are left - thus no need for minus
                %sign to distinguish them from closed gaps
                for iSeed = 1 : length(rowsNotEmpty)
                    trackSeedConnect(trackSeedConnect(:,2) == -rowsNotEmpty(...
                        iSeed),2) = iSeed;
                    trackSeedConnect(trackSeedConnect(:,3) == -rowsNotEmpty(...
                        iSeed),3) = iSeed;
                end
                
                %determine new "seedLength"
                seedLength = length(rowsNotEmpty);
                
                %store the sequence of events of this track
                seqOfEvents = [segmentStartTime ones(seedLength,1) ...
                    (1:seedLength)' trackSeedConnect(:,3); ...
                    segmentEndTime 2*ones(seedLength,1) ...
                    (1:seedLength)' trackSeedConnect(:,2)];
                
                %sort sequence of events in ascending order of time
                [tmp,indxOrder] = sort(seqOfEvents(:,1));
                seqOfEvents = seqOfEvents(indxOrder,:);
                
                %add 1 to the times of merges
                indx = find(~isnan(seqOfEvents(:,4)) & seqOfEvents(:,2) == 2);
                seqOfEvents(indx,1) = seqOfEvents(indx,1) + 1;
                
                %find the frame where the compound track starts and the frames
                %where it ends
                frameStart = seqOfEvents(1,1);
                frameEnd   = seqOfEvents(end,1);
                
                %store final tracks, removing frames before anything happens and
                %after everything happens
                tracksFinal(iTrack).tracksFeatIndxCG = tracksFeatIndxCG(:,...
                    frameStart:frameEnd);
                tracksFinal(iTrack).tracksCoordAmpCG = tracksCoordAmpCG(:,...
                    8*(frameStart-1)+1:8*frameEnd);
                tracksFinal(iTrack).seqOfEvents = seqOfEvents;
                
            end %(for iTrack = 1 : numTracksCG)
            close(f)
            %%
            
            a = arrayfun(@(x) size(x.tracksFeatIndxCG,2),tracksFinal);
            longtracksFinal = tracksFinal(a>20);% remove all tracks shorter than 20 frames
            
            R.setTracks(longtracksFinal,pos);
        end
        
        
        
        function setEpitheliumFitParams(R,pos)
            CorneaCells = R.getCorneaCellsLbl(pos);
            
            for i=1:numel(CorneaCells)
                distScore = CorneaCells{i}.Centroids(:,3)-CorneaCells{i}.TopoZ;
                [h, Xbins] = histcounts(distScore,100,'Normalization', 'probability');
                Xbins = (Xbins(2:end)+Xbins(1:end-1))/2;
                if i==1
                    plot(Xbins,h);
                    shg
                    title('select gaussian area for lower (epithelium) peak')
                    pause;
                    J =InAxes;
                else
                    J = logical((Xbins>(BETA(2)-2*BETA(3))).*(Xbins<(BETA(2)+BETA(3))));
                end
                hToFit = h(J);
                XtoFit = Xbins(J);
                
                %
                posit = median(XtoFit)
                stdev = std(XtoFit)
                amp = max(hToFit)*sqrt(2*pi)*stdev
                BETA0 = [amp posit stdev];
                [BETA,RESNORM,RESIDUAL,EXITFLAG] = lsqcurvefit(@GaussianFit, BETA0 ,XtoFit, hToFit,[0 -inf 0], [inf inf inf]);
                x = min(Xbins):0.1:max(Xbins);
                if EXITFLAG>=0;
                    plot(Xbins, h, '-.', x, GaussianFit(BETA, x));
                    figure(gcf)
                end
                set(gca,'xlim',[-100,200],'ylim',[0,0.2]);
                drawnow;
                
                CorneaCells{i}.FitParam.func = @(x) GaussianFit(BETA,x);
                CorneaCells{i}.FitParam.cumFunc = @(x) (1+erf((x-BETA(2))/(sqrt(2)*BETA(3))))/2;
                CorneaCells{i}.FitParam.stdev = BETA(3);
                CorneaCells{i}.FitParam.posit = BETA(2);
                CorneaCells{i}.FitParam.amp = BETA(1);
                
                
                
                
                %Take upper and lower limits to keep as epithelium cells. +/-2\sigma
                lowLim = CorneaCells{i}.FitParam.posit-2*CorneaCells{i}.FitParam.stdev;
                highLim = CorneaCells{i}.FitParam.posit+2*CorneaCells{i}.FitParam.stdev;
                CorneaCells{i}.Jepi = ((distScore<highLim).*(distScore>lowLim));
                CorneaCells{i}.epiScore = 1-CorneaCells{i}.FitParam.cumFunc(distScore);
                CorneaCells{i}.epiScore(~CorneaCells{i}.Jepi)=NaN;
                
                CorneaCells{i}.Jepi = find(CorneaCells{i}.Jepi);
            end
        end
        
        function refineEpitheliumSurface(R,pos)
            %
            % Try to refine manifold...
            CorneaCells = R.getCorneaCellsLbl(pos);
            pth = R.pth;
            MD = Metadata(pth);
            imageSize = R.PIVlbl{1}.ImageDims;
            parfor i=1:numel(CorneaCells)
                i
                MD = Metadata(pth);
                Centroids = CorneaCells{i}.Centroids;
                Intensities = CorneaCells{i}.Intensities;
                CC = CorneaCells{i}.CC;
                Jepi = CorneaCells{i}.Jepi;
                
                %get function that defines the whole eputhelium
                F = scatteredInterpolant(Centroids(Jepi,1),Centroids(Jepi,2),Centroids(Jepi,3),'linear','nearest');
                
                
                
                Tforms = MD.getSpecificMetadata('driftTform','Position',pos, 'frame', i);
                dY = Tforms{1}(7);
                dX = Tforms{1}(8);
                
                pad = 100;
                x = (-pad+1:imageSize(2)+pad)+dY;
                y = (-pad+1:imageSize(1)+pad)+dX;
                
                [xx, yy] = meshgrid(x,y);
                zz = F(xx,yy);
                zz = imgaussfilt(zz,50);
                %
                TopoZ = interp2(xx,yy,zz,Centroids(:,1),Centroids(:,2),'nearest');
                %aha!
                DT = delaunayTriangulation([Centroids(:,1),Centroids(:,2),TopoZ]);
                [~, DistFromManifold] = nearestNeighbor(DT,Centroids(:,:));
                DistFromManifold = DistFromManifold.*sign(Centroids(:,3)-TopoZ);
                
                [~,J] = sort(DistFromManifold); %Sort by distance from smoothed epi map
                TopoZ = TopoZ(J);
                Centroids = Centroids(J,:);
                Intensities = Intensities(J);
                CC.PixelIdxList = CC.PixelIdxList(J);
                CorneaCells{i}.Centroids = Centroids;
                CorneaCells{i}.Intensities = Intensities;
                CorneaCells{i}.TopoZ = TopoZ;
                CorneaCells{i}.CC = CC;
            end
            R.setCorneaCellsLbl(CorneaCells,pos)
        end
        
        %
        
        
        
        function R = merge(Rvec,varargin)
            % TODO: extend support for merging object with different
            % TimeVec structures, for now I'm assuming they are the same!
            arg.prefix='%g_';
            arg = parseVarargin(varargin,arg);
            arg.Results = MultiPositionSingleCellWoundResults;
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
        
        
        
        function h = plotCellVsTime(R,pos)
            CorneaCells = R.getCorneaCellsLbl(pos);
            plot((0:(numel(R.CorneaCellLbls{1})-1))/2,cellfun(@(x) x.num, CorneaCells));
            ylabel('#Cells')
            xlabel('Time(h)')
            shg
        end
        
        
        
        
        %plotting
        function HeatMapData(R, dataname,pos, frame, varargin)
            arg.clims = [-15, 15];
            arg = parseVarargin(varargin,arg);
            
            a = R.getData(dataname,pos);
            a = a(:,frame);
                      if min(arg.clims)>=0
                          colormap(magma())
                      else
            colormap(makeColorMap([0.6 0 0.6],[0 0 0],[0.8 0.8 0]))
                      end;
            imagesc(unique(R.PIVlbl{1}.X), unique(R.PIVlbl{1}.Y), reshape(a,numel(unique(R.PIVlbl{1}.X)),numel(unique(R.PIVlbl{1}.Y)))',arg.clims);
            set(gcf,'color','w');
            axis equal
            title(dataname)
            colorbar
            set(gca,'ydir','normal')
            shg;
        end
        
        
        function h = PlotDisp(R,pos, i, dt,varargin)
            
            %function to plot displacement bw frame i and i+dt
            CorneaCells = R.getCorneaCellsLbl(pos);
            CentroidsD = CorneaCells{i}.Centroids;
            CentroidsA = CorneaCells{i+dt}.Centroids;
            
            
            CM = double(CorneaCells{i}.Link12Mat);%continuous connectivity map
            for j=i+1:i+dt-1
                CM=CM*double(CorneaCells{j}.Link12Mat);
            end
            [indx1, indx2] = find(CM);
            
            [indx1, J] = sort(indx1);
            
            idxToPlot=J;
            
            
            indx2 = indx2(idxToPlot);
            CentroidsD = CentroidsD(indx1,:);
            CentroidsA = CentroidsA(indx2,:);
            
            range = 1:length(indx1);
            if nargin>4
                range = varargin{1};
                if range(end)>indx1(end)
                    range = range(1):indx1(end);
                    disp('range end out of bounds, drawing all points up to the total # of cells.')
                end
                if range(1)>indx1(end)
                    range = 1:length(indx1);
                    disp('range fully out of bounds, drawing all points.')
                end
            end;
            
            %tzeva = viridis(length(indx1));
            %scatter3(CentroidsD(range,1),CentroidsD(range,2),-CentroidsD(range,3),[],tzeva(range,:),'*');
            %hold on
            %scatter3(CentroidsA(range,1),CentroidsA(range,2),-CentroidsA(range,3),[],tzeva(range,:));
            dx=CentroidsA(range,1)-CentroidsD(range,1);
            dy=CentroidsA(range,2)-CentroidsD(range,2);
            dz=CentroidsA(range,3)-CentroidsD(range,3);
            
            h = quiver3(CentroidsD(range,1),CentroidsD(range,2),-CentroidsD(range,3),dx,dy,-dz, 0);
            
            
            set(gca,'xlim',[-200 3000],'ylim',[-200 3000],'CameraPositionMode','manual','CameraPosition',[-2.6364e+03 -1.4045e+04 2.1889e+03])
            %hold on
            shg
        end
        
        
        
        %% Some sanity checks
        function allTracks = allTrackMatrix(R,pos)
            longtracksFinal = R.getTracks(pos);
            frames = R.Frames;
            allTracks = zeros(size(longtracksFinal,1),numel(frames));
            for i=1:size(longtracksFinal,1);
                trackStart = longtracksFinal(i).seqOfEvents(1,1);
                trackEnd = longtracksFinal(i).seqOfEvents(end,1);
                allTracks(i,trackStart:trackEnd) = longtracksFinal(i).tracksFeatIndxCG;
            end
        end
        
        function h = plotFracInTracks(R,pos)
            h = plot(sum(allTrackMatrix(R,pos)>0)'./cellfun(@(x) x.num, R.getCorneaCellsLbl(pos)));
            xlabel('Frame')
            ylabel({'Fraction of detected cells' 'in an active track'})
            shg
        end
        
        function histTrackLength(R,pos)
            hist(sum(allTrackMatrix(R,pos)>0,2),100);
            xlabel('Length (frames)')
            ylabel('#')
            shg
        end
        
        
        
        
    end
    
    
    
end