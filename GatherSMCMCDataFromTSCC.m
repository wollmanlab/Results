%% This function executes the command to return an array of result objects
% This function receives the runpth and array of idArray
function [RtmpStructArray] = GatherSMCMCDataFromTSCC(runpth,idArray)

% initialize the RtmpArray
RtmpStructArray = cell(length(idArray),1);

% read the results object
for i = 1:length(idArray)
    localPth = fullfile(runpth,idArray{i});
    if exist(fullfile(localPth,'Results.mat')) > 0
        results = MultiRunComputationalResults(localPth);
        RtmpStructArray{i} = toStruct(results);    
    end 
end 
end

