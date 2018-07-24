function PIVlbl = PIVConstructor(pth, Well)
MD=Metadata(pth); 
NucFrame=stkread(MD,'Channel','DeepBlue','Position',Well,'specific',1,'flatfieldcorrection',false);

PIVlbl = PIVLbl;
PIVlbl.Tvec= MD.getSpecificMetadata('TimestampFrame','Position',Well,'Channel','DeepBlue');
PIVlbl.ImageDims = [size(NucFrame,1),size(NucFrame,2)];

fpthpos = [pth filesep Well];
    Temp = 'PIV_Results';
    endTemp = '.jvc';
    endTempNo = '.jvc~';
    flist = dir(fpthpos);
    flist = {flist.name};

    ix1 = regexp(flist, endTemp);ix2 = regexp(flist, Temp);ix3 = regexp(flist, endTempNo);
    ix=~cellfun('isempty',ix1) & ~cellfun('isempty',ix2) & cellfun('isempty',ix3) ;
    flistDo = flist(ix);
    PIVlbl.flist = flistDo;
    PIVlbl.pth = pth;
    PIVlbl.PosName = Well;
    PIVlbl.setWoundLbl(WoundTracker_v3(pth, Well),Well)
end
