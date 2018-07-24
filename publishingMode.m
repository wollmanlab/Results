function publishingMode(state)

switch state
    case 'on'
        setappdata(0,'publishing',true)
    case 'off'
        setappdata(0,'publishing',false)
    otherwise
        warning('Unknown punlishing mode, please use on/off only')
end