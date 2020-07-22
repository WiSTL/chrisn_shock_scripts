function [options] = getOptions(options,inputArgs,varargin)
    % gets default options for a function, replaces with inputArgs inputs if they are present
    % biafra ahanonu
    % 2014.02.17 [22:21:49]
    %
    % inputs
    %   options - structure with options as fieldnames
    %   inputArgs - varargin containing name-value pairs passed from parent function.
 
    % list of valid options to accept, simple way to deal with illegal user input
    validOptions = fieldnames(options);
 
    % loop over each input name-value pair, check whether name is valid and overwrite fieldname in options structure.
    for i = 1:2:length(inputArgs)
        val = inputArgs{i};
        if ischar(val)
            % allow input of an options structure that overwrites existing fieldnames with its own, for increased flexibility
            if strcmp('options',val)
                inputOptions = inputArgs{i+1};
                [options] = mirrorRightStruct(inputOptions,options);
            elseif ~isempty(strmatch(val,validOptions))
                options.(val) = inputArgs{i+1};
            end
        else
            continue;
        end
    end
 
function [pullStruct] = mirrorRightStruct(pushStruct,pullStruct)
    % overwrites fields in pullStruct with those in pushStruct, other pullStruct fields rename intact
    % more generally, copies fields in pushStruct into pullStruct, if there is an overlap in field names, pushStruct overwrites.
    pushNames = fieldnames(pushStruct);
    for name = 1:length(pushNames)
        iName = pushNames{name};
        pullStruct.(iName) = pushStruct.(iName);
    end