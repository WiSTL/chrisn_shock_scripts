function [outputColormap] = customColormap(colorList,varargin)
        % creates a custom colormap
        % biafra ahanonu
        % started: 2014.01.03
        % inputs
                % colorList - a cell array containing vectors with RGB values, e.g. {[1 1 1], [0 0 1],[1 0 0]}
        % variable inputs
                % nPoints - numeric value specifying the number of points between each color value in the gradient
        % outputs
                %
 
        % changelog
                % 2014.05.09 - commented to make more readable
        % TODO
                %
 
        %========================
        options.nPoints = 50;
        % get options
        options = getOptions(options,varargin);
        % display(options)
        % unpack options into current workspace
        % fn=fieldnames(options);
        % for i=1:length(fn)
        %       eval([fn{i} '=options.' fn{i} ';']);
        % end
        %========================
 
        try
                % check that colorList was input, else return default map
                if isempty(colorList)
                        colorList = {[1 1 1], [0 0 1],[1 0 0]};
                end
                nColors = length(colorList);
                redMap = [];
                greenMap = [];
                blueMap = [];
                % loop over each color in the list and append its values to the RGB map values
                for i=1:(nColors-1)
                        % linspace is used to create an even gradient between the values specified
                        redMap = [redMap linspace(colorList{i}(1),colorList{i+1}(1),options.nPoints)];
                        greenMap = [greenMap linspace(colorList{i}(2),colorList{i+1}(2),options.nPoints)];
                        blueMap = [blueMap linspace(colorList{i}(3),colorList{i+1}(3),options.nPoints)];
                end
                % concatenate the RGB map values to create a custom colormap matrix
                outputColormap = [redMap', greenMap', blueMap'];
        catch err
                display(repmat('@',1,7))
                disp(getReport(err,'extended','hyperlinks','on'));
                display(repmat('@',1,7))
        end