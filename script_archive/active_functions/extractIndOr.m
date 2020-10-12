function [indices, origins] = extractIndOr(directory,filename) 
    % extractIndOr extracts the indices of the dataset to be analyzed and
    % the origin location in meters for each trial. To be used in analyzing
    % experimental data retrieved using ImageJ. The indexes matter because
    % we may not want to analyze all of the trials in a dataset. 
    %
    % Inputs:
    %   directory :- the directory location for the data to be analyzed.
    % 
    % Outputs: 
    %   indices :- the indexes of the data files to be analyzed. Corresponds
    %   to the trial number and video number. 
    %
    %   origins :- the x and y coordinates of the origin of the workspace as
    %   measured in ImageJ. These locations will be used to translate the
    %   swimmer trajectories.
    
    data = readtable( sprintf('%s/%s-indor.csv',directory,filename) ) ; 
    data = table2array( data(1:end,:) ) ; 
    
    indices = data(:,1) ; 
    origins = data(:,2:3) ; 
end