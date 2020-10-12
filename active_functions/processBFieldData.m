function processBFieldData(directory,filenameGeneric,nfiles)
    % processBFieldData pulls in all data collected under a particular
    % filenameGeneric format and calculates the average and standard
    % deviation for every point in those files
    % 
    % Inputs: 
    %   directory :- directory name
    %
    %   filenameGeneric :- the prefix of the series of files to analyze
    %
    %   nfiles :- the number of files to analyze
    %
    % Outputs: 
    %   A single stats file output to the same experimental data file. 
    
    len = 14 * 11 ; % number of data points per file
    data = zeros(len,nfiles) ; % contains values from all files, each column is a different file
    
    for i = 1:nfiles
        temp = readtable( sprintf('%s/%s_r%.0f.csv',directory,filenameGeneric,i) ) ; 
        temp = table2array( temp(1:end,2:end) ) ; 
        data(:,i) = reshape(temp,[len,1]) ; 
    end
    
    avgB = mean(data,2) ; % take mean of every row
    stddevB = std(data,0,2) ;  % take std dev of every row, N-1 normalization
    
    % save data to file
    headerA = ['mean B (G) ', 'std dev B (G) \n'] ; 
    A = [avgB, stddevB]' ; 
    outputname = sprintf('%s_stats',filenameGeneric) ; 
    
    fid = fopen( sprintf('%s/%s.csv',directory,outputname), 'w' ) ; 
    fprintf( fid, headerA ) ; 
    fprintf(fid, '%.3f %.3f \n', A) ; 
    fclose(fid) ; 

end