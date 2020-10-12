function [X,Y,avgdata] = plotBFieldData(directory,filenameGeneric,nfiles,rawOn,nfig)
    % plotBFieldData plots all the Bfield data for one configuration on a
    % 3D plot
    
    len = 14 * 11 ; % number of data points per file
    data = zeros(len,nfiles) ; % contains values from all files, each column is a different file
    
    % setup x and y coordinates
    xgv = -5:1:5 ;  
    ygv = -6:1:7 ; 
    [X,Y] = meshgrid(xgv,ygv) ; 
    
    X = reshape(X,[len,1]) ; 
    X = X * 1E-2 ;  % [m]
    Y = reshape(Y,[len,1]) ; 
    Y = Y * 1E-2 ;  % [m]
    
    if rawOn == 1 
        for i = 1:nfiles
            temp = readtable( sprintf('%s/%s_r%.0f.csv',directory,filenameGeneric,i) ) ; 
            temp = table2array( temp(1:end,2:end) ) ; 
            data(:,i) = reshape(temp,[len,1]) ; 
        end
                
        figure(1)
        for i = 1:nfiles
            plot3(X,Y,data(:,i),'o') ; 
            hold on
        end
    end
    
    avgdata = readtable( sprintf('%s/%s_stats.csv',directory,filenameGeneric) ) ;
    avgdata = table2array( avgdata(1:end,1:2) ) ;  
    avgdata = avgdata * 1E-4 ; % convert from G to T
        
    posErr = 0.005*ones(len,1) ; % position error is 0.005 m
    
    % formatting
    markerS = 10 ; 
    lineW = 3 ; 
    
    figure(nfig)
    plot3(X,Y,avgdata(:,1),'xr','MarkerSize',markerS) ; 
    hold on
    plot3([X,X]',[Y,Y]',[-avgdata(:,2),avgdata(:,2)]'+avgdata(:,1)','-r','LineWidth',lineW,'HandleVisibility','off') ; 
    hold on
    plot3([-posErr,posErr]'+X',[Y,Y]',[avgdata(:,1),avgdata(:,1)]','-r','LineWidth',lineW,'HandleVisibility','off') ; 
    hold on
    plot3([X,X]',[-posErr,posErr]'+Y',[avgdata(:,1),avgdata(:,1)]','-r','LineWidth',lineW,'HandleVisibility','off') ; 
    hold on
end