function [meanposx,meanposy,stdposx,stdposy] = processExpData(indices, origins, directory, filename, outputname, fignum)

    time = {} ;        
    posx = {} ; 
    posy = {} ; 
    velx = {} ; 
    vely = {} ; 
    accx = {} ; 
    accy = {} ; 

    figure(fignum) 

    % import data
    for i = 1:length(indices)

        data = readtable( sprintf('%s/%s%0.0f_spotstats.csv',directory,filename,indices(i)) ) ; 
        data = table2array(data(:,3:end)) ; 
        ex = data(:,3) ; % [m]
        ey = data(:,4) ; % [m]

        % shift x and y over by coordinate
        ox = origins(i,1)  ; % [m]
        oy = origins(i,2)  ; % [m]
        posx{i} = ex - ox  ;
        posy{i} = oy - ey  ;
        et = data(:,6) ; 
        time{i} = et ; 

        % instantaneous velocity
        vx = diff(ex) ./ diff(et) ; % [m/s]
        vy = diff(ey) ./ diff(et) ; 

        % add zeros to pad velocity arrays
        vx = [0 ; vx] ; 
        vy = [0 ; vy] ; 
        velx{i} = vx ; 
        vely{i} = vy ;

        % acceleration
        ax = diff(vx) ./ diff(et) ; % [m/s2]
        ay = diff(vy) ./ diff(et) ; 

        % add zeros to pad acceleration arrays
        ax = [0 ; ax] ; 
        ay = [0 ; ay] ; 
        accx{i} = ax ; 
        accy{i} = ay ; 

        % plot position
    %     subplot(1,3,1)
        if i < length(indices)
            handleVis = 'off' ; 
        elseif i == length(indices) 
            handleVis = 'on' ; 
        end
        
        plot(posx{i},posy{i},'-k','LineWidth',1,'HandleVisibility',handleVis)
        hold on 
        title('Individual Swimmer Trajectories (n = 8)','interpreter','latex','FontSize',14); 
        xlabel('X [m]','interpreter','latex','FontSize',12) ; 
        ylabel('Y [m]','interpreter','latex','FontSize',12) ; 

    end

    try
        % calculate the mean for each data point in the cell arrays
        posxdata = cell2mat(posx) ; 
        posydata = cell2mat(posy) ; 

        % calculate averages across rows
        meanposx = mean(posxdata,2) ; 
        meanposy = mean(posydata,2) ; 

        stdposx = std(posxdata,0,2) ; 
        stdposy = std(posydata,0,2) ; 
    catch
        warning('There must be the same number of entries in every trajectory file. Using interpolation method instead. May introduce additional error.') ; 
        
        % interpolate data so all arrays are the same length
        % find the trial with the most data points
        A = cellfun(@length,time) ; % find the number of points in each cell in the time cell array
        [~,ind] = max(A) ; % ind is the index of the longest trial

        % get the x-position data points
        Xq = posx{ind} ; 
        Xq = Xq + 1E-4*rand([max(A),1]) ; % add small amount of randomization to Xq to allow interp1 to run on unique values

        % new y-position data
        Yinterp = zeros( max(A), length(indices) ) ; % as many rows as longest cell and as many columns as trials
        Yinterp(:,ind) = posy{ind} ; % the column for the longest trial is already known

        % use the x-position data as the query points for all other arrays
        for i = 1:length(indices) 
            if i == ind
                i = i + 1 ;  
            end

            posx{i} = posx{i} + 1E-4*rand([size(posx{i},1),1]) ; 
            posy{i} = posy{i} + 1E-4*rand([size(posy{i},1),1]) ; 

            % use interp1 to get new y-position data for all other arrays
            Yinterp(:,i) = interp1(posx{i},posy{i},Xq,'linear') ; 

        end

        % find average values
        meanposx = Xq ; 
        stdposx = std(meanposx,0,2) ; 
        meanposy = mean(Yinterp,2) ;
        stdposy = std(Yinterp,0,2) ; 
    end
    
    % save data to file
    headerA = ['meanposx (m) ', 'meanposy(m) ', 'stddevposx(m)','stddevposy(m) \n'] ; 
    
    A = [meanposx, meanposy, stdposx, stdposy]' ; 
    
    fid = fopen( sprintf('%s/%s.csv',directory,outputname), 'w' ) ; 
    fprintf( fid, headerA ) ; 
    fprintf(fid, '%.3f %.3f %.3f %.3f \n', A) ; 
    fclose(fid) ; 

    % add averaged trajectory with std dev to figure 
    figure(fignum)
    errorbar(meanposx,meanposy,stdposy,stdposy,'-r')
end