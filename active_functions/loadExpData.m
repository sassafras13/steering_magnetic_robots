function loadExpData(indices, origins, directory, filename, type)
    %loadExpData loads experimental data to the workspace as cell arrays
    
    if type == 1 
    
        time = {} ;        
        posx = {} ; 
        posy = {} ; 
        velx = {} ; 
        vely = {} ; 
        accx = {} ; 
        accy = {} ; 

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
            
            % save data to workspace
            assignin('base','time',time) ;
            assignin('base','posx',posx) ; 
            assignin('base','posy',posy) ; 
            assignin('base','velx',velx) ; 
            assignin('base','vely',vely) ; 
            assignin('base','accx',accx) ; 
            assignin('base','accy',accy) ;

        end
        
    elseif type == 0 
        data = readtable( sprintf('%s/%s.csv',directory,filename) ) ;
        data = table2array(data(1:end,1:4)) ; 
        
        meanposx = data(:,1) ; 
        meanposy = data(:,2) ; 
        stdposx = data(:,3) ; 
        stdposy = data(:,4) ; 
        assignin('base','meanposx',meanposx) ; 
        assignin('base','meanposy',meanposy) ; 
        assignin('base','stdposx',stdposx) ; 
        assignin('base','stdposy',stdposy) ; 
    end   
 

end