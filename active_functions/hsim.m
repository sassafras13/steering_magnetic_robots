function [x,y] = hsim(niter,waypointT,y0,mag,mass,w,L,ct,Barray,gBarray)
    cn = 2*ct ; 
    yarray = {} ; 
    timearray = {} ; 

    for i = 1:niter

        % time 
        T0 = waypointT(i) ; 
        Tf = waypointT(i+1) ; 
        Tspan = [T0 Tf] ; 

        options = odeset('RelTol',1E-3,'AbsTol',1E-6,'Stats','on','OutputFcn',...
            @(t,y,flag) myOutputFcnCoil(t,y0,flag,mag,mass,w,L,ct,cn,Barray,gBarray)) ;
        
        tic 
        [t,y] = ode15s(@(Tspan, y0) singleDynamicsCoil(Tspan,y0,mag,...
            mass,w,L,ct,cn,Barray,gBarray), Tspan, y0, options) ; 
        toc

        yarray{i} = y ; 
        timearray{i} = t ; 

        y0 = [y(end,1) ; y(end,2) ; y(end,3) ; 0 ; 0 ; 0] ;

    end

    % then query the results using the experimental data as query points (with
    % some random variation to prevent duplicate query points)

    % convert all cells to one long matrix
    ydata = [] ; 

    for i = 1:niter
        ydata = [ydata ; yarray{i}(:,:)] ; 
    end

    x = ydata(:,1) ; % just pick out the x and y data
    y = ydata(:,2) ; 

%     yq = interp1(x,y,xq) ; 
%     
%     % replace NaN with 0
%     yq(isnan(yq)) = 0 ; 
%     
% 
%     % h = [xq , yq] ; % hypothesis data is x and y data for the full time range
%     h = yq ;

end