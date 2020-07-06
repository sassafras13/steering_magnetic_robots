function theta = linRegX(expdata,coils,mu0,mass,mag,volume,alpha,thetaseed) 

    % input values
    x = expdata(:,1) ; 
    y = expdata(:,2) ; 
    xdot = expdata(:,3) ; 
    ydot = expdata(:,4) ;
    
    % experimental output values
    xddot = expdata(:,5) ; 
    yddot = expdata(:,6) ; 
    
    % convergence value 
    eps = 1E-5 ; 
    
    % initialize values
    thetaOld = 10*thetaseed ; 
    theta = thetaseed ; 
    count = 0 ; 
    
    % gradient descent
    % run until convergence on theta
    while abs(thetaOld - theta) > eps
        
        close all ; 
        
        count = count + 1 
        thetaOld = theta ; 
        
        % calculate cost function for a range of values on either side of theta
        thetaRange = (0.5*theta):(0.1*theta):(1.5*theta) ; 
        cost = zeros(length(thetaRange),1) ; 
        
        tic
        for i = 1:length(thetaRange) 
                        
            % array to contain hypothesis values
            h = zeros(length(x),1) ; % [xddot , yddot] 
            
            % calculate hypothesis for all inputs and current value of theta
            for j = 1:length(x) 
                [xddoti, ~] = hypothesis(x(j),y(j),xdot(j),ydot(j),...
                    coils,mu0,mass,mag,volume,thetaRange(i)) ; 
                h(j,1) = xddoti ; 
            end
        
            cost(i) = (1/(2 * length(x))) * sum((h - xddot) .^2)  
            
            figure(1)
            title('blue is actual data, other is hypothesis') ; 
            hold on
            plot(xddot,'ob') 
            hold on
            plot(h,'o')
            hold on
            
        end
        toc
        
        figure(2)
        xlabel('theta') ; ylabel('cost') ; 
        hold on
        plot(thetaRange,cost,'-') ; 
        hold on
       
        
        % calculate slope of cost function at current value of theta
        dJdthetaRange = diff(cost) ./ diff(thetaRange)' ; 
        [~, closestIndex] = min(abs(thetaRange' - theta)) ; 
        dJdtheta = dJdthetaRange(closestIndex) ; 
                
        % update value of theta
        theta = thetaOld - alpha * dJdtheta  
        
    end
        
end