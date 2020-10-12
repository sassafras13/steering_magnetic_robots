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
    eps = 1E-3 ; 
    
    % initialize values
    thetaOld = 10*thetaseed ; 
    theta = thetaseed ; 
    
    % gradient descent
    % run until convergence on theta
    while abs(thetaOld - theta) > eps
        
        thetaOld = theta ; 
        
        % calculate cost function for a range of values on either side of theta
        thetaRange = (0.5*theta):(0.1*theta):(1.5*theta) ; 
        costx = zeros(length(thetaRange),1) ; 
        costy = zeros(length(thetaRange),1) ; 
        
        for i = 1:length(thetaRange) 
                        
            % array to contain hypothesis values
            h = zeros(length(x),2) ; % [xddot , yddot] 
            
            % calculate hypothesis for all inputs and current value of theta
            for j = 1:length(x) 
                [xddoti, yddoti] = hypothesis(x(j),y(j),xdot(j),ydot(j),...
                    coils,mu0,mass,mag,volume,thetaRange(i)) ; 
                h(j,1) = xddoti ; 
                h(j,2) = yddoti ; 
            end
        
            costx(i) = (1/(2 * length(x))) * sum((h(:,1) - xddot) .^2) ; 
            costy(i) = (1/(2 * length(x))) * sum((h(:,2) - yddot) .^2) ; 
        end

        % calculate slope of cost function at current value of theta
        dJdthetaRangex = diff(costx) ./ diff(thetaRange) ; 
        dJdtheta = dJdthetaRange(find(thetaRange == theta)) ; 
                
        % update value of theta
        theta = thetaOld - alpha * dJdtheta ; 
        
    end
        
end