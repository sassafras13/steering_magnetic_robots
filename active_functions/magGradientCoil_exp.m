function [x_gB, y_gB, gBx, gBy, gB] = magGradientCoil_exp(X, Y, avgdata)

    % reshape X, Y and avgdata into grid shape
    X = reshape(X,[14,11]) ; 
    Y = reshape(Y,[14,11]) ; 
    avgdata = reshape(avgdata(:,1),[14,11]) ; 

    % avg data contains both mean and std dev, only need mean
    B = avgdata ; % [T]
    
    % the values of x and y for the meshgrid are the same as the input 
    x_gB = X ; % [m]
    y_gB = Y ; % [m]
    
    % these are just approximations I do not actually have data for this!
    [gBx, gBy] = gradient(B, X(1, :), Y(:, 1)) ; % [T]
    
    gB = sqrt( gBx.^2 + gBy.^2 ) ; % [T]
    
end