function [x_B, y_B, Bx, By, B] = magFieldCoil_exp(X, Y, avgdata)

    % reshape X, Y and avgdata into grid shape
    X = reshape(X,[14,11]) ; 
    Y = reshape(Y,[14,11]) ; 
    avgdata = reshape(avgdata(:,1),[14,11]) ; 

    % avg data contains both mean and std dev, only need mean
    B = avgdata ; % [T]
    
    % the values of x and y for the meshgrid are the same as the input 
    x_B = X ; % [m]
    y_B = Y ; % [m]
    
    % these are just approximations I do not actually have data for this!
    Bx = B ; % [T]
    By = -1*B ; % [T]
    
end