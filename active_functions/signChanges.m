function [gBx,gBy] = signChanges(gBx,gBy,gx,gy,I1,I2)
    % signChanges Change the signs of the gradient components gBx and gBy
    % to reflect where in the workspace the gradients were calculated. This
    % function uses information about the quadrant that the observation
    % point is located in, as well as the direction of current flowing
    % through each of the coils. Negative current is pointing N to the 
    % left. Positive current is pointing N to the right. 
    %
    % Inputs: 
    %   gBx, gBy :- the x and y components of the gradient at all points in
    %   the workspace
    %
    %   gx, gy :- the x and y meshgrid values that specify the locations of
    %   all the calculation points of the gradient
    %
    %   I1, I2 :- the current through each coil. I1 is left coil, I2 is 
    %   right coil.
    %
    % Outputs: 
    %   gBx, gBy :- the x and y components of the gradient at all points in
    %   the workspace, with their signs adjusted for gradient direction. 

    for i = 1:(size(gBx,1)*size(gBx,2))

        % switch: quadrant (quadrant is a function that returns quadrant as [1,4]
        % based on x and y coordinates of current point) 
        quadrant = checkQuad(gx(i),gy(i)) ; 
        switch quadrant
            % case 1: if I2 is positive, then flip the signs on both x and y components
            % if I2 is negative then do nothing
            case 1
                if (I2 > 0)
                    gBx(i) = -gBx(i) ; 
                    gBy(i) = -gBy(i) ; 
                end
            % case 2: if I1 is positive, then flip the signs on both x and y components
            % if I2 is negative then do nothing    
            case 2
                if (I1 > 0) 
                    gBx(i) = -gBx(i) ; 
                    gBy(i) = -gBy(i) ; 
                end
            % case 3: if I1 is positive, then do nothing
            % if I1 is negative then flip the signs on both x and y components
            case 3
                if (I1 < 0)
                    gBx(i) = -gBx(i) ; 
                    gBy(i) = -gBy(i) ; 
                end
            % case 4: if I2 is positive then do nothing
            % if I2 is negative then flip the signs on both x and y components
            case 4
                if (I2 < 0) 
                    gBx(i) = -gBx(i) ; 
                    gBy(i) = -gBy(i) ; 
                end
        end
    end     
end