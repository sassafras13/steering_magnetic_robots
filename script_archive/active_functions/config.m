function [coeff1, coeff2] = config(I1, I2)
    % config Calculates the coefficients of the current values for use in
    % magFieldCoil
    %
    % Inputs: 
    %   I1, I2 :- currents to coils 1 and 2
    %
    % Outputs: 
    %   coeff1, coeff2 :- binary values indicating current direction as +1
    %   or -1 depending on inputs
    
    if I1 > 0
        coeff1 = 1 ; 
    else
        coeff1 = -1 ; 
    end
    
    if I2 > 0
        coeff2 = 1 ; 
    else
        coeff2 = -1 ; 
    end
end