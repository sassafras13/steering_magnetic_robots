function [xddot, yddot] = hypothesis(x,y,xdot,ydot,coils,mu0,mass,mag,volume,theta) 
    
    [gBx, gBy, gB] = magGradientCoilPoint(coils,mu0,x,y)      ; 
    gBangle = atan2(gBy,gBx) ;
    
    velAngle = atan2(ydot, xdot) ; 
    if (velAngle - gBangle) < pi/4 
        c = -1 ; 
    else
        c = 1 ; 
    end
    
    xddot = (1/mass) * ( (mag * volume * gB * cos(gBangle))...
        + (theta * c * sqrt(xdot^2 + ydot^2)) ) * cos(gBangle) ; 
    
    yddot = (1/mass) * ( (mag * volume * gB * cos(gBangle))...
        + (theta * c * sqrt(xdot^2 + ydot^2)) ) * sin(gBangle) ; 
end