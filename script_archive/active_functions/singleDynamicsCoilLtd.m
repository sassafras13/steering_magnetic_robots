function qdot = singleDynamicsCoilLtd(t,q,mag,mass,w,L,ct,cn,Barray,gBarray) 
    
%     disp('time') ; 
%     disp(t); 
    
    %%%%%%%%%%%%%%%%%
    % state variables
    %%%%%%%%%%%%%%%%% 
    x = q(1) ; 
    y = q(2) ; 
    xdot = q(3) ; 
    ydot = q(4) ; 
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    % magnetic gradient pull 
    %%%%%%%%%%%%%%%%%%%%%%%%
    
    % prepare data for interp2 function
    x_gB = gBarray(:,:,1) ; 
    y_gB = gBarray(:,:,2) ; 
    gBx = gBarray(:,:,3) ;  
    gBy = gBarray(:,:,4) ; 
    gB = gBarray(:,:,5) ; 
    
    % interpolate
    gBx_i = interp2(x_gB,y_gB,gBx,x,y) ; 
    gBy_i = interp2(x_gB,y_gB,gBy,x,y) ; 
    gB_i = interp2(x_gB,y_gB,gB,x,y) ; 
    
    gBangle = atan2(gBy_i,gBx_i) ;     
    Fbx = abs(mag*gB_i*cos(gBangle)) * cos(gBangle) ;
    Fby = abs(mag*gB_i*cos(gBangle)) * sin(gBangle) ; 
    
    %%%%%%%%%%%%%%%%%%%
    % hydrodynamic drag
    %%%%%%%%%%%%%%%%%%%
    % check if angle of motion matches angle of gradient to within epsilon;
    % this will tell us if the drag is negative (angles agree) or positive
    % (angles oppose)
    eps = pi/4 ; % [rad]
    velAngle = atan2(ydot, xdot) ; % [angle of velocity]
    if (velAngle - gBangle) < eps 
        c = -1 ; 
    else
        c = 1 ; 
    end
    
    % magnitude of drag force
    magFh = cn * sqrt( (xdot^2) + (ydot^2) ) ; 
    
    % drag forces in x and y
    Fhx = c * magFh * cos(gBangle) ;
    Fhy = c * magFh * sin(gBangle) ;    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % accelerations in x and y
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    xddot = (1/mass) * (Fbx + Fhx) ; 
    yddot = (1/mass) * (Fby + Fhy) ; 

    %%%%%%%%%%%%%%%%%%%%%%
    % return accelerations 
    %%%%%%%%%%%%%%%%%%%%%%
    qdot = [xdot ; ydot ; xddot ; yddot ] ;   

end