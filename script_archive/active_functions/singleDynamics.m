function qdot = singleDynamicsCoil(t,q,mu0,gBmax,mag,m,z,I,a,nturns,w,L,h,mus)

    %%%%%%%%%%%%%%%%%
    % state variables
    %%%%%%%%%%%%%%%%%
    x = q(1) ; 
    y = q(2) ; 
%     theta = q(3) ; 
%     xdot = q(4) ; 
%     ydot = q(5) ; 
%     thetadot = q(6) ; 
    xdot = q(3) ; 
    ydot = q(4) ;
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    % magnetic gradient pull 
    %%%%%%%%%%%%%%%%%%%%%%%%
    [gBx, gBy, gB] = magGradientElectromagnetSinglePoint(mu0,gBmax,z,I,a,nturns,x,y)    ; 
    
    % interpolate values at a point
%     gBxq = interp2(gx,gy,gBx,x,y) ; 
%     gByq = interp2(gx,gy,gBy,x,y) ; 
%     gBq = interp2(gx,gy,gB,x,y) ; 
    
    gBangle = atan2(gBy,gBx) ; 
%     Fbangle = gBangle - theta ; 
%     Fb = abs(mag*gB*cos(Fbangle)) ;
    volume = w * L * h ; 
    Fb = abs(mag*volume*gB*cos(gBangle)) ; 
    
    %%%%%%%%%%%%%%%%%%%
    % hydrodynamic drag
    %%%%%%%%%%%%%%%%%%%
    % check if angle of motion matches angle of gradient to within epsilon;
    % this will tell us if the drag is negative (angles agree) or positive
    % (angles oppose)
    velAngle = atan2(ydot, xdot) ; 
    if (velAngle - gBangle) < pi/4 
        c = -1 ; 
    else
        c = 1 ; 
    end
    
    % net velocity
    vdot = sqrt((xdot^2) + (ydot^2)) ; 
    Fh = c*abs(6*pi*mus*vdot*w) ; 
%     Fh = abs(drag*vdot)  
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % accelerations in x and y
    %%%%%%%%%%%%%%%%%%%%%%%%%%
%     if Fb > Fh
%         vddot = (1/m)* (Fb - Fh) ; 
%     elseif Fh > Fb
%         vddot = 0 ; 
%     end
    vddot = (1/m) * (Fb + Fh) ; 
%     vddot = (1/m) * Fb ; 
    xddot = vddot*cos(gBangle) ; 
    yddot = vddot*sin(gBangle) ; 
    
%     %%%%%%%%%%%%%%%%%%%%%%%%%%
%     % hydrodynamic drag torque
%     %%%%%%%%%%%%%%%%%%%%%%%%%%
%     % tangential and normal vectors 
% %     that = [cos(theta) ; sin(theta)] ; 
%     nhat = [-sin(theta) ; cos(theta)] ; 
%     
%     % drag coefficients
%     ct = 6*pi*mus*w ; 
%     cn = 2*ct ; 
%     
%     % velocity of center of link
%     p0dot = [xdot ; ydot] ; 
%     
%     % torque 
%     Th = -( ((L^2)/4) * cn * p0dot' * nhat ) - ( ((L^3)/12) * cn * thetadot )  ; 
%         
%     %%%%%%%%%%%%%%%%%
%     % magnetic torque
%     %%%%%%%%%%%%%%%%%
%     % magnetic field
%     [gx2, gy2, ~, ~, B, ~] = magFieldCoil(coils,mu0,positionArray,Bmax) ; 
%     
%     % interpolate values at a point
%     Bq = interp2(gx2,gy2,B,x,y) ; 
%     
%     Tbangle = gBangle - theta ; 
%     Tb = mag*Bq*sin(Tbangle) ; 
%     
%     %%%%%%%%%%%%%%%%%%%%%%%
%     % acceleration in theta
%     %%%%%%%%%%%%%%%%%%%%%%%
%     % moment inertia of a rectangle rotating about z axis
%     I = (1/12) * m * ((w^2) + (L^2)) ; 
%     thetaddot = (1/I) * (Tb - Th) ;
    
    %%%%%%%%%%%%%%%%%%%%%%
    % return accelerations 
    %%%%%%%%%%%%%%%%%%%%%%
%     qdot = [xdot ; ydot ; thetadot ; xddot ; yddot ; thetaddot] ;     
    qdot = [xdot ; ydot ; xddot ; yddot] ;     

end