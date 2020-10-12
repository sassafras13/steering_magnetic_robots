function qdot = singleDynamicsCoil(t,q,mag,mass,w,L,ct,cn,Barray,gBarray) 
    
%     disp('time') ; 
%     disp(t); 
    
    %%%%%%%%%%%%%%%%%
    % state variables
    %%%%%%%%%%%%%%%%% 
    x = q(1) ; 
    y = q(2) ; 
    theta = q(3)  ;
    xdot = q(4) ; 
    ydot = q(5) ; 
    thetadot = q(6) ; 
    
    if or( (size(Barray,1) < 2), (size(gBarray,1) < 2) )
        qdot = zeros(6,1) ; 
    else
    
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
        Fbx = abs(mag*gB_i*cos(theta - gBangle)) * cos(gBangle) ;
        Fby = abs(mag*gB_i*cos(theta - gBangle)) * sin(gBangle) ; 

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

        % tangential and normal vectors 
        that = [cos(theta) ; sin(theta)] ; 
        nhat = [-sin(theta) ; cos(theta)] ; 

        % forces in tangential and normal directions
        Fht = ct * [xdot ; ydot]' * that ; 
        Fhn = cn * [xdot ; ydot]' * nhat ; 

        % magnitude of drag force
        magFh = sqrt( (Fht^2) + (Fhn^2) ) ; 

        % drag forces in x and y
        Fhx = c * magFh * cos(gBangle) ;
        Fhy = c * magFh * sin(gBangle) ;    

        %%%%%%%%%%%%%%%%%%%%%%%%%%
        % accelerations in x and y
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        xddot = (1/mass) * (Fbx + Fhx) ; 
        yddot = (1/mass) * (Fby + Fhy) ; 
    %     fprintf('Fbx = %.4f ; Fby = %.4f ; Fhx = %.4f ; Fhy = %.4f \n',Fbx,Fby,Fhx,Fhy); 

        %%%%%%%%%%%%%%%%%%%%%%%%%%
        % hydrodynamic drag torque
        %%%%%%%%%%%%%%%%%%%%%%%%%%

        % velocity of center of link
        p0dot = [xdot ; ydot] ; 

        % torque 
        Th = (-( ((L^2)/4) * cn * p0dot' * nhat ) - ( ((L^3)/12) * cn * thetadot ))   ;

        %%%%%%%%%%%%%%%%%
        % magnetic torque
        %%%%%%%%%%%%%%%%%

        % prepare data for interp2 function
        x_B = Barray(:,:,1) ; 
        y_B = Barray(:,:,2) ; 
        Bx = Barray(:,:,3) ;  
        By = Barray(:,:,4) ; 
        B = Barray(:,:,5) ; 

        % interpolate
        Bx_i = interp2(x_B,y_B,Bx,x,abs(y)) ; 
        By_i = interp2(x_B,y_B,By,x,abs(y)) ; 
        B_i = interp2(x_B,y_B,B,x,abs(y)) ; 

        Bangle = atan2(By_i,Bx_i) ;     
        Tb = mag * B_i * sin(Bangle - theta)  ;

        %%%%%%%%%%%%%%%%%%%%%%%
        % acceleration in theta
        %%%%%%%%%%%%%%%%%%%%%%%
        % moment inertia of a rectangle rotating about z axis
        I = (1/12) * mass * ((w^2) + (L^2)) ; 

        thetaddot = (1/I) * (Tb + Th) ;

        %%%%%%%%%%%%%%%%%%%%%%
        % return accelerations 
        %%%%%%%%%%%%%%%%%%%%%%
        qdot = [xdot ; ydot ; thetadot ; xddot ; yddot ; thetaddot] ;   
    end

end