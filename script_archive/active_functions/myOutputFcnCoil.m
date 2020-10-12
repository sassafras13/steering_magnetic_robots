function status = myOutputFcnCoil(t,y0,flag,mag,mass,w,L,ct,cn,Barray,gBarray) 

    persistent Banglearray tarray gBanglearray Tb Th thetaddot
    
    switch flag
        case 'init'
            x = y0(1) ; 
            y = y0(2) ; 
            theta = y0(3) ; 
            xdot = y0(4) ; 
            ydot = y0(5) ; 
            thetadot = y0(6) ; 
            
            tarray = t(1) ; 
            
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
            gBanglearray = gBangle ; 
            

            Fbx = abs(mag*gB_i*cos(theta - gBangle)) * cos(gBangle) ;
            Fby = abs(mag*gB_i*cos(theta - gBangle)) * sin(gBangle) ; 

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

            % tangential and normal vectors 
            that = [cos(theta) ; sin(theta)] ; 
            nhat = [-sin(theta) ; cos(theta)] ; 

            Fht = ct * [xdot ; ydot]' * that ; 
            Fhn = cn * [xdot ; ydot]' * nhat ; 

            magFh = sqrt( (Fht^2) + (Fhn^2) ) ; 

            Fhx = c * magFh * cos(gBangle) ;
            Fhy = c * magFh * sin(gBangle) ;     

            %%%%%%%%%%%%%%%%%%%%%%%%%%
            % accelerations in x and y
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            xddot = (1/mass) * (Fbx + Fhx) ; 
            yddot = (1/mass) * (Fby + Fhy) ; 

            %%%%%%%%%%%%%%%%%%%%%%%%%%
            % hydrodynamic drag torque
            %%%%%%%%%%%%%%%%%%%%%%%%%%

            % velocity of center of link
            p0dot = [xdot ; ydot] ; 

            % torque 
            Th = (1E0) * (-( ((L^2)/4) * cn * p0dot' * nhat ) - ( ((L^3)/12) * cn * thetadot ) )  ;

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
            Banglearray = Bangle ; 

            Tb = mag * B_i * sin(Bangle - theta)  ;

            %%%%%%%%%%%%%%%%%%%%%%%
            % acceleration in theta
            %%%%%%%%%%%%%%%%%%%%%%%
            % moment inertia of a rectangle rotating about z axis
            I = (1/12) * mass * ((w^2) + (L^2)) ; 
            
            thetaddot = (1/I) * (Tb - Th) ;
         
%%            
        case ''
            x = y0(1) ; 
            y = y0(2) ; 
            theta = y0(3) ; 
            xdot = y0(4) ; 
            ydot = y0(5) ; 
            thetadot = y0(6) ; 
            
            tarray = [tarray, t(1)] ; 
            
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
            gBanglearray = [gBanglearray, gBangle] ; 
            
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

            % tangential and normal vectors 
            that = [cos(theta) ; sin(theta)] ; 
            nhat = [-sin(theta) ; cos(theta)] ; 

            Fht = ct * [xdot ; ydot]' * that ; 
            Fhn = cn * [xdot ; ydot]' * nhat ; 

            magFh = sqrt( (Fht^2) + (Fhn^2) ) ; 

            Fhx = c * magFh * cos(gBangle) ;
            Fhy = c * magFh * sin(gBangle) ;   

            %%%%%%%%%%%%%%%%%%%%%%%%%%
            % hydrodynamic drag torque
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            % velocity of center of link
            p0dot = [xdot ; ydot] ; 

            % torque 
            Th = [Th, (1E0) * (-( ((L^2)/4) * cn * p0dot' * nhat ) - ( ((L^3)/12) * cn * thetadot ))]  ; 

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
            Banglearray = [Banglearray, Bangle] ; 

            Tb = [Tb, mag * B_i * sin(Bangle - theta)] ;  
            
            % moment inertia of a rectangle rotating about z axis
            I = (1/12) * mass * ((w^2) + (L^2)) ; 
            
            thetaddot = [thetaddot, (1/I) * (Tb - Th)] ;
            
        case 'done'
            assignin('base','Bangle',Banglearray) ; % save data to workspace
            assignin('base','gBangle',gBanglearray) ; 
            assignin('base','tarray',tarray) ; 
            assignin('base','Tb',Tb) ; 
            assignin('base','Th',Th) ; 
            assignin('base','thetaddot',thetaddot) ; 
    end
    
    status = 0 ; 
end