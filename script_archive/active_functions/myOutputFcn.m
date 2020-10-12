function status = myOutputFcn(t,y0,flag,mu0,gBmax,z,I,a,nturns,mag,w,L,h,mus)

    persistent Fb Fh

    switch flag
        case 'init'
            x = y0(1)  ;
            y = y0(2)  ;
            xdot = y0(3) ; 
            ydot = y0(4) ; 
            
            [gBx, gBy, gB] = magGradientElectromagnetSinglePoint(mu0,gBmax,z,I,a,nturns,x,y)    ; 
            gBangle = atan2(gBy,gBx) ; 
            volume = w * L * h ;
            Fb = abs(mag*volume*gB*cos(gBangle))  ; 
            
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
            
        case ''
            x = y0(1)  ;
            y = y0(2) ;
            xdot = y0(3) ; 
            ydot = y0(4) ; 
            
            [gBx, gBy, gB] = magGradientElectromagnetSinglePoint(mu0,gBmax,z,I,a,nturns,x,y)    ; 
            gBangle = atan2(gBy,gBx) ; 
            volume = w * L * h ;
            Fb = [Fb, abs(mag*volume*gB*cos(gBangle))]   ;
            
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
            Fh = [Fh, c*abs(6*pi*mus*vdot*w)] 
            
        case 'done'
            assignin('base','Fb',Fb) ; % save data to workspace
            assignin('base','Fh',Fh) ; 
    end
    
    status = 0 ; 
end