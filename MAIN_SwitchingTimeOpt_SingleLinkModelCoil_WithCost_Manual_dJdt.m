% MAIN_SwitchingTimeOpt_SingleLinkModelCoil.m
% Purpose: Implement switching time optimization for macroswimmer gradient
% steering project. Includes a terminal cost to the cost function and
% applies weights to running and terminal costs. 
% Author: Emma Benjaminson
% References: 
%
% [1] Johnson, E and Murphey, T. "Second-Order Switching Time Optimization
% for Nonlinear Time-Varying Dynamic Systems." IEEE Trans on Auto. Control,
% vol 56(8), August 2011. 
%% 

clc ; clear all ; close all ; 

%% Directories

localDir = fileparts(mfilename('fullpath')) ;
restoredefaultpath ;
addpath(fullfile(localDir, 'active_functions')) ;

%% Setup Variables

% load macroswimmer data
run VARIABLES.m ; 
load('MAT_gBarray.mat') ; 
load('MAT_Barray.mat') ; 

% define perturbation amount epsilon 
epsi = [1E-4 ; 1E-4 ; 1E-6 ; 1E-5] ; % perturbation of state variables 
alpha = 5E3 ; % step size
alpha_eps = 1E2 ; % amount to increase alpha

gamma = 0.90 ; % update parameter for refining alpha
gammaBIG = 0.75 ; % update parameter for when the current tau values resulted in a bad trajectory

% define initial guess at tau values
tau = [0,80,110,180,350] ; % [s]
tau_eps = 1 ; % [s] perturbation value for tau
tau_per = tau ; 
tau_per(2:end) = tau_per(2:end) + tau_eps ; 
tau_old = 0 ; 

% provide control inputs u_i 
u = [1.9, 0 ; 
     0 , 1.9 ; 
     1.9, 0 ; 
     0 , 0 ]; % [A] currents to coils
    
% start location
y0start = [-0.024 , 0.04 , pi/2 , 0 , 0 , 0]' ; 
    
% goal location 
xd = [-0.01 ; 0.03 ; 0 ; 0] ; 

% create desired trajectories
nwaypoints = 2 ; 
xwaypoints = [y0start(1:2)' ; ones(nwaypoints,2) ; xd(1:2)'] ; 

% set the waypoints (i.e. points where trajectory changes direction)
for i = 1:nwaypoints
    if rem(i,2) == 1 % if i is odd
        xwaypoints(i+1,1) = xd(1) ; 
        xwaypoints(i+1,2) = y0start(2) - ( i * (y0start(2) - xd(2)) / (nwaypoints + 1) ) ; 
    elseif rem(i,2) == 0 
        xwaypoints(i+1,1) = y0start(1) ; 
        xwaypoints(i+1,2) = y0start(2) - ( i * (y0start(2) - xd(2)) / (nwaypoints + 1) ) ; 
    end
end
    
npoints = 100 ; 
xd_t = [] ; 
td = [] ; 

% fill in points along the trajectory between waypoints
for i = 1:(nwaypoints+1)
    xd_x = linspace(xwaypoints(i,1),xwaypoints(i+1,1),npoints) ; 
    xd_y = linspace(xwaypoints(i,2),xwaypoints(i+1,2),npoints) ; 
    xd_t = [xd_t ; [xd_x', xd_y']] ; 
    
    td_i = linspace(tau(i),tau(i+1),npoints) ; 
    td = [td , td_i] ; 
end

xd_t = [xd_t , zeros(size(xd_t,1),2) ] ; 

% plot desired trajectory to check
figure(1)
title('desired trajectory') ; 
plot(xd(1),xd(2),'og','MarkerSize',5,'MarkerFaceColor','g') ; 
hold on 
plot(xd_t(:,1),xd_t(:,2),'-b') ; 
hold on

% define weights for cost function
w1 = [1,1,0,0] ; % weight for running cost
w2 = [1,1,0,0] ; % weight for final cost

% initialize derivatives of cost function
dJdtau = zeros((size(tau,2)-1),1) ; 

% for saving old versions of the trajectory
ynom_old = 0 ; 
tnom_old = 0 ; 

%% MAIN

% count iterations
iter = 1 ; 
niter = 1001 ; 

% save all the norms dJdtau
normArray = zeros(niter,1) ; 

% continue to update tau_i until convergence condition is met: 
% while norm([dJdtau]) > 1E-4 
while iter < niter
    
    % forward simulate the dynamics for current values of tau
    % OutputFcn saves time array and qdot 
    
    tic
    
    %%%%% NORMAL %%%%%
    disp('Calculate normal') ;
    ynom = [] ; 
    tnom = [] ; 
    
    % start location
    y0 = y0start ; 
    
    for i = 1:3
        lastwarn('') ;
        tspan = [tau(i) tau(i+1)] ; 
        options = odeset('Stats','off','OutputFcn',...
            @(t,y,flag) OutputFcn(t,y0,flag,tau,u,0,i,mag,mass,w,L,ct,cn,Barray,gBarray)) ; 
        [tnom_i,ynom_i] = ode15s(@(t,y) dynFunc(t,y,tau,u,mag,mass,w,L,ct,cn,Barray,gBarray), tspan, y0, options) ; 
        
        [warnMsg, warnId] = lastwarn ; 
        if ~isempty(warnMsg)
            ynom = ynom_old ; 
            tnom = tnom_old ; 
            tau = tau_old ; 
            alpha = alpha * gammaBIG ; 
            y0 = [ynom(end,1) ; ynom(end,2) ; ynom(end,3) ; 0 ; 0 ; 0 ] ;
            i = 4 ; 
            break ; 
        else
            ynom = [ynom ; ynom_i] ;  
            tnom = [tnom ; tnom_i] ;  
            y0 = [ynom_i(end,1) ; ynom_i(end,2) ; ynom_i(end,3) ; 0 ; 0 ; 0 ] ;
        end
         
    end
    toc
    
    if isempty(ynom) == 1 
        ynom = ynom_old ; 
        alpha = alpha * gammaBIG ; 
    end
        
    if iter == 1 
        width = 3 ; 
        color = [0.5 0.2 0.6] ; 
    else
        width = 2 ; 
        color = [(0.5 - 0.05*(iter-1)) 0.2 0.6] ;
    end
    
    if or( iter == 1, (rem(iter,10) == 0) ) 
        figure(1)
        title('simulated trajectories') ; 
        plot(ynom(:,1),ynom(:,2),'-','LineWidth',width,'Color',color) ; % trajectory
        hold on
        plot(ynom(1,1),ynom(1,2),'og') ; % start point
        hold on 
        plot(ynom(end,1),ynom(end,2),'o','MarkerSize',5,'MarkerFaceColor',color) ; % end point
        hold on
    end
    
    % update before truncating
    ynom_old = ynom ; 
    tnom_old = tnom ; 
    
    % truncate ynom to 4 state variables
    ynom = ynom(:,[1,2,4,5]) ; 
    ydot_0 = [ydot_0_1([1,2,4,5],:)' ;
              ydot_0_2([1,2,4,5],:)' ; 
              ydot_0_3([1,2,4,5],:)' ]; 
%               ydot_0_4([1,2,4,5],:)' ] ; 
          
    % concatenate time array
    timeArray_0 = [timeArray_0_1' ; 
                   timeArray_0_2' ; 
                   timeArray_0_3' ]; 
%                    timeArray_0_4' ] ; 

               
    % calculate perturbed values 
    tic
    
    %%%%% PERTURBED %%%%%
    disp('Calculate perturbed values') ;
    yper = [] ; 
    tper = [] ; 
    
    % start location
    y0 = y0start ; 
    
    for i = 1:3
        lastwarn('') ;
        tspan = [tau_per(i) tau_per(i+1)] ; 
        options = odeset('Stats','off','OutputFcn',...
            @(t,y,flag) OutputFcn(t,y0,flag,tau_per,u,1,i,mag,mass,w,L,ct,cn,Barray,gBarray)) ; 
        [tper_i,yper_i] = ode15s(@(t,y) dynFunc(t,y,tau_per,u,mag,mass,w,L,ct,cn,Barray,gBarray), tspan, y0, options) ; 
        
        [warnMsg, warnId] = lastwarn ; 
        if ~isempty(warnMsg)
            yper = yper_old ; 
            tper = tper_old ; 
            tau_per = tau_per_old ; 
            alpha = alpha * gammaBIG ; 
            y0 = [yper(end,1) ; yper(end,2) ; yper(end,3) ; 0 ; 0 ; 0 ] ;
            i = 4 ; 
            break ; 
        else
            yper = [yper ; yper_i] ;  
            tper = [tper ; tper_i] ;  
            y0 = [yper_i(end,1) ; yper_i(end,2) ; yper_i(end,3) ; 0 ; 0 ; 0 ] ;
        end
         
    end
    toc
    
    if isempty(yper) == 1 
        yper = yper_old ; 
        alpha = alpha * gammaBIG ; 
    end
    
    % truncate ynom to 4 state variables
    yper = yper(:,[1,2,4,5]) ; 
    y_per_dot_1 = [ydot_1_1([1,2,4,5],:)' ;
              ydot_1_2([1,2,4,5],:)' ; 
              ydot_1_3([1,2,4,5],:)' ]; 
%               ydot_1_4([1,2,4,5],:)' ] ; 
          
    % concatenate time array
    timeArray_per_1 = [timeArray_1_1' ; 
                   timeArray_1_2' ; 
                   timeArray_1_3' ]; 
%                    timeArray_1_4' ] ; 

    
    %%%%%%%
    % calculate Xi for each transition
    Xi = zeros((size(tau,2)-1),size(ynom,2)) ; 
    
    for i = 2:size(tau,2)
        tau_i = tau(i) ; 
        Xi(i-1,:) = XiFunc(timeArray_0,ydot_0,tau_i) ; 
    end

    % find dJdtau by backwards integrating psi from tf to tau_i
    tspan_psi = [tau(end-1), tau(2)] ; 
    
    % psi0 (initial condition for psi) is equal to derivative of the
    % terminal cost m(x,tf)
%     psi0 = dmdx(ynom,xd,w2,epsi) ; 
    psi0 = zeros(4,1) ; 
    
    % backwards integrate psi_dot
    options = odeset('RelTol',1E0,'AbsTol',1E0) ; 
    [t_psi,psi] = ode15s(@(t_psi,psi) psiFunc(t_psi,psi,tnom,ynom,td,...
        xd_t,timeArray_0,ydot_0,tper,yper,timeArray_per_1,y_per_dot_1,w1,...
        mag,mass,w,L,ct,cn,Barray,gBarray,tau,epsi), tspan_psi, psi0, options) ; 

    % interpolate psi values at tau time points
    psi_array = zeros((size(tau,2)-1),size(psi0,1)) ; 
    
    % calculate dJdtau and psi_array
    for i = 1:(size(tau,2)-1)
        psi_array(i,:) = interp1(t_psi,psi,tau(i+1)) ; 
        psi_array(isnan(psi_array)) = 0 ; 
        dJdtau(i,1) = psi_array(i,:) * (Xi(i,:)') ; 
    end

    % perform steepest descent to update tau_i
    dJdtau
    tau 
    H = eye(size(tau,2)-1) ; 
    tau_old = tau ; 
    
    % choose an effective value of alpha
    % calculate the current value of J
    tspan = [tau(1) tau(end)] ; 
    J0 = 0 ; 
    [t,J] = ode15s(@(t,J) Jfunc(t,J,tnom,ynom,td,xd_t,w1), tspan, J0) ;
    J = sum(J) ; 
    
    Jnew = 1E9 ; 
    i_alpha = 1 ; 
    
    % loop through until we find a good value of alpha
    while Jnew > J
    
        % calculate new tau values
        tau_new = zeros( 1, (size(dJdtau,1)+1) ) ; 
        tau_new(2:end) = tau(2:end) - ((alpha*inv(H)*dJdtau)')  

        % forward simulate the new dynamics with these new tau values

        %%%%% NORMAL WITH NEW TAU VALUES %%%%%
        disp('Calculate normal with new tau values') ;
        ynom_new = [] ; 
        tnom_new = [] ; 

        % start location
        y0 = y0start ; 

        for i = 1:3
            tspan = [tau_new(i) tau_new(i+1)] ; 
            options = odeset('Stats','off','OutputFcn',...
                @(t,y,flag) OutputFcn(t,y0,flag,tau,u,0,i,mag,mass,w,L,ct,cn,Barray,gBarray)) ; 
            [tnom_i,ynom_i] = ode15s(@(t,y) dynFunc(t,y,tau,u,mag,mass,w,L,ct,cn,Barray,gBarray), tspan, y0, options) ; 
            ynom_new = [ynom_new ; ynom_i] ;  
            tnom_new = [tnom_new ; tnom_i] ;  
            y0 = [ynom_i(end,1) ; ynom_i(end,2) ; ynom_i(end,3) ; 0 ; 0 ; 0 ] ;
        end
        
        % truncate to 4 state variables
        ynom_new = ynom_new(:,[1,2,4,5]) ;

        % calculate the new value of J
        tspan = [tau_new(1) tau_new(end-1)] ; 
        J0 = 0 ; 
        [t,Jnew] = ode15s(@(t,J) Jfunc(t,J,tnom_new,ynom_new,td,xd_t,w1), tspan, J0) ;
        Jnew = sum(Jnew)  
        
        % if new J is greater than the current value of J, decrease alpha by
        % alpha ^i
        if Jnew > J
%             alpha = alpha^(-i_alpha) ; 
            alpha = alpha * gamma ; 
        % if new J is less than the current value of J, increase alpha by alpha
        % + epsilon
        elseif Jnew <= J
            alpha = alpha + alpha_eps ; 
        end
       
        i_alpha = i_alpha + 1  
    end
    
    
    tau(2:end) = tau(2:end) - ((alpha*inv(H)*dJdtau)') 
    normArray(iter,1) = norm(dJdtau) ; 
    
    if or( iter == 1, (rem(iter,10) == 0) ) 
        figure(2)
        title('dJdtau') ;
        hold on
        for i = 1:size(dJdtau,1)
            plot(iter,dJdtau(i),'-ok') 
        end
        hold on

        figure(3)
        title('tau')  ;
        hold on
        for i = 1:size(tau,2) 
            plot(iter,tau(i),'-ok') 
        end
        hold on 
    end
    
    iter = iter + 1 
    
end
    
norm(dJdtau)

tau

figure(4)
plot(normArray)
xlabel('Iteration') ; ylabel('Norm of dJdtau') ; 

%% Function Definitions

function ydot = dynFunc(t,y,tau,u,mag,mass,w,L,ct,cn,Barray,gBarray)

    if size(tau,2) > 1 

        % find out which mode we are in based on current time and apply
        % appropriate dynamics 
        if and( (t >= tau(1)), (t <= tau(2)) ) % there must be a better way to do this than hard-coding in the indices of tau?
            u1 = u(1,1) ; 
            u2 = u(1,2) ; 
            Barray_i = Barray{1} ; 
            gBarray_i = gBarray{1} ; 
        elseif and( (t > tau(2)), (t <= tau(3)) )
            u1 = u(2,1) ; 
            u2 = u(2,2) ; 
            Barray_i = Barray{2} ; 
            gBarray_i = gBarray{2} ; 
        elseif and( (t > tau(3)), (t <= tau(4)) )
            u1 = u(3,1) ; 
            u2 = u(3,2) ; 
            Barray_i = Barray{1} ; 
            gBarray_i = gBarray{1} ; 
        elseif and( (t > tau(4)), (t <= tau(5)) )
            u1 = u(4,1) ; 
            u2 = u(4,2) ; 
            Barray_i = 0 ; 
            gBarray_i = 0 ; 
        end
        
    elseif size(tau,2) == 1
        u1 = u(1,1) ; 
        u2 = u(1,2) ; 
    end
        
    ydot = singleDynamicsCoil(t,y,mag,mass,w,L,ct,cn,Barray_i,gBarray_i) ; 
    
end

function status = OutputFcn(t,y0,flag,tau,u,i,ii,mag,mass,w,L,ct,cn,Barray,gBarray)

    persistent timeArray ydot
    
    switch flag
        case 'init' 
            ydot = dynFunc(t(1),y0,tau,u,mag,mass,w,L,ct,cn,Barray,gBarray) ; 
            timeArray = t(1) ;
        case ''
            ydot = [ydot, dynFunc(t(1),y0,tau,u,mag,mass,w,L,ct,cn,Barray,gBarray)] ; 
            timeArray = [timeArray, t(1)] ; 
        case 'done'
            assignin('base',sprintf('ydot_%.0f_%.0f',i,ii),ydot) ; % save data to workspace
            assignin('base',sprintf('timeArray_%.0f_%.0f',i,ii),timeArray) ; 
    end

    status = 0 ; 
end

function Xi = XiFunc(t,ydot,tau_i)
    [~,k] = min( abs(t-tau_i) )  ; % find the index of the time value closest to t_psi
    if k < 3
        xdot_i = ydot(k,:) ; % state at time tau
    else
        xdot_i = ydot( (k-2),: ) ; 
    end
    
    ki = k + 2; 
    if ki > size(ydot,1)
        ki = k ; 
    end
    
    xdot_i1 = ydot(ki,:) ; 
    Xi = xdot_i - xdot_i1 ; 
    
end

    
function psidot = psiFunc(t_psi,psi,tnom,ynom,td,xd_t,tdotnom,ydotnom,tper,yper,tdotper,ydotper,w1,mag,mass,w,L,ct,cn,Barray,gBarray,tau,epsi)
    
    % use central diff method to find dldx
    dldx = dldxFunc(t_psi,tnom,ynom,td,xd_t,tper,yper,w1)  
%     dldx = zeros(4,1) ; 
    
    % use forward diff method to find dfdx
%     dfdx = dfdxFunc(t_psi,tnom,ynom,tdotnom,ydotnom,tper,yper,tdotper,ydotper,mag,mass,w,L,ct,cn,Barray,gBarray,tau)  
    dfdx = dfdxFunc(t_psi,tnom,ynom,tdotnom,ydotnom,epsi,mag,mass,w,L,ct,cn,Barray,gBarray,tau)  

    t_psi 
    psi 
    % calculate psidot
    psidot = -dldx' - (psi'*dfdx) ;
    psidot = psidot' 
end

function dldx = dldxFunc(t_psi,tnom,ynom,td,xd_t,tper,yper,w1)
    [~,k] = min( abs(tnom-t_psi) )  ; % find the index of the time value closest to t_psi
    xnom = ynom(k,:) ; % this is the nominal value of the state variables at time t_psi 
    
    % find current point in xd_t for calculating cost function
    [~,kd] = min( abs(td-t_psi) ) ; 
    xd_i = xd_t(kd,:) ; 

    % find the current perturbed value
    [~,kper] = min( abs(tper-t_psi) ) ; 
    xper = yper(kper,:) ; 
    
    % need to perturb the costFunc about each state variable individually
    % calculate dl 
    % calculate dldx
    % return vector dldx
    dldx = zeros( size(ynom,2),1 ) ; 
    lnom = costFunc(xnom,xd_i',w1) ; 

    for i = 1:size(ynom,2)
        xper_i = xnom ; 
        xper_i(i) = xper(i) ; 
        lper = costFunc(xper,xd_i',w1) ; 
        dldx(i,:) = ( lper - lnom ) ./ ( xper(i) - xnom(i) )  ; 
    end
    
    dldx(isnan(dldx)) = 0 ; 
    
end

function l = costFunc(x,xd,w) 
    l =  ( w .* (x' - xd)' * (x' - xd) ) ; 
end

function dfdx = dfdxFunc(t_psi,tnom,ynom,tdotnom,ydotnom,epsi,mag,mass,w,L,ct,cn,Barray,gBarray,tau)

    % input correct Barray and gBarray
    if and( (t_psi >= tau(1)), (t_psi <= tau(2)) ) % there must be a better way to do this than hard-coding in the indices of tau?
        Barray_i = Barray{1} ; 
        gBarray_i = gBarray{1} ; 
    elseif and( (t_psi > tau(2)), (t_psi <= tau(3)) )
        Barray_i = Barray{2} ; 
        gBarray_i = gBarray{2} ; 
    elseif and( (t_psi > tau(3)), (t_psi <= tau(4)) )
        Barray_i = Barray{1} ; 
        gBarray_i = gBarray{1} ; 
    elseif and( (t_psi > tau(4)), (t_psi <= tau(5)) )
        Barray_i = 0 ; 
        gBarray_i = 0 ; 
        dfdx = zeros( size(ynom,2), size(ynom,2) ) ;
        return ; 
    end

    dfdx = zeros( size(ynom,2), size(ynom,2) ) ; 

    % find the index k for the nominal time array
    [~,knom] = min( abs(tnom-t_psi) )  ; % find the index of the time value closest to t_psi 
    xnom = ynom(knom,:) ; 

    for i = 1:size(ydotnom,2)
        
        [~,kdotnom] = min( abs(tdotnom-t_psi) ) ; 
        
        if kdotnom < 1
            kdotnom = 1 ; 
        elseif kdotnom > size(tdotnom,2)
            kdotnom = size(tdotnom,2) ; 
        end
        
        fnom = ydotnom(:,kdotnom) ; 
        fnom_i = fnom(i) ; 
        
        for j = 1:size(ydotnom,2)
            
            xper = xnom ; 
            xper(1,i) = xnom(1,i) + epsi(j) ; 
            
            fper = singleDynamicsCoilLtd(t_psi,xper,mag,mass,w,L,ct,cn,Barray_i,gBarray_i) ; 
            fper_j = fper(j) ; 
            
            dfdx(i,j) = (fper_j - fnom_i) / epsi(j) ; 
        end
    end
end

% function dfdx = dfdxFunc(t_psi,tnom,ynom,tdotnom,ydotnom,tper,yper,tdotper,ydotper,mag,mass,w,L,ct,cn,Barray,gBarray,tau)
% 
%     dfdx = zeros( size(ynom,2), size(ynom,2) ) ; 
% 
%     % find the index k for the nominal time array
%     [~,knom] = min( abs(tnom-t_psi) )  ; % find the index of the time value closest to t_psi 
%     xnom = ynom(knom,:) ; 
%     
%     [~,kdotnom] = min( abs(tdotnom-t_psi) ) ; 
%     xdotnom = ydotnom(kdotnom,:) ; 
%     
%     % find the index for the perturbed value
%     [~,kper] = min( abs(tper-t_psi) ) ; 
%     xper = yper(kper,:) ; 
%     
%     [~,kdotper] = min( abs(tdotper-t_psi) ) ;
%     xdotper = ydotper(kdotper,:) ; 
% 
%     for i = 1:size(ydotnom,2)        
%         
%         fnom_i = xdotnom(i) ; 
%         fper_i = xdotper(i) ; 
% 
%         for j = 1:size(ydotnom,2)
%             xnom_j = xnom(j) ; 
%             xper_j = xper(j) ; 
%             
%             dfdx(i,j) = (fper_i - fnom_i) / (xper_j - xnom_j) ; 
%         end
%     end
% end

function psi0 = dmdx(ynom,xd,w2,epsi)

    % find final index of ynom
    k = size(ynom,1) ;
    xnom = ynom(k,:) ; % this is the nominal value of the state variables at time t_psi  
    
    % need to perturb the costFunc about each state variable individually
    % calculate dl 
    % calculate dldx
    % return vector dldx
    psi0 = zeros( size(ynom,2),1 ) ; 
    
    for i = 1:size(ynom,2)
        xper = xnom ; 
        xper(1,i) = xnom(1,i) + epsi(i) ; 
                
        mnom = costFunc(xnom,xd,w2) ; 
        mper = costFunc(xper,xd,w2) ; 
        
        psi0(i,:) = ( mnom - mper ) ./ epsi(i) ; 
    end
    
    psi0(isnan(psi0)) = 0 ; 
end

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

end

function Jdot = Jfunc(t,J,tnom,ynom,td,xd_t,w)
    % find x and xd at time t
    [~,k_i] = min( abs(tnom-t) )  ; % find the index of the time value closest to t_psi 
    x_i = ynom(k_i,:) ; 
    
    [~,kd_i] = min( abs(td-t) )  ; % find the index of the time value closest to t_psi 
    xd_i = xd_t(kd_i,:) ; 

    % calculate the cost function at time t
    Jdot = costFunc(x_i, xd_i', w) ; 
end