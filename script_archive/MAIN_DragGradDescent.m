% MAIN_DragGradDescent.m
% Purpose: Use gradient descent to calculate the drag coefficient for the
% macroswimmer.
% Author: Emma Benjaminson

%% 

clc ; clear all ; close all ; 

%% Directories

localDir = fileparts(mfilename('fullpath')) ;
restoredefaultpath ;
addpath(fullfile(localDir, 'active_functions')) ;

%% Variables

% load macroswimmer data
run VAR_SingleLinkModelCoil.m ; 
load('gBarray.mat') ; 
load('Barray.mat') ; 

% load experimental data
directory = sprintf('%s/TP12',localDir) ; 
% directory = '/home/biorobotics/microswimmers/Gradient Steering/TP12' ; 
filename = 'mpB' ; 
[indices, origins] = extractIndOr(directory,filename) ; 
% filename = 'tp12_mpA_iter' ; 
filename = 'mpB-CompositeData' ; 
loadExpData(indices, origins, directory, filename, 0)

%%%%%%%%%% USER DEFINED VALUE %%%%%%%%%%%%%%%%%%%%%%%%%%%
% select experimental data
% index = 1 ; % index of the trial to plot

% waypoint times - when the swimmer changes direction
% waypointT = [0,225,470,670] ; % [s] 
waypointT = [0,170] ; % [s] 

% n = size(posx{index},1) ; % number of points 
% n = 47 ;

% number of primitives in trajectory
niter = 1 ; 

% initial conditions (USER DEFINES THETA!)
% y0 = [x ; y ; theta ; xdot ; ydot ; thetadot] ; 
y0 = [meanposx(1,1) ; meanposy(1,1) ; pi/2 ; 0 ; 0 ; 0 ]  ;

% control inputs
I1 = [0] ; % [A] current
I2 = [1.9] ; % [A] current
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Gradient Descent Loop

theta_new = 1E-1 ; % starting guess for theta
theta_old = 0 ; % this parameter will store old values for gradient descent update

i = 1 ; % iteration counter
min_iter = 300 ; % minimum number of iterations

% tuning parameters
eps_perturb = 5E-3; % perturbation amount for central diff method [TUNE]
alpha = 0.8 ; % learning rate [TUNE] NOTE (02/25/2020) TRY ALPHA > 0.4
gamma = 0.99 ; 
epsTheta = 1E-6 ; % threshold for when while loop will end [TUNE]

% plot values of theta vs iteration counter
figure(3)
plot(0,theta_new,'ob') 
hold on

% gradient descent while loop
% while and( (abs(theta_old - theta_new) > epsTheta), (i < min_iter) )
while i < min_iter

    theta_old = theta_new ; % save new value of theta to be used in this iteration
    
    %%%%%%%%%%%%
    % Hypothesis
    %%%%%%%%%%%%
    
    disp('\n RUNNING CALCULATION FOR HYPOTHESIS. \n') ; 

    % run ode15s for the same time span as the experimental data that we have,
    % and the current drag coefficient
    [x,y] = hsim(niter,waypointT,y0,mag,mass,w,L,theta_old,Barray{2},gBarray{2}) ;
    
    % plot hypothesis and variations
    figure(1)
    plot(x,y,'-k')
    hold on

    %%%%%%%%%%%%%%%
    % Cost Function
    %%%%%%%%%%%%%%%
    J = costFunc([x,y],[meanposx,meanposy]) ; 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Derivative of Cost Function
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % use central difference method
    % need to calculate h for two different values of theta slightly perturbed
    % from current theta then divide by 2 * delta theta

    disp('\n RUNNING CALCULATION FOR H PLUS. \n') ; 
    
    % plus
    theta_plus = theta_old + eps_perturb ; 
    [xplus,yplus] = hsim(niter,waypointT,y0,mag,mass,w,L,theta_plus,Barray{2},gBarray{2}) ; 
    J_plus = costFunc([xplus,yplus],[meanposx,meanposy]) ; 
    
    disp('\n RUNNING CALCULATION FOR H MINUS. \n') ; 
    
    % minus
    theta_minus = theta_old - eps_perturb ; 
    [xminus,yminus]  = hsim(niter,waypointT,y0,mag,mass,w,L,theta_minus,Barray{2},gBarray{2}) ; 
    J_minus = costFunc([xminus,yminus],[meanposx,meanposy]) ; 
    
    % plot variations on hypothesis
    figure(1)
    plot(xminus,yminus,'-r') 
    hold on 
    plot(xplus,yplus,'-b')
    hold on
    
    % plot cost function values
    figure(2)
    plot(i,J,'ok')
    hold on
    plot(i,J_plus,'ob') 
    hold on
    plot(i,J_minus,'or') 

    dJdtheta = (J_plus - J_minus) / (2*eps_perturb) ; 
    
    %%%%%%%%%%%%%%
    % Choose Alpha
    %%%%%%%%%%%%%%
    Jnew = 1E9 ; 
    iter = 2 ; 
    
    while Jnew > J
        
        % calculate new value of theta
        theta_new = theta_old - alpha*dJdtheta  
        
        % forward simulate dynamics
        disp('Calculating dynamics to update alpha') ; 
        [x,y] = hsim(niter,waypointT,y0,mag,mass,w,L,theta_new,Barray{2},gBarray{2}) ;
        
        % calculate new cost function
        Jnew = costFunc([x,y],[meanposx,meanposy]) ; 
        
        if Jnew > J
            alpha = alpha^iter
        end
            
        iter = iter + 1  
    end 
        
    %%%%%%%%%%%%%%
    % Update Theta
    %%%%%%%%%%%%%%
    theta_new = theta_old - alpha*dJdtheta  
    
    % plot theta value
    figure(3)
    plot(i,theta_new,'ok')
    hold on
    
    i = i + 1  % update counter
    
end

% add titles
figure(1)
title('Trajectories') 

figure(2) 
title('J values') ; 

figure(3) 
title('theta values') ; 
