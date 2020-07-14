% MAIN_PlotDynamics.m
% Purpose: Plot the dynamics data from simulation in various
% configurations.
% Author: Emma Benjaminson
% References: 

%% 

clc ; clear all ; close all ; 

%% Directories

localDir = fileparts(mfilename('fullpath')) ;
restoredefaultpath ;
addpath(fullfile(localDir, 'active_functions')) ;

%% Load Data



% %% FIGURE 1: plot trajectory, velocity and rotation 
% 
% %%%%%%%%%%%%
% % trajectory
% %%%%%%%%%%%%
% figure(1)
% subplot(2,2,1)
% plot(posx{index},posy{index},'-ob')
% hold on
% 
% figure(1)
% subplot(2,2,1)
% 
% for i = 1:niter
%     plot(yarray{i}(:,1),yarray{i}(:,2),'-or') 
%     hold on 
% end
% 
% xlabel('X (m)','interpreter','latex') ; 
% ylabel('Y (m)','interpreter','latex') ; 
% title('Trajectory','interpreter','latex') ; 
% legend('Exp Data','Model','interpreter','latex') ; 
% 
% %%%%%%%
% % theta
% %%%%%%%
% figure(1)
% subplot(2,2,2)
% plot(t,y(:,3)) 
% hold on
% xlabel('Time (s)','interpreter','latex') ; 
% ylabel('$$\theta$$ (rad)','interpreter','latex') ; 
% title('$$\theta$$','interpreter','latex') ; 
% legend('Model','interpreter','latex') ; 
% 
% %%%%%%%%%%
% % velocity
% %%%%%%%%%%
% % experimental data
% ex = posx{index} ; 
% ey = posy{index} ; 
% et = time{index} ; 
% 
% % instantaneous velocity
% vx = diff(ex) ./ diff(et) ; % [m/s]
% vy = diff(ey) ./ diff(et) ; 
% 
% % add zeros to pad velocity arrays
% vx = [0 ; vx] ; 
% vy = [0 ; vy] ; 
% 
% figure(1)
% subplot(2,2,3) 
% plot(et(1:47),vx(1:47),'--xg') ; % experimental x velocity
% hold on
% plot(et(1:47),vy(1:47),'--og') ; % experimental y velocity
% hold on
% plot(t,y(:,4),'-xb') ; % model x velocity
% hold on
% plot(t,y(:,5),'-ob') ; % model y velocity
% hold on 
% xlabel('Time (s)','interpreter','latex') ; 
% ylabel('Velocity (m/s)','interpreter','latex') ; 
% title('Velocities') ; 
% legend('Exp vx','Exp vy','Model vx','Model vy','interpreter','latex','Location','northwest') ; 
% 
% %%%%%%%%%%
% % thetadot
% %%%%%%%%%%
% figure(1) 
% subplot(2,2,4) 
% plot(t,y(:,6),'-b') ; 
% hold on 
% xlabel('Time (s)','interpreter','latex') ; 
% ylabel('$$\dot{\theta}$$ (rad/s)','interpreter','latex') ; 
% title('$$\dot{\theta}$$','interpreter','latex') ; 
% legend('Model', 'interpreter','latex') ; 
% 
% %% FIGURE 2: trajectory only
% 
% figure(2) 
% plot(posx{index},posy{index},'-ob')
% hold on
% plot(yarray{1}(:,1),yarray{1}(:,2),'-r') 
% hold on
% plot(yarray{2}(:,1),yarray{2}(:,2),'-r','HandleVisibility','off') 
% hold on
% plot(yarray{3}(:,1),yarray{3}(:,2),'-r','HandleVisibility','off') 
% hold on
% xlabel('X (m)','interpreter','latex') ; 
% ylabel('Y (m)','interpreter','latex') ; 
% title('Trajectory','interpreter','latex','FontSize',14) ; 
% legend('Exp Data','Model','interpreter','latex') ;

% %% FIGURE 3: velocity only
% 
% figure(3)
% plot(et,vx,'--xg') ; % experimental x velocity
% hold on
% plot(et,vy,'--og') ; % experimental y velocity
% hold on
% plot(timearray{1},yarray{1}(:,4),'-xb') ; % model x velocity
% hold on
% plot(timearray{1},yarray{1}(:,5),'-ob') ; % model y velocity
% hold on 
% plot(timearray{2},yarray{2}(:,4),'-xb') ; % model x velocity
% hold on
% plot(timearray{2},yarray{2}(:,5),'-ob') ; % model y velocity
% hold on 
% plot(timearray{3},yarray{3}(:,4),'-xb') ; % model x velocity
% hold on
% plot(timearray{3},yarray{3}(:,5),'-ob') ; % model y velocity
% hold on 
% xlabel('Time (s)','interpreter','latex') ; 
% ylabel('Velocity (m/s)','interpreter','latex') ; 
% title('Velocities') ; 
% legend('Exp vx','Exp vy','Model vx','Model vy','interpreter','latex','Location','northwest') ; 
% 
%% FIGURE 4: Theta and Thetadot

% figure(4)
% 
% yyaxis left
% plot(timearray{1},yarray{1}(:,3),'-xb') ; % model theta position
% hold on
% plot(timearray{2},yarray{2}(:,3),'-xb','HandleVisibility','off') ; % model theta position
% hold on
% plot(timearray{3},yarray{3}(:,3),'-xb','HandleVisibility','off') ; % model theta position
% hold on
% ylabel('$$\theta$$ (rad)','interpreter','latex') ; 
% 
% yyaxis right
% plot(timearray{1},yarray{1}(:,6),'-or') ; % model theta velocity
% hold on
% plot(timearray{2},yarray{2}(:,6),'-or','HandleVisibility','off') ; % model theta velocity
% hold on
% plot(timearray{3},yarray{3}(:,6),'-or','HandleVisibility','off') ; % model theta velocity
% hold on
% xlabel('Time (s)','interpreter','latex') ; 
% ylabel('$$\dot{\theta}$$ (rad/s)','interpreter','latex') ; 
% title('Rotational Motion','interpreter','latex') ; 
% legend('$$\theta$$','$$\dot{\theta}$$', 'interpreter','latex') ; 

% %% FIGURE 5: Swimmer with orientation overlaid (USED TO GENERATE FIGURE 1 IN IROS 2020)
% load('MAT_SingleLinkModelCoil_GOOD') ; 
% nfig = 5 ; 
% 
% % indices of where to plot the swimmer orientation
% yind1 = [1,50,80,80*5,640,720,800] ; 
% yind2 = [50,50*5,350,400,500,550] ; 
% yind3 = [25,50,100,550,600,700] ; 
% 
% % import s1-CompositeData for plotting experimental data
% data = readtable('/home/emma/Documents/Research/TP11/Sequence 1/5 sec_frame Composites/s1-CompositeData.csv') ; 
% data = table2array( data(1:end, 1:(end-1)) ) ; 
% 
% % values for plotting experimental data
% meanposx = data(:,1)*1E-3 ; 
% meanposy = data(:,2)*1E-3 ; 
% stdposx = data(:,7)*1E-3 ; 
% stdposy = data(:,8)*1E-3 ; 
% ntraj = 3 ; 
% wp = [1,26,51,73] ; % include start and end points
% % color = '#90aad1' ; 
% color = '#52639e' ; 
% 
% % gradient plotting
% xmingrad = -0.04 ; % [m]
% xmaxgrad = 0.04 ; % [m]
% ymingrad = 0.015 ; % [m]
% ymaxgrad = 0.046 ; % [m]
% positionArrayGrad = [dx ; xmingrad ; xmaxgrad ; ymingrad ; ymaxgrad] ; 
% 
% % calculate gradient
% coils1 = [1.9 ; 0 ; a ; d ; nturns] ;
% n = 10 ; % size of position vector
% [gx2, gy2, gBx, gBy, ~] = magGradientCoil(coils1,mu0,positionArrayGrad,n,gBmax) ; 
% 
% 
% figure(nfig) 
% plotExpData(meanposx, meanposy, stdposx, stdposy, nfig, ntraj, wp, color) ; 
% hold on 
% quiver(gx2,gy2,gBx,gBy,1,'-k') 
% hold on
% plot(yarray{1}(:,1),yarray{1}(:,2),'-r','LineWidth',2) 
% hold on
% plot(yarray{2}(:,1),yarray{2}(:,2),'-r','HandleVisibility','off','LineWidth',2) 
% hold on
% plot(yarray{3}(:,1),yarray{3}(:,2),'-r','HandleVisibility','off','LineWidth',2) 
% hold on
% plotSwimmer(nfig,yarray{1},yind1,1) ; 
% hold on
% plotSwimmer(nfig,yarray{2},yind2,0) ; 
% hold on
% plotSwimmer(nfig,yarray{3},yind3,1) ; 
% hold on
% xlabel('Z (m)','interpreter','latex','FontSize',32) ; 
% ylabel('X (m)','interpreter','latex','FontSize',32) ; 
% % title('Trajectory','interpreter','latex','FontSize',14) ; 
% legend('Experimental Data','Gradient Field','Model','interpreter','latex','FontSize',32) ;
% legend boxon
% axis([-0.03,xmaxgrad,ymingrad,ymaxgrad]) ; 

% % FIGURE 6: Swimmer with orientation overlaid (USED TO GENERATE FIGURE 5 IN IROS 2020)
% load('MAT_SingleLinkModelCoil_mpA') ; 
% nfig = 6 ; 
% 
% indices of where to plot the swimmer orientation
% yind1 = [1,10,20,30,800,1000,1500,2000] ; 
% 
% % import s1-CompositeData for plotting experimental data
% data = readtable('/home/emma/repos/microswimmers/Gradient Steering/TP12/mpA-CompositeData.csv') ; 
% data = table2array( data(1:end, 1:(end-1)) ) ; 
% 
% color = '#52639e' ; 
% ntraj = 1 ; 
% wp = [1,34] ; 
% 
% gradient plotting
% xmingrad = -0.04 ; % [m]
% xmaxgrad = 0.04 ; % [m]
% ymingrad = 0.015 ; % [m]
% ymaxgrad = 0.046 ; % [m]
% positionArrayGrad = [dx ; xmingrad ; xmaxgrad ; ymingrad ; ymaxgrad] ; 
% 
% calculate gradient
% coils1 = [1.9 ; 0 ; a ; d ; nturns] ;
% n = 10 ; % size of position vector
% [gx2, gy2, gBx, gBy, ~] = magGradientCoil(coils1,mu0,positionArrayGrad,n,gBmax) ; 
% 
% 
% figure(nfig) 
% plotExpData(meanposx, meanposy, stdposx, stdposy, nfig, ntraj, wp, color) ; 
% hold on 
% quiver(gx2,gy2,gBx,gBy,1,'-k') 
% hold on
% plot(yarray{1}(:,1),yarray{1}(:,2),'-r','LineWidth',2) 
% hold on
% plotSwimmer(nfig,yarray{1},yind1,1) ; 
% hold on
% ax = gca ; 
% ax.FontSize = 14 ; 
% xlabel('Z (m)','interpreter','latex','FontSize',32) ; 
% ylabel('X (m)','interpreter','latex','FontSize',32) ; 
% title('Motion Primitive Data','interpreter','latex','FontSize',32) ; 
% legend('Experimental Data','Gradient Field','Model','interpreter','latex','FontSize',20) ;
% legend boxon
% grid on
% axis([-0.03,xmaxgrad,ymingrad,ymaxgrad]) ; 

%% FIGURE 7: Swimmer with orientation overlaid (USED TO GENERATE FIGURE 6 IN IROS 2020)

load('MAT_SingleLinkModelCoil_case1') ; 
nfig = 7 ; 

% indices of where to plot the swimmer orientation
yind1 = [10,20,30,500,700] ; 
yind2 = [400,750,950] ; 
yind3 = [1,600,660,680] ; 

% % import s1-CompositeData for plotting experimental data
% data = readtable('/home/emma/repos/microswimmers/Gradient Steering/TP12/mpA-CompositeData.csv') ; 
% data = table2array( data(1:end, 1:(end-1)) ) ; 

color = '#52639e' ; 
ntraj = 1 ; 
wp = [1,13,22,40] ; 

% gradient plotting
xmingrad = -0.04 ; % [m]
xmaxgrad = 0.04 ; % [m]
ymingrad = 0.015 ; % [m]
ymaxgrad = 0.046 ; % [m]
positionArrayGrad = [dx ; xmingrad ; xmaxgrad ; ymingrad ; ymaxgrad] ; 

% calculate gradient
% coils1 = [1.9 ; 0 ; a ; d ; nturns] ;
% n = 10 ; % size of position vector
% [gx2, gy2, gBx, gBy, ~] = magGradientCoil(coils1,mu0,positionArrayGrad,n,gBmax) ; 


figure(nfig) 
plotExpData(meanposx, meanposy, stdposx, stdposy, nfig, ntraj, wp, color) ; 
hold on 
% quiver(gx2,gy2,gBx,gBy,1,'-k') 
% hold on
plot(yarray{1}(:,1),yarray{1}(:,2),'-r','LineWidth',2) 
hold on
% plot(yarray{2}(:,1),yarray{2}(:,2),'-r','HandleVisibility','off','LineWidth',2) 
% hold on
% plot(yarray{3}(:,1),yarray{3}(:,2),'-r','HandleVisibility','off','LineWidth',2) 
% hold on
% plotSwimmer(nfig,yarray{1},yind1,1) ; 
% hold on
% plotSwimmer(nfig,yarray{2},yind2,0) ; 
% hold on
% plotSwimmer(nfig,yarray{3},yind3,1) ; 
% hold on
ax = gca ; 
ax.FontSize = 14 ; 
xlabel('Z (m)','interpreter','latex','FontSize',32) ; 
ylabel('X (m)','interpreter','latex','FontSize',32) ; 
title('Case 1: Results','interpreter','latex','FontSize',32) ; 
legend('Experimental Data','Model','interpreter','latex','FontSize',20) ;
legend boxon
grid on
% axis([-0.03,xmaxgrad,ymingrad,ymaxgrad]) ; 

% %% FIGURE 8: Swimmer with orientation overlaid (USED TO GENERATE FIGURE 1 IN IROS 2020)
% 
% load('MAT_SingleLinkModelCoil_case1') ; 
% nfig = 1 ; 
% 
% % indices of where to plot the swimmer orientation
% yind1 = [10,20,30,500,700] ; 
% yind2 = [400,750,950] ; 
% yind3 = [1,600,660,680] ; 
% 
% 
% % % import s1-CompositeData for plotting experimental data
% % data = readtable('/home/emma/repos/microswimmers/Gradient Steering/TP12/mpA-CompositeData.csv') ; 
% % data = table2array( data(1:end, 1:(end-1)) ) ; 
% 
% color = '#52639e' ; 
% ntraj = 1 ; 
% wp = [1,13,22,40] ; 
% 
% % gradient plotting
% xmingrad = -0.04 ; % [m]
% xmaxgrad = 0.04 ; % [m]
% ymingrad = 0.015 ; % [m]
% ymaxgrad = 0.046 ; % [m]
% positionArrayGrad = [dx ; xmingrad ; xmaxgrad ; ymingrad ; ymaxgrad] ; 
% 
% % calculate gradient
% coils1 = [1.9 ; 0 ; a ; d ; nturns] ;
% n = 10 ; % size of position vector
% [gx2, gy2, gBx, gBy, ~] = magGradientCoil(coils1,mu0,positionArrayGrad,n,gBmax) ; 
% 
% 
% figure(nfig) 
% plotExpData(meanposx(1:wp(2),1), meanposy(1:wp(2),1), stdposx(1:wp(2),1), stdposy(1:wp(2),1), nfig, ntraj, wp, color) ; 
% hold on 
% quiver(gx2,gy2,gBx,gBy,1,'-k') 
% hold on
% plot(yarray{1}(:,1),yarray{1}(:,2),'-r','LineWidth',2) 
% hold on
% plot(yarray{2}(:,1),yarray{2}(:,2),'-r','HandleVisibility','off','LineWidth',2) 
% hold on
% plot(yarray{3}(:,1),yarray{3}(:,2),'-r','HandleVisibility','off','LineWidth',2) 
% hold on
% plotSwimmer(nfig,yarray{1},yind1,1) ; 
% hold on
% plotSwimmer(nfig,yarray{2},yind2,0) ; 
% hold on
% plotSwimmer(nfig,yarray{3},yind3,1) ; 
% hold on
% ax = gca ; 
% ax.FontSize = 14 ; 
% xlabel('Z (m)','interpreter','latex','FontSize',32) ; 
% ylabel('X (m)','interpreter','latex','FontSize',32) ; 
% title('Motion Primitive 1 for $$0 < t < t_1$$','interpreter','latex','FontSize',32) ; 
% legend('Experimental Data','Gradient Field','Model','interpreter','latex','FontSize',20) ;
% legend boxon
% grid on
% axis([-0.03,0.03,ymingrad,ymaxgrad]) ; 
% 
% %% FIGURE 9: Swimmer with orientation overlaid (USED TO GENERATE FIGURE 1 IN IROS 2020)
% 
% load('MAT_SingleLinkModelCoil_case1') ; 
% nfig = 9 ; 
% 
% % indices of where to plot the swimmer orientation
% yind1 = [10,20,30,500,700] ; 
% yind2 = [400,750,950] ; 
% yind3 = [1,600,660,680] ; 
% 
% % % import s1-CompositeData for plotting experimental data
% % data = readtable('/home/emma/repos/microswimmers/Gradient Steering/TP12/mpA-CompositeData.csv') ; 
% % data = table2array( data(1:end, 1:(end-1)) ) ; 
% 
% color = '#52639e' ; 
% ntraj = 2 ; 
% wp = [1,13,22,40] ; 
% 
% % gradient plotting
% xmingrad = -0.04 ; % [m]
% xmaxgrad = 0.04 ; % [m]
% ymingrad = 0.015 ; % [m]
% ymaxgrad = 0.046 ; % [m]
% positionArrayGrad = [dx ; xmingrad ; xmaxgrad ; ymingrad ; ymaxgrad] ; 
% 
% % calculate gradient
% coils1 = [0 ; 1.9 ; a ; d ; nturns] ;
% n = 10 ; % size of position vector
% [gx2, gy2, gBx, gBy, ~] = magGradientCoil(coils1,mu0,positionArrayGrad,n,gBmax) ; 
% 
% 
% figure(nfig) 
% plotExpData(meanposx(1:wp(3),1), meanposy(1:wp(3),1), stdposx(1:wp(3),1), stdposy(1:wp(3),1), nfig, ntraj, wp, color) ; 
% hold on 
% quiver(gx2,gy2,gBx,gBy,1,'-k') 
% hold on
% plot(yarray{1}(:,1),yarray{1}(:,2),'-r','LineWidth',2) 
% hold on
% plot(yarray{2}(:,1),yarray{2}(:,2),'-r','HandleVisibility','off','LineWidth',2) 
% hold on
% % plot(yarray{3}(:,1),yarray{3}(:,2),'-r','HandleVisibility','off','LineWidth',2) 
% % hold on
% plotSwimmer(nfig,yarray{1},yind1,1) ; 
% hold on
% plotSwimmer(nfig,yarray{2},yind2,0) ; 
% hold on
% % plotSwimmer(nfig,yarray{3},yind3,1) ; 
% % hold on
% ax = gca ; 
% ax.FontSize = 14 ; 
% xlabel('Z (m)','interpreter','latex','FontSize',32) ; 
% ylabel('X (m)','interpreter','latex','FontSize',32) ; 
% title('Motion Primitive 2 for $$t_1 < t < t_2$$','interpreter','latex','FontSize',32) ; 
% legend('Experimental Data','Gradient Field','Model','interpreter','latex','FontSize',20) ;
% legend boxon
% grid on
% axis([-0.03,0.03,ymingrad,ymaxgrad]) ; 
% 
% %% FIGURE 10: Swimmer with orientation overlaid (USED TO GENERATE PART 1 INSET IN IROS 2020 VIDEO)
% load('MAT_SingleLinkModelCoil_mpA') ; 
% nfig = 7 ; 
% 
% % indices of where to plot the swimmer orientation
% yind1 = [10,20,800,1000, 1500,2000] ; 
% 
% % % import s1-CompositeData for plotting experimental data
% % data = readtable('/home/emma/repos/microswimmers/Gradient Steering/TP12/mpA-CompositeData.csv') ; 
% % data = table2array( data(1:end, 1:(end-1)) ) ; 
% 
% color = '#52639e' ; 
% ntraj = 1 ; 
% wp = [1,34] ; 
% 
% % gradient plotting
% xmingrad = -0.04 ; % [m]
% xmaxgrad = 0.04 ; % [m]
% ymingrad = 0.015 ; % [m]
% ymaxgrad = 0.046 ; % [m]
% positionArrayGrad = [dx ; xmingrad ; xmaxgrad ; ymingrad ; ymaxgrad] ; 
% 
% % calculate gradient
% coils1 = [1.9 ; 0 ; a ; d ; nturns] ;
% n = 10 ; % size of position vector
% [gx2, gy2, gBx, gBy, ~] = magGradientCoil(coils1,mu0,positionArrayGrad,n,gBmax) ; 
% 
% 
% figure(nfig) 
% box on
% plotExpData(meanposx(1:34,1), meanposy(1:34,1), stdposx(1:34,1), stdposy(1:34,1), nfig, ntraj, wp, color) ; 
% hold on 
% quiver(gx2,gy2,gBx,gBy,1,'-k') 
% hold on
% plot(yarray{1}(1:end,1),yarray{1}(1:end,2),'-r','LineWidth',2) 
% hold on
% plotSwimmer(nfig,yarray{1},yind1,1) ; 
% hold on
% ax = gca ; 
% ax.FontSize = 14 ; 
% xlabel('Z (m)','interpreter','latex','FontSize',32) ; 
% ylabel('X (m)','interpreter','latex','FontSize',32) ; 
% % title('Trajectory','interpreter','latex','FontSize',14) ; 
% legend('Experimental Data','Gradient Field','Model','interpreter','latex','FontSize',20) ;
% legend boxon
% grid on
% 
% % axis([-0.03,xmaxgrad,ymingrad,ymaxgrad]) ; 
% 
% %% FIGURE 1: Swimmer with orientation overlaid (USED TO GENERATE PART 2 IN IROS 2020 VIDEO)
% 
% load('MAT_SingleLinkModelCoil_case1') ; 
% nfig = 1 ; 
% 
% % indices of where to plot the swimmer orientation
% yind1 = [10,20,30,700] ; 
% yind2 = [200] %,750] ; %,950] ; 
% yind3 = [] %,680] ; 
% 
% % % import s1-CompositeData for plotting experimental data
% % data = readtable('/home/emma/repos/microswimmers/Gradient Steering/TP12/mpA-CompositeData.csv') ; 
% % data = table2array( data(1:end, 1:(end-1)) ) ; 
% 
% color = '#52639e' ; 
% ntraj = 2 ; 
% wp = [1,13,16] ;
% % wp = [1,13,22,40] ; 
% 
% % gradient plotting
% xmingrad = -0.06 ; % [m]
% xmaxgrad = 0.04 ; % [m]
% ymingrad = 0.015 ; % [m]
% ymaxgrad = 0.046 ; % [m]
% positionArrayGrad = [dx ; xmingrad ; xmaxgrad ; ymingrad ; ymaxgrad] ; 
% 
% % calculate gradient
% coils1 = [0 ; 1.9 ; a ; d ; nturns] ;
% n = 10 ; % size of position vector
% [gx2, gy2, gBx, gBy, ~] = magGradientCoil(coils1,mu0,positionArrayGrad,n,gBmax) ; 
% 
% 
% figure(nfig) 
% plotExpData(meanposx(1:16,1), meanposy(1:16,1), stdposx(1:16,1), stdposy(1:16,1), nfig, ntraj, wp, color) ; 
% hold on 
% quiver(gx2,gy2,gBx,gBy,1,'-k') 
% hold on
% plot(yarray{1}(1:end,1),yarray{1}(1:end,2),'-r','LineWidth',2) 
% hold on
% plot(yarray{2}(1:500,1),yarray{2}(1:500,2),'-r','HandleVisibility','off','LineWidth',2) 
% hold on
% % plot(yarray{3}(1:400,1),yarray{3}(1:400,2),'-r','HandleVisibility','off','LineWidth',2) 
% % hold on
% plotSwimmer(nfig,yarray{1},yind1,1) ; 
% hold on
% plotSwimmer(nfig,yarray{2},yind2,0) ; 
% hold on
% plotSwimmer(nfig,yarray{3},yind3,1) ; 
% hold on
% ax = gca ; 
% ax.FontSize = 14 ; 
% xlabel('Z (m)','interpreter','latex','FontSize',32) ; 
% ylabel('X (m)','interpreter','latex','FontSize',32) ; 
% legend('Experimental Data','Gradient Field','Model','interpreter','latex','FontSize',20) ;
% legend boxon
% grid on
% axis([-0.04,0.05,-0.05,0.05]) ; 

% %% FIGURE 6: Plot MP1
% 
% nfig = 6 ; 
% 
% coils1 = [1.9 ; 0 ; a ; d ; nturns] ;
% 
% % gradient plotting
% xmingrad = -0.04 ; % [m]
% xmaxgrad = 0.04 ; % [m]
% ymingrad = 0.022 ; % [m]
% ymaxgrad = 0.046 ; % [m]
% positionArrayGrad = [dx ; xmingrad ; xmaxgrad ; ymingrad ; ymaxgrad] ; 
% 
% % calculate gradient
% n = 10 ; % size of position vector
% [gx2, gy2, gBx, gBy, ~] = magGradientCoil(coils1,mu0,positionArrayGrad,n,gBmax) ; 
% 
% figure(nfig)
% quiver(gx2,gy2,gBx,gBy,3,'HandleVisibility','off') 
% hold on
% plot(posx{index}(1:45),posy{index}(1:45),'-ob') % experimental data
% hold on
% plot(yarray{1}(:,1),yarray{1}(:,2),'-r') 
% hold on
% plotSwimmer(nfig,yarray{1},yind1,1) ; 
% hold on
% axis([xmingrad,xmaxgrad,ymingrad,ymaxgrad]) ; 
% xlabel('X (m)','interpreter','latex') ; 
% ylabel('Y (m)','interpreter','latex') ; 
% title('Motion Primitive 1','interpreter','latex','FontSize',14) ; 
% legend('Exp Data','Model','interpreter','latex') ;
% 
% %% FIGURE 7: Plot MP2
% 
% nfig = 7 ; 
% 
% coils1 = [0 ; 1.9 ; a ; d ; nturns] ;
% 
% % gradient plotting
% xmingrad = -0.03 ; % [m]
% xmaxgrad = 0.04 ; % [m]
% ymingrad = 0.014 ; % [m]
% ymaxgrad = 0.028 ; % [m]
% positionArrayGrad = [dx ; xmingrad ; xmaxgrad ; ymingrad ; ymaxgrad] ; 
% 
% % calculate gradient
% n = 10 ; % size of position vector
% [gx2, gy2, gBx, gBy, ~] = magGradientCoil(coils1,mu0,positionArrayGrad,n,gBmax) ; 
% 
% figure(nfig)
% quiver(gx2,gy2,gBx,gBy,1,'HandleVisibility','off') 
% hold on
% plot(posx{index}(45:94),posy{index}(45:94),'-ob') % experimental data
% hold on
% plot(yarray{2}(:,1),yarray{2}(:,2),'-r') 
% hold on
% plotSwimmer(nfig,yarray{2},yind2,0) ; 
% hold on
% axis([xmingrad,xmaxgrad,ymingrad,ymaxgrad]) ; 
% xlabel('X (m)','interpreter','latex') ; 
% ylabel('Y (m)','interpreter','latex') ; 
% title('Motion Primitive 2','interpreter','latex','FontSize',14) ; 
% legend('Exp Data','Model','interpreter','latex') ;
% 
% %% FIGURE 8: Plot MP3
% 
% nfig = 8 ; 
% 
% coils1 = [1.9 ; 0 ; a ; d ; nturns] ;
% 
% % gradient plotting
% xmingrad = -0.03 ; % [m]
% xmaxgrad = 0.05 ; % [m]
% ymingrad = 0.005 ; % [m]
% ymaxgrad = 0.02 ; % [m]
% positionArrayGrad = [dx ; xmingrad ; xmaxgrad ; ymingrad ; ymaxgrad] ; 
% 
% % calculate gradient
% n = 10 ; % size of position vector
% [gx2, gy2, gBx, gBy, ~] = magGradientCoil(coils1,mu0,positionArrayGrad,n,gBmax) ; 
% 
% figure(nfig)
% quiver(gx2,gy2,gBx,gBy,1,'HandleVisibility','off') 
% hold on
% plot(posx{index}(94:134),posy{index}(94:134),'-ob') % experimental data
% hold on
% plot(yarray{3}(:,1),yarray{3}(:,2),'-r') 
% hold on
% plotSwimmer(nfig,yarray{3},yind3,0) ; 
% hold on
% axis([xmingrad,xmaxgrad,ymingrad,ymaxgrad]) ; 
% xlabel('X (m)','interpreter','latex') ; 
% ylabel('Y (m)','interpreter','latex') ; 
% title('Motion Primitive 3','interpreter','latex','FontSize',14) ; 
% legend('Exp Data','Model','interpreter','latex') ;

%% PLOT RESULTS
% 
% load('MAT_SingleLinkModelCoil_case3') ; 
% 
% nfig = 10 ; 
% 
% figure(nfig) 
% plot(yarray{1}(:,1),yarray{1}(:,2),'-r','LineWidth',2) 
% hold on
% plot(yarray{2}(:,1),yarray{2}(:,2),'-r','LineWidth',2,'HandleVisibility','off') 
% hold on
% plot(yarray{3}(:,1),yarray{3}(:,2),'-r','LineWidth',2,'HandleVisibility','off') 
% hold on
% plot(yarray{4}(:,1),yarray{4}(:,2),'-r','LineWidth',2,'HandleVisibility','off') 
% hold on
% plot(yarray{5}(:,1),yarray{5}(:,2),'-r','LineWidth',2,'HandleVisibility','off') 
% hold on
% plot(yarray{6}(:,1),yarray{6}(:,2),'-r','LineWidth',2,'HandleVisibility','off') 
% hold on
% plot(yarray{7}(:,1),yarray{7}(:,2),'-r','LineWidth',2,'HandleVisibility','off') 
% hold on
% ax = gca ; 
% ax.FontSize = 14 ; 
% xlabel('Z (m)','interpreter','latex','FontSize',32) ; 
% ylabel('X (m)','interpreter','latex','FontSize',32) ; 
% title('Case 2: Results','interpreter','latex','FontSize',32) ; 
% legend('Model','interpreter','latex','FontSize',20) ;
% legend boxon
% grid on