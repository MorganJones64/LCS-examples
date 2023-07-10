
% Example code of Finite-Time-Lyapunov Exponent analysis for an Unsteady Double-Gyre Flow
% Derived from ME 564 - Mechanical Engineering Analysis:
% Lecturer Steven Brunton, University of Washington
% (http://faculty.washington.edu/sbrunton/me564/) 
% (https://cassyni.com/events/PZ37zqJFMgX6XS26KWKR6V)
% Author: Steve Brunton, Morgan Jones
%% Initiate constants for integration and dynamical system
clear all
close all
addpath('./customcolormap')
tstart = tic;
A = 0.1;    % parameters from Shadden 2005 Physica D
eps = 0.25;
omega = 2*pi/10;  % frequency of gyre oscillations
dt =0.025;  % timestep
T = 15;     % duration of integration (higher number = finer ridges)
int = 'f'; %''f for forward integration, 'b' for backward integration;
%% Part 1 - Compute one FTLE field: Initialize grid of particles through vector field
vs=0.6;
dx = .025; % Grid Resolution (smaller number = finer ridges)
xvec = 0:dx:2;
yvec = 0:dx:1;
[x0,y0] = meshgrid(xvec,yvec);  % grid of particles
yIC(1,:,:) = x0';
yIC(2,:,:) = y0';

% Intiate Vector Field and Grid of Particles
figure(1)
subplot(2,1,1)
dy = doublegyreVEC(0,yIC,A,eps,omega);
% plot the vector field using 'quiver'
quiver(yIC(1,1:4:end,1:4:end),yIC(2,1:4:end,1:4:end),dy(1,1:4:end,1:4:end),dy(2,1:4:end,1:4:end),vs,'Color','#ababab','linewidth',1.8);
axis([0 2 0 1]), drawnow
subplot(2,1,2)
% plot initial conditions
plot(yIC(1,:),yIC(2,:),'r.','LineWidth',2,'MarkerSize',4)
axis([0 2 0 1]), drawnow
set(gcf,'Position',[100 100 600 600])
set(gcf,'color','w')
%% Part 2 - Compute trajectory 
if int == 'f' % if the integration is set in forward time
    sgn = 1; % use positive time integration
else % if integration is set in backwards time
    sgn = -1; % use negative time integration
end
yin = yIC;
for i=flip(1:T/dt)
    figure(1)
    time = i*dt;
    subplot(2,1,1)
    dy = doublegyreVEC(time,yIC,A,eps,omega);
    % evolve vector field in time
    quiver(yIC(1,1:4:end,1:4:end),yIC(2,1:4:end,1:4:end),dy(1,1:4:end,1:4:end),dy(2,1:4:end,1:4:end),vs,'Color','#ababab','linewidth',1.8);
    axis([0 2 0 1])
    drawnow
    % compute trajectory
    yout = rk4singlestep(@(t,y)doublegyreVEC(t,y,A,eps,omega),sgn*dt,sgn*time,yin);
    yin = yout;
    
    % plot trajectory
    subplot(2,1,2)
    plot(yout(1,:),yout(2,:),'r.','LineWidth',2,'MarkerSize',4)
    axis([0 2 0 1])     
end

% reshape 3-dim array into 2-dim array
xT = reshape(yout(1,:,:),length(xvec),length(yvec));
yT = reshape(yout(2,:,:),length(xvec),length(yvec));

%% Part 3 -  Compute the finite-time Lyapunov exponent (sigma)
% Finite difference to compute the gradient
[dxTdx0,dxTdy0] = gradient(xT,dx,dx);
[dyTdx0,dyTdy0] = gradient(yT,dx,dx);
if int == 'f'
    mycolormap = customcolormap(linspace(0,1,6),{'#a9212d','#b8412a','#ca6827','#edb121','#f5d586','#ffffff'});
else
    mycolormap = customcolormap(linspace(0,1,6),{'#132d8d','#0d449d','#0758ab','#2295e0','#a9d9f8','#ffffff'});
end
% compute sigma: large sigma indicates large mixing!
for i=1:length(xvec)
    for j=1:length(yvec)
        D(1,1) = dxTdx0(i,j);
        D(1,2) = dxTdy0(i,j);
        D(2,1) = dyTdx0(i,j);
        D(2,2) = dyTdy0(i,j);
        sigma(i,j) = (1/T)*sqrt(max(eig(D'*D)));
    end
end

figure
set(gcf,'Position',[100 100 600 264.15])
contourf(x0',y0',sigma,80,'LineStyle','none')
axis([0 2 0 1])
colorbar
colormap(mycolormap)
clim([0,3])
set(gcf,'color','w')

%% Part 4 - Save several time series of sigma values for a video (runs longer)
clear all
A = 0.1;    % parameters from Shadden 2005 Physica D
eps = 0.25;
omega = 2*pi/10;  % frequency of gyre oscillations

dx = .01; % Increased grid resolution for finer ridges
xvec = 0:dx:2;
yvec = 0:dx:1;
[x0,y0] = meshgrid(xvec,yvec);  % grid of particles
yIC(1,:,:) = x0';
yIC(2,:,:) = y0';

dt =0.025;  % timestep
T = 15;     % duration of integration
int = 'f'; %''f for forward, 'b' for backward;
nSigma = 100; % number of FTLE snapshots

for r = 1:nSigma
    if int == 'f'
        sgn = 1;
        tVec = r:(T/dt)+(r-1); %time vector for forward integration
    else
        sgn = -1;
        tVec = flip(r:(T/dt)+(r-1)); %time vector for backwards integration
    end
    
    yin = yIC;
    for i=tVec
        time = i*dt;
        yout = rk4singlestep(@(t,y)doublegyreVEC(t,y,A,eps,omega),sgn*dt,sgn*time,yin);
        yin = yout;   
    end
    % reshape 3-dim array into 2-dim array (final positions)
    xT = reshape(yout(1,:,:),length(xvec),length(yvec));
    yT = reshape(yout(2,:,:),length(xvec),length(yvec));
    
    % compute flowmap gradient
    [dxTdx0,dxTdy0] = gradient(xT,dx,dx); 
    [dyTdx0,dyTdy0] = gradient(yT,dx,dx);
    
    %compute sigma values at each grid point for one time-step
    for i=1:length(xvec)
        for j=1:length(yvec)
            D(1,1) = dxTdx0(i,j);
            D(1,2) = dxTdy0(i,j);
            D(2,1) = dyTdx0(i,j);
            D(2,2) = dyTdy0(i,j);
            sigma(i,j,r) = (1/T)*sqrt(max(eig(D'*D)));
        end
    end
end

if int == 'f'
%save('ForwardFTLEgyreFlow.mat')
else
%save('BackwardFTLEgyreFlow.mat')
end
%% Part 5 - Video of Unsteady FTLE
% Finite difference to compute the gradient

if int == 'f'
    mycolormap = customcolormap(linspace(0,1,6),{'#a9212d','#b8412a','#ca6827','#edb121','#f5d586','#ffffff'});
else
    mycolormap = customcolormap(linspace(0,1,6),{'#132d8d','#0d449d','#0758ab','#2295e0','#a9d9f8','#ffffff'});
end
%compute sigma: large sigma indicates large mixing!
figure
set(gcf,'Position',[100 100 600 264.15])
for i=1:nSigma
    contourf(x0',y0',sigma(:,:,i),80,'LineStyle','none')
    axis([0 2 0 1])
    axis off
    colorbar
    colormap(mycolormap)
    clim([0,3])
    set(gca,'FontSize',14)
    set(gcf,'color','w')
    drawnow
end

%% Part 6 - Run a patch of particles + see attractiv/repulsive manifold

color = 'k';
if int == 'f'
    sgn = 1;
    mycolormap = customcolormap(linspace(0,1,6),{'#a9212d','#b8412a','#ca6827','#edb121','#f5d586','#ffffff'});
else
    sgn = -1;
    mycolormap = customcolormap(linspace(0,1,6),{'#132d8d','#0d449d','#0758ab','#2295e0','#a9d9f8','#ffffff'});
end

%initiate particle patch
dx = .025;
xvec = 1.0:dx:1.4; % location of particles from x=0 to 2
yvec = 0.3:dx:0.8; % location of particles from y=0 to 1
[x0t,y0t] = meshgrid(xvec,yvec);  % grid of particles
yICp(1,:,:) = x0t';
yICp(2,:,:) = y0t';
yin = yICp;

figure
set(gcf,'Position',[100 100 600 264.15])
for i=1:nSigma
    time = i*dt;     
    yout = rk4singlestep(@(t,y)doublegyreVEC(t,y,A,eps,omega),dt,time,yin);
    yin = yout;
    contourf(x0',y0',sigma(:,:,i),80,'LineStyle','none')
    hold on
    scatter(yout(1,:),yout(2,:),15,'o','Color',color,'MarkerEdgeColor','w','MarkerFaceColor','k')
    axis([0 2 0 1])
    set(gcf,'color','w')
    colormap(mycolormap)
    clim([0,2])
    xticklabels([])
    yticklabels([])
    hold off
    drawnow
end
