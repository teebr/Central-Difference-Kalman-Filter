ccc

%Numerical example based on the Van Der Pol example from R2016b UKF Demo

%Set up filter properties:
x0 = [5;-2];
Q = diag([0.02 , 0.1]);
R = 0.2;
P0 = 10*Q;
dt = 0.05;
mu = 1;

%generate true continuous data
t = (0:dt:10)';
[~,xTrue]=ode45(@(t,x) [x(2);-x(1) + mu*(1-x(1)^2)*x(2)],t,[2; 0]); %note the different ICs
yTrue = xTrue(:,1);
rng(1); % Fix the random number generator for reproducible results
yMeas = yTrue .* (1+sqrt(R)*randn(size(yTrue))); % sqrt(R): Standard deviation of noise

%run Simulink model with CDKF block
sim VDP_Sim_2014b
xp = squeeze(xp)';

%plots:
figure
hAx = subplot(2,1,1);
title(sprintf('Van Der Pol estimator with measurement\nnoise and incorrect initial states'))

plot(t,yTrue)
plot(t,yp)
plot(t,yMeas)
ylabel('x_1')
hLeg = legend({'True','Filter Estimate','Measured'},'Location','Best');
%toggle lines when the legend is clicked on:
ihFcn = @(h,e) set(e.Peer,'Visible',lower(regexprep(e.Peer.Visible,{'ff','n'},{'N','ff'},'once')));
hLeg.ItemHitFcn = ihFcn;

hAx(2,1) = subplot(2,1,2);
plot(t,xTrue(:,2))
plot(t,xp(:,2))
hAx(2).ColorOrderIndex = 1;
ylabel('x_2')
xlabel('Time [s]')
linkaxes(hAx,'x')