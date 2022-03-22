
% PredPreySim
% Script simulating the interaction between predator and prey, described
% by the Lotka-Volterra ordinary differential equation.
% The Lotka-Volterra ODE is defined in the Matlab function PredPreyODE,
% which must be available in the same folder as this script.
% See PredPreyODE help text for more information.
% Time interval for the simulation
t0 = 0;
t1 = 40;
timeinterval = [t0 t1];
% Initial conditions, e.g. number of predator and prey at time zero
pred0 = 500;
prey0 = 1000;
y0 = [prey0; pred0]; % Vector y0. Note, the order is critical. 
% Set parameters alfa, beta, delta, gamma
alfa = 5;   % Reproduction rate of prey
beta = 0.02;  % Mortality rate of predator per prey
delta = 0.02; % Reproduction rate of predator per prey
gamma = 0.4;  % Mortality rate of predator
% Here you can do different settings for Matlabs ODE solvers (see help text
% help odeset for details). 
% Note, this is not necessary and only used when you'd like to change the 
% default settings. Normally you just leave it out.
% The code 'setting=odeset' result in no change.
% setting = odeset;
% You can try other setting though, like
 setting = odeset('OutputFcn','odeplot','RelTol', 1e-2);
% The relative tolerance will be changed to 10^-2 and the ode-solver will
% plot the solution automatically
% Call ode45. The parametern setting is optional. Used only when changing default 
% settings. 
[t, y] = ode45(@(t,y) PredPreyODE(t,y,alfa,beta,delta,gamma), timeinterval, y0, setting);
% Plot the solution if the plot-setting is not on (setting.OutputFcn is
% empty)
if isempty(setting.OutputFcn)
  plot(t,y);    
end  
  xlabel('time')
  ylabel('Number of predator and prey');
  title('Predator, prey simulation');
  legend('Number of prey','Number of predator');
