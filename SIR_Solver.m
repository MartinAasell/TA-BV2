t0 = 0;
t1 = 120;
timeinterval = [t0 t1];
beta = 0.3;
v = 1/7;
N = 1000;
mu = 0.002/365;
I0 = 5;
S0 = N - I0;
R0 = 0;

y0 = [S0; I0; R0];

setting = odeset('OutputFcn','odeplot','RelTol', 1e-2);
[t, y] = ode45(@(t,y) SIRmodel(t,y,N,mu,beta,v), timeinterval, y0, setting);
% Plot the solution if the plot-setting is not on (setting.OutputFcn is
% empty)
if isempty(setting.OutputFcn)
  plot(t,y);    
end  
  xlabel('time')
  ylabel('Number of predator and prey');
  title('Predator, prey simulation');
  legend('S','I', 'R');