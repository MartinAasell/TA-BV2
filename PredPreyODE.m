function yout = PredPreyODE(t,y,alfa, beta, delta, gamma)
% yout = PredPreyODE(t,y,alfa,beta,delta,gamma);
% The Lotka-Volterra differential equation    
%    y'(t) = alfa*y1     - beta*y1*y2
%            delta*y1*y2 - gamma*y2
% yout is a vector: y=[y1; y2], y1 is the number of prey (i.e. rabbits), 
%                            y2 is the number of predator (i.e. foxes)
%
% alfa, beta, delta, gamma are positive real parameters describing the 
% interaction of the two species:
%   alfa  - Reproduction rate of prey
%   beta  - Mortality rate of predator per prey
%   delta - Reproduction rate of predator per prey
%   gamma - Mortality rate of predator
yout = [ alfa*y(1)       - beta*y(1)*y(2);
         delta*y(1)*y(2) - gamma*y(2)    ];