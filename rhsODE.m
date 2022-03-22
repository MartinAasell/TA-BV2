function y = rhsODE(t,y)
% y = rhs=DE(t,y)
% Calculate the function f(t,y) = y-0.5e^(t/2)*sin(5t) + 5e^(t/2)*cos(5t)
y = y - 0.5*exp(t/2)*sin(5*t) + 5*exp(t/2)*cos(5*t);
end