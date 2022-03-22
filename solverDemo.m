function [ ] = solverDemo(func, tspan, y0)
% solverDemo(@odefunc, tspan, y0)
% Demo program that illustrates the ODE-solvers Euler's method, Heuns
% method and Classical Runge-Kutta.
% The program solves the problem y'=func(t,y) on the interval tspan=[t0 t_n]
% The right-hand-side f(t,y) must be defined in a Matlab
% function, and return y corresponding to f(t,y).
%
% The solution is plotted step-by-step and is compared with the
% ode45-solution, which can be seen as 'exact' here. The step-by-step
% plotting get switched off if step size h is small leading to many
% intervals (it is inconvenient to plot many intervals step-by-step).
% Also, step-by-step plotting is switched off if 'All methods' is chosen.
% Only scalar problems allowed.
% Stefan Pålsson, 2018
maxinterval = 50; % Maximum number of intervals for plotting
if (length(y0) > 1)
    disp('Sorry, solverDemo do not accept systems of ODEs.');
    disp('Program interupted');
    return
end
% The accurat ('exact') solution, ode45
[texact, yexact] = ode45(func, tspan, y0);
h = OpenDialogBox('First, give stepsize h (h < 1)');
if isnan(h)
    disp('I can''t run without stepsize. Program interupted,');
    return
end
% Set h = interval length if h larger than interval
if h > (tspan(2)-tspan(1))
    h = tspan(2)-tspan(1);
end
if (h >= 1)
    disp('Stepsize h too big! Program interuppted.')
    return;
end
proceed = 1;
while proceed   
    choice = menu('Choose method/methods','Euler''s method','Heun''s method',...
        'Classical Runge-Kutta','All methods','Proceed but change stepsize h', ...
        'Stop program');
    method={};
    same_h = 1;
    switch choice
        case 1
            [t, y] = euler(func, tspan, y0, h);
            method{1} = 'Euler Method';
        case 2
            [t, y] = heun(func, tspan, y0, h);
            method{1} = 'Heuns method';
        case 3
            [t, y] = RK4(func, tspan, y0, h);
            method{1} = 'Classical Runge-Kutta';
        case 4
            [~, y_euler] = euler(func, tspan, y0, h);
            [~, y_heun] = heun(func, tspan, y0, h);
            [t, y_RK] = RK4(func, tspan, y0, h);
            y = [y_euler y_heun y_RK];
            method={'Euler Method' 'Heuns method' 'Classical Runge-Kutta'...
                'All methods'};
        case 5
            h = OpenDialogBox('You''d like to change stepsize h');
            if isnan(h)
                disp('OK, you did not change stepsize h.');
            end       
            if h > (tspan(2)-tspan(1))
                h = tspan(2)-tspan(1);
            end 
            same_h = 0;
        case 6
            proceed = 0;
            close all;
    end
    
    % Plot solution
    if proceed && same_h
        draw(t, y, texact, yexact, h, maxinterval, method)
    end
    
end
end   % End main
% -------------------
% Internal functions
% -------------------
% ----------------
% Euler method
%
function [t, yout] = euler(func, tspan, y0, h)
% [t,y] = euler(@f, tspan, y0, h)
% Using Eulers method to solve y'=f(t,y) on interval tspan with initial 
% value y0 and time step h. tspan = [t0 t1] when interval is t0 to t1. 
% f must be a function handle.
% Make y0 a colun vector if it is a row vector
y0 = y0(:); 
% Create time vector
t = (tspan(1):h:tspan(2))';
h_i = h*ones(length(t),1);
if t(end) ~= tspan(2)
    h_i = [h_i; tspan(2)-t(end)];
    t = [t; tspan(2)];   
end
yout = zeros(length(y0), length(t)); % Preassign space for efficiency
yout(:,1) = y0;   % y0 in first column
for i = 2:length(t)
    yout(:,i) = yout(:,i-1)+h_i(i)*(func(t(i-1),yout(:,i-1)));
end
% Flip yout
yout = yout';
end
% ----------------
% Heuns method
%
function [t, yout] = heun(func, tspan, y0, h)
% [t,y] = heun(@f, tspan, y0, h)
% Using Heuns method to solve y'=f(t,y) on interval tspan with initial 
% value y0 and time step h. tspan = [t0 t1] when interval is t0 to t1. 
% f must be a function handle.
% Make y0 a column vector if it is a row vector
y0 = y0(:); 
% Create time vector
t = (tspan(1):h:tspan(2))';
h_i = h*ones(length(t),1);
if t(end) ~= tspan(2)
    h_i = [h_i; tspan(2)-t(end)];
    t = [t; tspan(2)];
end
yout = zeros(length(y0), length(t)); % Preassign space for efficiency
yout(:,1) = y0;   % y0 in first column
for i = 2:length(t)
    k1 = func(t(i-1),yout(:,i-1));
    k2 = func(t(i-1)+h,yout(:,i-1)+h_i(i)*k1);
    yout(:,i) = yout(:,i-1)+0.5*h_i(i)*(k1+k2);
end
% Flip yout
yout = yout';
end
% ---------------------
% Classical Runge-Kutta
function [t, yout] = RK4(func, tspan, y0, h)
% [t,y] = RK4(@f, tspan, y0, h)
% Using Classical Runge Kutta to solve y'=f(t,y) on interval tspan with 
% initial value y0 and time step h. tspan = [t0 t1] when interval is t0 
% to t1. f must be a function handle.
% Make y0 a column vector if it is a row vector
y0 = y0(:); 
% Create time vector
t = (tspan(1):h:tspan(2))';
h_i = h*ones(length(t),1);
if t(end) ~= tspan(2)
    h_i = [h_i; tspan(2)-t(end)];
    t = [t; tspan(2)];
end
yout = zeros(length(y0),length(t)); % Preassign space for efficiency
yout(:,1) = y0;   % y0 in first column
for i = 2:length(t)
    k1 = func(t(i-1),yout(:,i-1));
    k2 = func(t(i-1)+0.5*h_i(i), yout(:,i-1)+0.5*h_i(i)*k1);
    k3 = func(t(i-1)+0.5*h_i(i), yout(:,i-1)+0.5*h_i(i)*k2);
    k4 = func(t(i-1)+h_i(i),yout(:,i-1)+k3*h_i(i));
    yout(:,i) = yout(:,i-1)+(h_i(i)/6)*(k1+ 2*k2 + 2*k3 + k4);
end
% Flip yout
yout = yout';
end 
% ---------------------
% Plotting
%
function [] = draw(t, y, texact, yexact, h, maxinterval, method)
len_method = length(method);
% Close same method figure windows if any, and open a new one
if isempty(findobj('Name',method{len_method}))== false
    close(method{len_method});
end
figh=figure('Name',method{len_method});
pos = get(figh,'position');
set(figh,'position',[pos(1:2)/2 pos(3:4)*1.2])
% Plot the ode45-solution
hh1 = plot(texact, yexact,'g--','linewidth',1);
title(['Stepzise = ',num2str(h),' => ',num2str(length(t)),' intervals']);
xlabel('t'); ylabel('y')
axlar = axis;
legendtext(1) = {'Exact solution (ode45-solution)'};
xlabel('t'); ylabel('y');
hold on;
n = length(t);
if len_method == 1  % If it's just one method, plot step by step
    % Switch off step by step plotting due to too many intervals
    if n > maxinterval
        rita = 0;
        str = ['Wont plot step by step! Stepsize=',num2str(h),...
            ' gives too many intervals'];
    else
        rita = 1;
        str=['Push any key on keyboard to continue with a new ',...
            method{1},' step'];
    end
    textposx = axlar(1)+(axlar(2)-axlar(1))/30;
    textposy = axlar(3)+(axlar(4)-axlar(3))/10;
    h_text = text(textposx, textposy,str, 'FontSize',12,'Color', 'red');
    disp('--');
    disp(str);
    
    % Draw solution step by step
    if rita
        axis([axlar(1) axlar(2) min([min(y) axlar(3)]) max([max(y) axlar(4)])]);
        axis manual
        pause on
        plot(t(1),y(1),'r*');
        for i = 2:n
            if i<=n, pause; end
            line([t(i-1) t(i)], [y(i-1) y(i)],'linewidth',1.5);
            plot(t(i),y(i),'r*');
        end
        delete(h_text);
    end
end
hh2=plot(t,y,'linewidth',1.5);
%str = [ method{len_method} ', stepsize ',num2str(h)];
% Change legend depending on one method or all methods (len_method==3)
if len_method == 4
    legendtext = [legendtext method{1:3}];
    handlevect = [];
else
    legendtext = [legendtext method];
    handlevect =[hh1 hh2];
end
legend(handlevect, legendtext, 'FontSize',12);
hold off
end
function h = OpenDialogBox(txt1)
% h = OpenDialogBox(txt)
% Dialogbox for choosing stepsize h
% The text in the dialog box
txt2 = 'h=';
txt = {sprintf([txt1 '\n' txt2])};      % Put dialog text together
%txt = {sprintf([txt1 '\n' txt2 '\n\n' txt3])};
dlg_title = 'Stepsize h'; % Titel in dialog box
num_lines = [1 50];       % Size of the dialog bow
%options.Interpreter = 'latex';
options.Interpreter = 'none';
options.WindowStyle = 'normal';
options.Resize = 'off';
answer = inputdlg(txt, dlg_title, num_lines, {''}, options);
if isempty(answer)
    h = NaN;
    return
end
h = str2double(answer{:});
end % function OpenDialogBox
