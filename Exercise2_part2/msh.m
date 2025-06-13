clear all
close all

global k Re

lwid = 2;
% axisFontSize = 12;
%% Parameters of the problem
% k = 0.5;
k = 1; % nu
Re = 10;

% problem solved: 
% d/dt(q1)  = [0 1  ](q1) 
%     (q2)    [k^2 0](q2)

%% Conditions at t = 0
q10 = 1;
q20 = 1;

q0 =[q10;
    q20];

%% Range in time
t0 = 0;
tf = 50;
%% Solve ODE using Runge-Kutta 4th - 5th order
[t, q] = ode45(@integra, [t0 tf], q0);

%% Plot solution
figure
plot(t, q(:, 1), 'k-', t, q(:, 2), 'b--');
xlabel('t')
ylabel('q_1(t), q_2(t)')
legend('q_1','q_2')