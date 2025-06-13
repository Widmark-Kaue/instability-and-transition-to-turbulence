function integ = integra(~,x)
global k Re

% solve system of ODE's given by
% d/dt(q1)  = [0 1  ](q1) 
%     (q2)    [k^2 0](q2)

% A = [0 1;
%     -k^2 0];

A = [1/100 - 1/Re 0;
        k -2/Re];

integ = A*x;