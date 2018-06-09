function [x_integ] = integral_DA2(x,t,x_int_0,method)
% [x_integ] = integral_DA2(x,t,x_int_0,method)
% 
% This function integrates the vector x over the time vector t, using
% Riemann sums
% 
% input
% x       : vector to integrate
% t       : time vector to integrate over
% x_int_0 : (optional) initial value of integral (if known)
% method  : (optional) integration method, 'euler' or 'rk'
%

x_integ = zeros(size(x,1),1);
if nargin < 3
    x_int_0 = 0;
end
x_integ(1,1) = x_int_0;

if (strcmp((method),'euler') || nargin < 4)
    for i = 2 : size(x,1) 
        h = t(i,1) - t(i-1,1);
        x_integ(i,1) = x_integ(i-1,1) + ((x(i) + x(i-1)) / 2) * h;
    end

% option 2 : RK
elseif strcmp((method),'rk')
    for i = 2 : size(x,1) 
        h = t(i,1) - t(i-1,1);
        k1 = h*((x(i) + x(i-1)) / 2);
        k2 = h*((x(i) + x(i-1) + 2*h/2) / 2);
        k3 = k2;
        k4 = h*((x(i) + x(i-1) + 2*h) / 2);
        x_integ(i,1) = x_integ(i-1,1) + k1/6 + k2/3 + k3/3 + k4/6;
        %((x(i) + x(i-1)) / 2) * h;
    end
% k1 = hf(xn; yn)
% k2 = hf(xn +h/2; yn + k1/2);
% k3 = hf(xn + h/2; yn + k2/2);
% k4 = hf(xn + h; yn + k3);
% yn+1 = yn + k1/6 + k2/3 + k3/3 + k4/6 + O(h5);

% [tt,yy] = ode23('testFct',TIME,phidot);

% % option 2 : RK
% elseif strcmp((method),'rkPr')
%     h = t(i,1) - t(i-1,1);
%     t_new = t(1,1):h/5:max(t);
%     h_new = h/5;
%     for i = 2 : size(x,1) 
%         k1 = h_new*((x(i) + x(i-1)) / 2);
%         k2 = h*((x(i) + x(i-1) + 2*h/2) / 2);
%         k3 = k2;
%         k4 = h*((x(i) + x(i-1) + 2*h) / 2);
%         x_integ(i,1) = x_integ(i-1,1) + k1/6 + k2/3 + k3/3 + k4/6;
%         %((x(i) + x(i-1)) / 2) * h;
%     end
end

return