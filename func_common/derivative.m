function [xdot] = derivative(x,t)
% [xdot] = derivative(x,t)
% This function differentiates the data vector x with respect to
% the specified time vector t using central differences
%

dt   = t(2,1) - t(1,1);
xdot = zeros(size(x));

for j = 1:size(x,2)
    %xdot(1,1) = 0; 
    for i = 2 : (length(t)-1)
        xdot(i,j) = 0.5 * (x(i+1,j) - x(i-1,j)) / (0.5*((t(i+1,1)-t(i-1,1))));
        if xdot(i,j) == Inf 
            xdot(i,j) = xdot(i-1,j);
        end
    end
end

xdot(length(t),:) = xdot((length(t)-1),:);
xdot(1,:) = xdot(2,:);% assume xdot(0)= xdot(1)
xdot = interpNan(xdot);

return