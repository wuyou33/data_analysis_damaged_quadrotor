function y = interpNan(x)

y = x;

nans = isnan(x);
vals = find(~nans);

nans([1:(min(vals)-1) (max(vals)+1):end]) = 0; % if first or last is NaN

if sum(nans)==0 
    y(isnan(x)) = 0;
else
    try
        %y(nans) = interp1(vals,x(vals),find(nans)); 
        y(nans) = interp1(vals,x(vals),find(nans),'line'); 
        y(isnan(y)) = 0; % but why is this happening??
    catch err
        display(err)
        return
    end
end

return
