
function [A_p52] = P52CtCq(aa,mu)
A_p52 = zeros(1,15);
A_p52 = [ones(size(mu)),mu, mu.^2, mu.^3, mu.^4, mu.^5, mu.*aa, mu.^2.*aa, mu.^3.*aa, mu.^4.*aa,...
         mu.*aa.^2, mu.^2.*aa.^2, mu.^3.*aa.^2,mu.*aa.^3, mu.^2.*aa.^3, mu.*aa.^4];
end