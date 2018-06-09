function [p,q,r] = Eulerdot2pqr_vector(phi0,theta0,psi0,time)


% phi = Frequency_Domain_Filter(phi0,500,dt);
phi = smooth(phi0,50);
% phit = [0;diff(phi)/dt];
phit0 = derivative(phi,time);
phit = smooth(phit0,50);

% theta = Frequency_Domain_Filter(theta0,50,dt);
theta = smooth(theta0,50);
% thetat = [0;diff(theta)/dt];
thetat0 = derivative(theta,time);
thetat = smooth(thetat0,50);

% psi = Frequency_Domain_Filter(psi0,50,dt);
psi = smooth(psi0,50);
% psit = [0;diff(psi)/dt];
psit0 = derivative(psi,time);
psit = smooth(psit0,50);

p = phit - sin(theta).*psit;
q = cos(phi).*thetat + sin(phi).*cos(theta).*thetat;
r = -sin(phi).*thetat + cos(phi).*cos(theta).*psit;

end