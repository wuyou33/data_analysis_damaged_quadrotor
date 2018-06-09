
function zs = Frequency_Domain_Filter(z,Wc,dt)
% Frequency Domain Filter
% z: raw measurements
% Wc: Cut-off frequency
% dt: sampling period
% Sihao Sun 12.8.2016

N = length(z);
Ns = ceil(N/4);
Wt = Wc + 2*pi^2/(Ns*dt);
g = zeros(1,2*Ns+1);
t = zeros(1,2*Ns+1);
i = 1;
for k=-Ns:Ns
    g(i) = pi/(2*k*dt) * (sin(Wt*k*dt) + sin(Wc*k*dt))/(pi^2-(Wt-Wc)^2*(k*dt)^2);
    t(i) = k*dt;
    if isnan(g(i))
        g(i) = g(i-1);
    end
    i = i+1;
end
g = g/sum(g);

if isrow(z)
    zl = fliplr(z);
else
    zl = flipud(z);
end
zr = zl;
za = [zl,z,zr];
zs = z;
for i =1:N
    zs(i) = sum(g.*za(N+1-Ns+i:N+1+Ns+i));
end

end



