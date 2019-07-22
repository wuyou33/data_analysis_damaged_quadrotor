w = -1:0.1:1; v = 0.5;
A0 = [w', w'.^2, w'.*abs(v'), w'.^3, w'.^2.*abs(v'), w'.^2.*v'.^2];
A1 = [w', w'.^2, w'.*abs(v'), w'.^3, w'.^2.*abs(v'), w'.*v'.^2];
A2 = [sign(w').*w', sign(w').*w'.^2, sign(w').*w'.*abs(v'), sign(w').*w'.^3, sign(w').*w'.^2.*abs(v'), sign(w').*w'.*v'.^2];
A3 = [sign(w').*abs(w)', sign(w').*w'.^2, sign(w').*abs(w)'.*abs(v'), sign(w').*abs(w)'.^3,sign(w').* abs(w)'.^2.*abs(v'), sign(w').*abs(w)'.*v'.^2];
K = [-1.17, -4.07, -0.964, -0.284, 2.35, 3.41]';

Cz0 = A0*K;
Cz1 = A1*K;
Cz2 = A2*K;
Cz3 = A3*K;

figure
plot(w,Cz0); 
hold on;
plot(w,Cz1); 
plot(w,Cz2);
plot(w,Cz3);