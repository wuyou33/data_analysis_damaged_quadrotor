%%
clear DATA_raw;

addpath('data')
addpath(genpath('E:\Data\_code_import_files'));

% DATA_raw{1} = load('17_23480-86200.mat');
% DATA_raw{2} = load('18_7240-88380.mat');
% DATA_raw{3} = load('19_27950-88220.mat');
% DATA_raw{1} = load('20_23850-87430.mat');
% DATA_raw{1} = load('75_35590-119800.mat');
% DATA_raw{1} = load('20_79690-87430.mat');
DATA_raw{1} = load('18_81820-88380.mat');
% DATA_raw{1} = load('75_114700-119800.mat');
% DATA_raw{1} = load('18_62130-69970.mat');
% DATA_raw{1} = load('9_25120-91680.mat');
% DATA_raw{1} = load('7_36850-96460.mat');
% DATA_raw{2} = load('32_29110-84740.mat');

Nc = 2;

DATA = DATA_raw;
names = fieldnames(DATA_raw{1});
for i = 1:length(DATA_raw)
    for j = 1:length(names)
        if ~strcmp(names{j},'DU')
            DATA{i}.(names{j}) = compress_matrix(DATA_raw{i}.(names{j})(DATA_raw{i}.DU(1):DATA_raw{i}.DU(2)),Nc);
        end
    end
end
fx = []; fy = []; fz = [];
mx = []; my = []; mz = [];
Fx = []; Fy = []; Fz = [];
Mx = []; My = []; Mz = [];
u  = []; v  = []; w  = [];
u1 = []; u2 = []; u3 = []; u4 = [];
v1 = []; v2 = []; v3 = []; v4 = [];
w1 = []; w2 = []; w3 = []; w4 = [];
omega1_bar = []; omega2_bar = []; omega3_bar = []; omega4_bar = [];
omega1 = []; omega2 = []; omega3 = []; omega4 = [];
alpha = [];  psi_mod = [];
for i = 1:length(DATA)
   fx = [fx ;DATA{i}.fx];
   fy = [fy ;DATA{i}.fy];
   fz = [fz ;DATA{i}.fz];
   Fx = [Fx ;DATA{i}.Fx];
   Fy = [Fy ;DATA{i}.Fy];
   Fz = [Fz ;DATA{i}.Fz];
   mx = [mx ;DATA{i}.mx];
   my = [my ;DATA{i}.my];
   mz = [mz ;DATA{i}.mz];   
   Mx = [Mx ;DATA{i}.Mx];
   My = [My ;DATA{i}.My];
   Mz = [Mz ;DATA{i}.Mz];
   u = [u ;DATA{i}.u];
   v = [v ;DATA{i}.v];
   w = [w ;DATA{i}.w];
   u1 = [u1 ;DATA{i}.u1];
   u2 = [u2 ;DATA{i}.u2];
   u3 = [u3 ;DATA{i}.u3];
   u4 = [u4 ;DATA{i}.u4];
   v1 = [v1 ;DATA{i}.v1];
   v2 = [v2 ;DATA{i}.v2];
   v3 = [v3 ;DATA{i}.v3];
   v4 = [v4 ;DATA{i}.v4];
   w1 = [w1 ;DATA{i}.w1];
   w2 = [w2 ;DATA{i}.w2];
   w3 = [w3 ;DATA{i}.w3];
   w4 = [w4 ;DATA{i}.w4];   
   omega1_bar = [omega1_bar ;DATA{i}.omega1_bar];
   omega2_bar = [omega2_bar ;DATA{i}.omega2_bar];
   omega3_bar = [omega3_bar ;DATA{i}.omega3_bar];
   omega4_bar = [omega4_bar ;DATA{i}.omega4_bar];
   omega1 = [omega1 ;DATA{i}.omega1];
   omega2 = [omega2 ;DATA{i}.omega2];
   omega3 = [omega3 ;DATA{i}.omega3];
   omega4 = [omega4 ;DATA{i}.omega4];
   alpha = [alpha ;DATA{i}.alpha];
   psi_mod = [psi_mod ;DATA{i}.psi_mod];
end
du = 1:length(fx);
%%
h = 30; %order of the Fourier series

psi_h1 = psi_mod(du) - 417; psi_h1 = psi_h1/57.3;
psi_h2 = psi_mod(du) - 308; psi_h2 = psi_h2/57.3;
psi_h3 = psi_mod(du) - 232; psi_h3 = psi_h3/57.3;
psi_h4 = psi_mod(du) - 128; psi_h4 = psi_h4/57.3;

psi_h1 = mod(psi_h1,2*pi);
psi_h2 = mod(psi_h2,2*pi);
psi_h3 = mod(psi_h3,2*pi);
psi_h4 = mod(psi_h4,2*pi);

%% batch model

load('Bebop2_guido_parameters.mat');
b = parameters.b;
l = parameters.l;
SL = [1 -1 -1 1];
SM = [1 1 -1 -1];
SN = [-1 1 -1 1];

AFX1 = [ones(length(du),1).*(u1.*omega1_bar+u2.*omega2_bar+u3.*omega3_bar+u4.*omega4_bar)];
AFX2 = [Fn(psi_h1,h,0,u1.*omega1_bar) + Fn(2*pi-psi_h2,h,0,u2.*omega2_bar)];
AFX3 = [Fn(psi_h3,h,0,u3.*omega3_bar) + Fn(2*pi-psi_h4,h,0,u4.*omega4_bar)];
AFY4 = [ones(length(du),1).*(v1.*omega1_bar+v2.*omega2_bar+v3.*omega3_bar+v4.*omega4_bar)];
AFY5 = [Fn(psi_h1,h,0,v1.*omega1_bar) + Fn(2*pi-psi_h2,h,0,v2.*omega2_bar)];
AFY6 = [Fn(psi_h3,h,0,v3.*omega3_bar) + Fn(2*pi-psi_h4,h,0,v4.*omega4_bar)];
AN1  = [ones(length(du),1).*(b*SL(1)*u1.*omega1_bar+SL(2)*b*u2.*omega2_bar+SL(3)*b*u3.*omega3_bar+SL(4)*b*u4.*omega4_bar)];
AN2  = [b*SL(1)*Fn(psi_h1,h,0,u1.*omega1_bar) + b*SL(2)*Fn(2*pi-psi_h2,h,0,u2.*omega2_bar)];
AN3  = [b*SL(3)*Fn(psi_h3,h,0,u3.*omega3_bar) + b*SL(4)*Fn(2*pi-psi_h4,h,0,u4.*omega4_bar)];
AN4  = [ones(length(du),1).*(l*SM(1)*v1.*omega1_bar+l*SM(2)*v2.*omega2_bar+l*SM(3)*v3.*omega3_bar+l*SM(4)*v4.*omega4_bar)];
AN5  = [l*SM(1)*Fn(psi_h1,h,0,v1.*omega1_bar) + l*SM(2)*Fn(2*pi-psi_h2,h,0,v2.*omega2_bar)];
AN6  = [l*SM(3)*Fn(psi_h3,h,0,v3.*omega3_bar) + l*SM(4)*Fn(2*pi-psi_h4,h,0,v4.*omega4_bar)];
AN7  = [SN(1)*P33(u1,v1,w1,omega1_bar) + SN(2)*P33(u2,v2,w2,omega2_bar) + SN(3)*P33(u3,v3,w3,omega3_bar) + SN(4)*P33(u4,v4,w4,omega4_bar)];

% AN7  = [P33(u,v,w)];

A_est = [AFX1    , AFX2    , AFX3    , zz(AFY4), zz(AFY5), zz(AFY6) , zz(AN7) ;
         zz(AFX1), zz(AFX2), zz(AFX3), AFY4    , AFY5    , AFY6    , zz(AN7);
         AN1     , AN2     , AN3     , AN4     , AN5     , AN6     , AN7];
% A_est = [AFX1    , AFX2    , AFX3    , zz(AFY4), zz(AFY5), zz(AFY6) ;
%          zz(AFX1), zz(AFX2), zz(AFX3), AFY4    , AFY5    , AFY6    ;
%          AN1     , AN2     , AN3     , AN4     , AN5     , AN6     ];

Z_est = [Fx(du); Fy(du); Mz(du)-mz(du)];
sigma = ones(size(Z_est));
sigma(1:length(du)*2) = 1e-3;
W = spdiags(sigma,0,length(sigma),length(sigma));
K = (A_est'*W*A_est)\A_est'*W*Z_est;
Y = A_est*K;
y = [fx;fy;mz];

id_k = cumsum([1,size(AFX1,2),size(AFX2,2),size(AFX3,2),size(AFY4,2),size(AFY5,2),size(AFY6,2),size(AN7,2)]);

k = cell(1,6);
for i =1:length(k)
   k{i} = K(id_k(i):id_k(i+1)-1); 
end

R2 = find_R2(Y,Z_est);
RMS = find_RMS(Y,Z_est);

fprintf('R2:\t%f\nRMS:\t%f\n\n',R2,RMS);

CC = corr(A_est);
[row,col] = find(abs(CC)>0.9);

% figure
% plot(row,col,'.');

figure
plot(Z_est(1:2*length(du))); hold on; 
plot(Y(1:2*length(du)));

figure
plot(Mz(du)-mz(du)); hold on;
plot(Y(2*length(du)+1:end));
%% submodel check without alpha batch
pp = 0:0.01:2*pi; pp = pp';
yx12 = Fn(pp,h,0,1)*k{2} + k{1};
yx34 = Fn(pp,h,0,1)*k{3} + k{1};

yy12 = Fn(pp,h,0,1)*k{5} + k{4};
yy34 = Fn(pp,h,0,1)*k{6} + k{4};

figure
subplot(2,1,1)
plot(pp*57.3,yx12); hold on;
plot(pp*57.3,yx34);
subplot(2,1,2)
plot(pp*57.3,yy12); hold on;
plot(pp*57.3,yy34);
return;
%%
% % lambda(psi) model
% Ax = [ones(size(u1(du))).*(u1(du).*omega1_bar(du) + u3(du).*omega3_bar(du)),...
%       P1n(psi_h1,6,1,u1(du).*omega1_bar(du)),...
%       P1n(psi_h3,6,1,u3(du).*omega3_bar(du))];
% Ay = [ones(size(v1(du))).*(v1(du).*omega1_bar(du) + v3(du).*omega3_bar(du)),...
%       P1n(psi_h1,6,1,v1(du).*omega1_bar(du)),...
%       P1n(psi_h3,6,1,v3(du).*omega3_bar(du))];

%periodic model
Ax = [ones(length(du),1).*(u1(du).*omega1_bar(du)+u2(du).*omega2_bar(du)+u3(du).*omega3_bar(du)+u4(du).*omega4_bar(du)),...
      Fn(psi_h1,h,0,u1(du).*omega1_bar(du)) + Fn(2*pi-psi_h2,h,0,u2(du).*omega2_bar(du)),...
      Fn(psi_h3,h,0,u3(du).*omega3_bar(du)) + Fn(2*pi-psi_h4,h,0,u4(du).*omega4_bar(du))];
        
Ay = [ones(length(du),1).*(v1(du).*omega1_bar(du)+v2(du).*omega2_bar(du)+v3(du).*omega3_bar(du)+v4(du).*omega4_bar(du)),...
      Fn(psi_h1,h,0,v1(du).*omega1_bar(du)) + Fn(2*pi-psi_h2,h,0,v2(du).*omega2_bar(du)),...
      Fn(psi_h3,h,0,v3(du).*omega3_bar(du)) + Fn(2*pi-psi_h4,h,0,v4(du).*omega4_bar(du))];


 
%periodic model with alpha
% Ax = [ones(length(du),1).*(u2(du).*omega2_bar(du)+u4(du).*omega4_bar(du)),...
%       Fn(psi_h2,h,0,u2(du).*omega2_bar(du),alpha(du)),...
%       Fn(psi_h4,h,0,u4(du).*omega4_bar(du),alpha(du))];
% 
% Ay = [ones(length(du),1).*(v2(du).*omega2_bar(du)+v4(du).*omega4_bar(du)),...
%       Fn(psi_h2,h,0,v2(du).*omega2_bar(du)),...
%       Fn(psi_h4,h,0,v4(du).*omega4_bar(du))];  
%   
zx = Fx(du);
zy = Fy(du);
% 
[kx] = lscov(Ax,zx);
[ky] = lscov(Ay,zy);
yx = Ax*kx;
yy = Ay*ky;

figure
plot(psi_mod(du),Fx(du),'.'); hold on;
plot(psi_mod(du),fx(du),'.');
plot(psi_mod(du),yx,'.');

% figure
% plot(psi_mod(du),Fy(du),'.'); hold on;
% plot(psi_mod(du),fy(du),'.');
% plot(psi_mod(du),yy,'.');

R2x = find_R2(yx,Fx(du));
RMSx = find_RMS(yx,Fx(du));

R2y = find_R2(yy,Fy(du));
RMSy = find_RMS(yy,Fy(du));

fprintf('R2:\t%f\t %f \nRMS:\t%f\t %f ',R2x,R2y, RMSx, RMSy);
%% submodel check with alpha
pp = 0:0.1:2*pi; pp = pp';
aa = -20:20; aa = aa';

[pp_grid,aa_grid] = meshgrid(pp,aa);

mm = (size(Ax,2)-1)/2;
kx1 = kx(2:mm+1); kx3 = kx(mm+2:end);
% ky1 = ky(2:mm+1); ky3 = ky(mm+2:end);

yx1 = zeros(length(aa),length(pp));
yx3 = zeros(length(aa),length(pp));

for i = 1:length(aa)
yx1(i,:) = Fn(pp,h,0,1,aa(i))*kx1;
yx3(i,:) = Fn(pp,h,0,1,aa(i))*kx3;
end


% figure
% plot(pp*57.3,A(:,1:mm)*kx1); hold on;
% plot(pp*57.3,A(:,mm+1:end)*kx3);
% 
% figure
% plot(pp*57.3,A(:,1:mm)*ky1); hold on;
% plot(pp*57.3,A(:,mm+1:end)*ky3);

figure
mesh(aa_grid,pp_grid,yx1);

figure
mesh(aa_grid,pp_grid,yx3);

%% submodel check without alpha
pp = 0:0.01:2*pi; pp = pp';

mm = (size(Ax,2)-1)/2;
kx1 = kx(2:mm+1); kx3 = kx(mm+2:end);
ky1 = ky(2:mm+1); ky3 = ky(mm+2:end);

yx1 = Fn(pp,h,0,1)*kx1 + kx(1);
yx3 = Fn(pp,h,0,1)*kx3 + kx(1);
yy1 = Fn(pp,h,0,1)*ky1 + ky(1);
yy3 = Fn(pp,h,0,1)*ky3 + ky(1);

% yx1 = Fn(pp,h,0,1)*kx(2:end);
figure
plot(pp*57.3,yx1); hold on;
plot(pp*57.3,yx3);

figure
plot(pp*57.3,yy1); hold on;
plot(pp*57.3,yy3);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
psi14 = psi_mod(du)/57.3;
psi14(psi14<pi) = 0; psi14(psi14>pi) = sin(psi14(psi14>pi) - pi);

psi23 = psi_mod(du)/57.3;
psi23(psi23>pi) = 0; psi23(psi23<pi) = sin(psi23);


Ax = [ones(size(u1(du))).*(u1(du).*omega1_bar(du) + u3(du).*omega3_bar(du)),...
      P1n(psi14,4,0,u1(du).*omega1_bar(du)),...
      P1n(psi23,4,0,u3(du).*omega3_bar(du))];
% Ax = [ones(size(u1(du))).*(u1(du).*omega1_bar(du) + u3(du).*omega3_bar(du)),...
%     (psi_mod(du)).*u1(du).*omega1_bar(du),(psi_mod(du)).^2.*u1(du).*omega1_bar(du),(psi_mod(du)).^3.*u1(du).*omega1_bar(du), (psi_mod(du)).^4.*u1(du).*omega1_bar(du),...
%     (psi_mod(du)).*u3(du).*omega3_bar(du),(psi_mod(du)).^2.*u3(du).*omega3_bar(du),(psi_mod(du)).^3.*u3(du).*omega3_bar(du),(psi_mod(du)).^4.*u3(du).*omega3_bar(du)];
% 
Ay = [ones(size(v1(du))).*(v1(du).*omega1_bar(du) + v3(du).*omega3_bar(du)),...
    (psi_mod(du)).*v1(du).*omega1_bar(du),(psi_mod(du)).^2.*v1(du).*omega1_bar(du),(psi_mod(du)).^3.*v1(du).*omega1_bar(du), (psi_mod(du)).^4.*v1(du).*omega1_bar(du),...
    (psi_mod(du)).*v3(du).*omega3_bar(du),(psi_mod(du)).^2.*v3(du).*omega3_bar(du),(psi_mod(du)).^3.*v3(du).*omega3_bar(du),(psi_mod(du)).^4.*v3(du).*omega3_bar(du)];
