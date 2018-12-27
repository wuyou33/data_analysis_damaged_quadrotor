%%
clear DATA_raw;

addpath('data')
addpath(genpath('E:\Data\_code_import_files'));

%DRF r<0
DATA_raw{1} = load('9_25120-91680.mat');
DATA_raw{2} = load('17_23480-86200.mat');
DATA_raw{3} = load('18_7240-88380.mat');
%DRF r>0
DATA_raw{4} = load('7_36850-96460.mat');
DATA_raw{5} = load('19_27950-88220.mat');
DATA_raw{6} = load('20_23850-87430.mat');
%SRF
DATA_raw{7} = load('32_29110-84740.mat'); %r<0 lf
DATA_raw{8} = load('23_30700-100200.mat'); %r>0 lb 
DATA_raw{9} = load('25_31500-92750.mat'); %r<0 rb
DATA_raw{10} = load('35_22600-53420.mat'); %r>0 rf

%Norminal
DATA_raw{11} = load('92_24380-90710.mat'); % psi = 0;
DATA_raw{12} = load('93_1-63210.mat'); % psi = -90;

% DATA_raw{7} = load('75_35590-119800.mat');
% DATA_raw{1} = load('20_79690-87430.mat');
% DATA_raw{1} = load('18_81820-88380.mat');
% DATA_raw{1} = load('75_114700-119800.mat');
% DATA_raw{1} = load('18_62130-69970.mat');
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
h = 5; %order of the Fourier series

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
R = parameters.R;
Area = pi*R^2;
rho = 1.225;
S = 4*b*l;

SL = [1 -1 -1 1];
SM = [1 1 -1 -1];
SN = [-1 1 -1 1];

dynhead1 = rho.*omega1.^2.*R^2;
dynhead2 = rho.*omega2.^2.*R^2;
dynhead3 = rho.*omega3.^2.*R^2;
dynhead4 = rho.*omega4.^2.*R^2;

mu1 = sqrt(u1.^2+v1.^2)./(omega1.*R); mu1(abs(mu1)==inf)=0;
mu2 = sqrt(u2.^2+v2.^2)./(omega2.*R); mu2(abs(mu2)==inf)=0;
mu3 = sqrt(u3.^2+v3.^2)./(omega3.*R); mu3(abs(mu3)==inf)=0;
mu4 = sqrt(u4.^2+v4.^2)./(omega4.*R); mu4(abs(mu4)==inf)=0;

mux1 = u1./(omega1.*R); mux1(abs(mux1)==inf)=0;
mux2 = u2./(omega2.*R); mux2(abs(mux2)==inf)=0;
mux3 = u3./(omega3.*R); mux3(abs(mux3)==inf)=0;
mux4 = u4./(omega4.*R); mux4(abs(mux4)==inf)=0;

muy1 = v1./(omega1.*R); muy1(abs(muy1)==inf)=0;
muy2 = v2./(omega2.*R); muy2(abs(muy2)==inf)=0;
muy3 = v3./(omega3.*R); muy3(abs(muy3)==inf)=0;
muy4 = v4./(omega4.*R); muy4(abs(muy4)==inf)=0;

lc1 = w1./(omega1.*R); lc1(abs(lc1)==inf)=0;
lc2 = w2./(omega2.*R); lc2(abs(lc2)==inf)=0; 
lc3 = w3./(omega3.*R); lc3(abs(lc3)==inf)=0;
lc4 = w4./(omega4.*R); lc4(abs(lc4)==inf)=0;

beta= calc_beta(u,v);
beta_abs = abs(atan(v./u));
alpha = atan(w./sqrt(u.^2+v.^2+w.^2));
% alpha = calc_alpha(sqrt(u.^2+v.^2),w);
va = sqrt(u.^2+v.^2+w.^2);

% AT0 = [P32(alpha,abs(beta),va.^2*rho*S)];
AT0 = [P33(u,abs(v),w)];
AT1 = [Fn(psi_h1,h,0,dynhead1*Area,mu1,lc1) + Fn(2*pi-psi_h2,h,0,dynhead2*Area,mu2,lc2) ...
         + Fn(psi_h3,h,0,dynhead3*Area,mu3,lc3) + Fn(2*pi-psi_h4,h,0,dynhead4*Area,mu4,lc4)];
AL1 = [Fn(psi_h1,h,0,dynhead1*Area*b*SL(1),mu1,lc1) + Fn(2*pi-psi_h2,h,0,dynhead2*Area*b*SL(2),mu2,lc2) ...
         + Fn(psi_h3,h,0,dynhead3*Area*b*SL(3),mu3,lc3) + Fn(2*pi-psi_h4,h,0,dynhead4*Area*b*SL(4),mu4,lc4)];     
AM1 = [Fn(psi_h1,h,0,dynhead1*Area*l*SM(1),mu1,lc1) + Fn(2*pi-psi_h2,h,0,dynhead2*Area*l*SM(2),mu2,lc2) ...
         + Fn(psi_h3,h,0,dynhead3*Area*l*SM(3),mu3,lc3) + Fn(2*pi-psi_h4,h,0,dynhead4*Area*l*SM(4),mu4,lc4)];     
AT1fr = [Fn(psi_h1,h,0,dynhead1*Area,mu1,lc1) + Fn(2*pi-psi_h2,h,0,dynhead2*Area,mu2,lc2), ...
         + Fn(psi_h3,h,0,dynhead3*Area,mu3,lc3) + Fn(2*pi-psi_h4,h,0,dynhead4*Area,mu4,lc4)];
AL1fr = [Fn(psi_h1,h,0,dynhead1*Area*b*SL(1),mu1,lc1) + Fn(2*pi-psi_h2,h,0,dynhead2*Area*b*SL(2),mu2,lc2), ...
         + Fn(psi_h3,h,0,dynhead3*Area*b*SL(3),mu3,lc3) + Fn(2*pi-psi_h4,h,0,dynhead4*Area*b*SL(4),mu4,lc4)];     
AM1fr = [Fn(psi_h1,h,0,dynhead1*Area*l*SM(1),mu1,lc1) + Fn(2*pi-psi_h2,h,0,dynhead2*Area*l*SM(2),mu2,lc2), ...
         + Fn(psi_h3,h,0,dynhead3*Area*l*SM(3),mu3,lc3) + Fn(2*pi-psi_h4,h,0,dynhead4*Area*l*SM(4),mu4,lc4)];
AL2 = [ -muy1 - muy2 - muy3 - muy4];
AL3 = [-SN(1)*mux1 - SN(2)*mux2 - SN(3)*mux3 - SN(4)*mux4];
AM2 = [ mux1 + mux2 + mux3 + mux4];
AM3 = [-SN(1)*muy1 - SN(2)*muy2 - SN(3)*muy3 - SN(4)*muy4];

% AL4 = [P32(alpha,sign(v).*beta_abs,va.^2*rho*S)];
% AM5 = [P32(sign(u).*alpha,beta_abs,va.^2*rho*S)];
AL4 = P33(u,abs(v),w,sign(v));
AM5 = P33(u,abs(v),w);

% frame T ; rotor dT dM 0.0424 kc 0.2366 ks 0.0634
% AA =    [{AT0}     {AT1} {zz(AL2)} {zz(AL3)};
%          {zz(AT0)} {AL1} {AL2}     {AL3}   ;
%          {zz(AT0)} {AM1} {AM2}     {AM3}   ]; 
     
% frame T + M ; rotor dT 0.037729
% AA =    [{AT0}     {AT1} {zz(AL4)} {zz(AM5)};
%          {zz(AT0)} {AL1} {AL4}     {zz(AM5)};
%          {zz(AT0)} {AM1} {zz(AL4)} {AM5}]; 
% % %      
% % % frame T + M ; rotor dT dM 0.037521 kc -0.0355 ks 0.0511
AA =    [{AT0}     {AT1} {zz(AL2)} {zz(AL3)} {zz(AL4)} {zz(AM5)};
         {zz(AT0)} {AL1} {AL2}     {AL3}     {AL4}     {zz(AM5)};
         {zz(AT0)} {AM1} {AM2}     {AM3}     {zz(AL4)} {AM5}]; 
     


A_est = cell2mat(AA);

Z_est = [-Fz(du)+fz(du);
         Mx(du)-mx(du);
         My(du)-my(du)];

sigma = ones(size(Z_est));
sigma(1:length(du)) = 5e-3;
W = spdiags(sigma,0,length(sigma),length(sigma));     
     
K = (A_est'*W*A_est)\A_est'*W*Z_est;
Y =  A_est*K;

Nk = size(AA,2);
ik = zeros(Nk,1); k_model = cell(Nk,1);
for i = 1:Nk
    if i == 1
        ik(i) = size(AA{1,i},2);
        k_model{i} = K(1:ik(i));
    else
        ik(i) = size(AA{1,i},2)+ik(i-1);
        k_model{i} = K(ik(i-1)+1:ik(i));
    end
end

R2 = find_R2(Y,Z_est);
RMS = find_RMS(Y,Z_est);

fprintf('R2:\t%f\nRMS:\t%f\n\n',R2,RMS);

CC = corr(A_est);
[row,col] = find(abs(CC)>0.90);
% 
figure
plot(row,col,'.');

figure
plot(Z_est); hold on; 
plot(Y);

Y_fz = 0; Y_mx = 0; Y_my = 0;
for i = 1:Nk
    Y_fz = Y_fz + AA{1,i}*k_model{i};
    Y_mx = Y_mx + AA{2,i}*k_model{i};
    Y_my = Y_my + AA{3,i}*k_model{i};
end
figure
subplot(3,1,1)
plot(-Fz(du)); hold on;
plot(Y_fz-fz(du)); 
plot(-fz(du),'linestyle',':');
subplot(3,1,2)
plot(Mx(du)); hold on;
plot(mx(du)+Y_mx); plot(mx(du),'linestyle','-.'); 
subplot(3,1,3)
plot(My(du)); hold on;
plot(my(du)+Y_my); plot(my(du),'linestyle','-.'); 
% save('E:\system identification\thrust_model\dCt_model_BB2','K','h');
% return;
%% submodel check same characteristics
pp = linspace(0,2*pi,100)';
mmu = linspace(0,0.17,30)'; llc = linspace(-0.07,0.03,30)';

[pp_grid,llc_grid,mmu_grid] = meshgrid(pp,llc,mmu);

dCt = zeros(size(pp_grid));

for i = 1:length(pp)
    for j = 1:length(llc)
        for k = 1:length(mmu)
            dCt(j,i,k) = Fn(pp(i),h,0,1,mmu(k),llc(j))*k_model{2};
        end
    end
end

xslice = [0:pi/6:2*pi];
yslice = [];
zslice = [];
figure;
slc = slice(pp_grid, llc_grid, mmu_grid, dCt, xslice, yslice, zslice);
for i = 1:size(slc)
    set(slc(i),'EdgeAlpha',0)
end
xlabel('psi'); ylabel('lambda_c');zlabel('mu');
cb = colorbar;
cb.Label.String = 'dCt';

%%
uu = linspace(-10,10,30)';
vv = linspace(0,2,2)'; ww = linspace(-5,5,30)';

[uu_grid,vv_grid,ww_grid] = meshgrid(uu,vv,ww);

Cz = zeros(size(uu_grid));
Cl = zeros(size(uu_grid));
Cm = zeros(size(uu_grid));
for i = 1:length(uu)
    for j = 1:length(vv)
        for k = 1:length(ww)
            Cz(j,i,k) = P33(uu(i),vv(j),ww(k))*k_model{1};
            Cl(j,i,k) = P33(uu(i),vv(j),ww(k))*k_model{end-1};
            Cm(j,i,k) = P33(uu(i),vv(j),ww(k))*k_model{end};
        end
    end
end

xslice = [];
yslice = [0 2];
zslice = [];
figure;
slc = slice(uu_grid, vv_grid, ww_grid, Cm, xslice, yslice, zslice);
for i = 1:size(slc)
    set(slc(i),'EdgeAlpha',0)
end
xlabel('u'); ylabel('v');zlabel('w');
cb = colorbar;
cb.Label.String = 'Cm';

return;
%%
aalpha = linspace(-pi/2,pi/2,20); 
bbeta  = linspace(-pi,pi,30);

CT = zeros(length(aalpha),length(bbeta));
CM = zeros(length(aalpha),length(bbeta));
CL = zeros(length(aalpha),length(bbeta));
for i = 1:length(aalpha)
    for j = 1:length(bbeta)
        CT(i,j) = P32(aalpha(i),abs(bbeta(j)),1)*k_model{1};
        CL(i,j) = P32(aalpha(i),bbeta(j),1)*k_model{end-1};
        CM(i,j) = P32(aalpha(i),bbeta(j),1)*k_model{end};      
    end
end
figure
subplot(3,1,1)
surf(aalpha,bbeta,CT'); xlabel('alpha'); ylabel('beta'); zlabel('Cz');
subplot(3,1,2)
surf(aalpha,bbeta,CL'); xlabel('alpha'); ylabel('beta'); zlabel('Cl');
subplot(3,1,3)
surf(aalpha,bbeta,CM'); xlabel('alpha'); ylabel('beta'); zlabel('Cm');

return;
%% submodel check different characteristics
pp = linspace(0,2*pi,100)';
mmu = linspace(0,0.17,30)'; llc = linspace(-0.07,0.03,30)';

[pp_grid,llc_grid,mmu_grid] = meshgrid(pp,llc,mmu);

dCt12 = zeros(size(pp_grid));
dCt34 = zeros(size(pp_grid));

for i = 1:length(pp)
    for j = 1:length(llc)
        for k = 1:length(mmu)
            dCt12(j,i,k) = Fn(pp(i),h,0,1,mmu(k),llc(j))*K1;
            dCt34(j,i,k) = Fn(pp(i),h,0,1,mmu(k),llc(j))*K2;
        end
    end
end

xslice = [];
yslice = [];
zslice = [0:0.02:0.20];
figure;
subplot(1,2,1)
slc = slice(pp_grid, llc_grid, mmu_grid, dCt12, xslice, yslice, zslice);
for i = 1:size(slc)
    set(slc(i),'EdgeAlpha',0)
end
xlabel('psi'); ylabel('lambda_c');zlabel('mu');
cb = colorbar;
cb.Label.String = 'dCt12';
subplot(1,2,2)
slc = slice(pp_grid, llc_grid, mmu_grid, dCt34, xslice, yslice, zslice);
for i = 1:size(slc)
    set(slc(i),'EdgeAlpha',0)
end
xlabel('psi'); ylabel('lambda_c');zlabel('mu');
cb = colorbar;
cb.Label.String = 'dCt34';

%% scattered plot
% figure; du = 52840:60670;
figure
scatter3(psi_h1(du),lc1(du),mu1(du),2,Fn(psi_h1(du),h,0,1,mu1(du),lc1(du))*K1); hold on;
scatter3(psi_h3(du),lc3(du),mu3(du),2,Fn(psi_h3(du),h,0,1,mu3(du),lc3(du))*K2)
cb = colorbar;
cb.Label.String = 'dCt';
xlabel('psi'); ylabel('lambda_c');zlabel('mu'); title('omega13')

% figure;  du = 114200:125700;
figure;
scatter3(psi_h2(du),lc2(du),mu2(du),2,Fn(psi_h2(du),h,0,1,mu2(du),lc2(du))*K1); hold on;
scatter3(psi_h4(du),lc4(du),mu4(du),2,Fn(psi_h4(du),h,0,1,mu4(du),lc4(du))*K2)
cb = colorbar;
cb.Label.String = 'dCt';
xlabel('psi'); ylabel('lambda_c');zlabel('mu');title('omega24');
return;