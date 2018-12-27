%%
clear DATA_raw;

addpath('data')
addpath(genpath('E:\Data\_code_import_files'));

%DRF r<0
DATA_raw{1} = load('9_9877-91680.mat');
DATA_raw{2} = load('17_8445-86200.mat');
DATA_raw{3} = load('18_6865-88380.mat');
% %DRF r>0
DATA_raw{4} = load('7_23650-96460.mat');
DATA_raw{5} = load('19_12350-88220.mat');
DATA_raw{6} = load('20_7174-87430.mat');
% %SRF
% DATA_raw{1} = load('32_29110-84740.mat'); %r<0 lf
% DATA_raw{1} = load('23_30700-100200.mat'); %r>0 lb 
% DATA_raw{3} = load('25_31500-92750.mat'); %r<0 rb
% DATA_raw{4} = load('35_22600-53420.mat'); %r>0 rf
% 
% %Norminal
DATA_raw{7} = load('92_24380-90710.mat'); % psi = 0;
DATA_raw{8} = load('93_1-63210.mat'); % psi = -90;
DATA_raw{9} = load('122_1-161100.mat'); % up and down

% DATA_raw{7} = load('75_35590-119800.mat');
% DATA_raw{1} = load('20_79690-87430.mat');
% DATA_raw{1} = load('18_81820-88380.mat');
% DATA_raw{1} = load('75_114700-119800.mat');
% DATA_raw{1} = load('18_62130-69970.mat');
Nc = 5;

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
h2 = 5; %order of the Fourier series for fx fy

% psi_h1 = psi_mod(du) - 417; psi_h1 = psi_h1/57.3;
% psi_h2 = psi_mod(du) - 308; psi_h2 = psi_h2/57.3;
% psi_h3 = psi_mod(du) - 232; psi_h3 = psi_h3/57.3;
% psi_h4 = psi_mod(du) - 128; psi_h4 = psi_h4/57.3;
% 
% psi_h1 = mod(psi_h1,2*pi);
% psi_h2 = mod(psi_h2,2*pi);
% psi_h3 = mod(psi_h3,2*pi);
% psi_h4 = mod(psi_h4,2*pi);

beta= calc_beta2(u,v)*57.3;
psi_h1 = beta(du) - 413; psi_h1 = psi_h1/57.3;
psi_h2 = beta(du) - 307; psi_h2 = psi_h2/57.3;
psi_h3 = beta(du) - 233; psi_h3 = psi_h3/57.3;
psi_h4 = beta(du) - 127; psi_h4 = psi_h4/57.3;

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

alpha1 = atan(lc1./mu1);
alpha2 = atan(lc2./mu2);
alpha3 = atan(lc3./mu3);
alpha4 = atan(lc4./mu4);

beta_abs = abs(atan(v./u));
alpha = asin(w./sqrt(u.^2+v.^2+w.^2));
% alpha = calc_alpha(sqrt(u.^2+v.^2),w);
va = sqrt(u.^2+v.^2+w.^2);

% AT0 = [P32(alpha,abs(beta),va.^2*rho*S)];
AFX0 = [u.^2*rho*S];
AFY0 = [v.^2*rho*S];
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

AFX6 = [ones(length(du),1).*(u1.*omega1+u2.*omega2+u3.*omega3+u4.*omega4)];
AFX7 = [Fn(psi_h1,h2,0,u1.*omega1) + Fn(2*pi-psi_h2,h2,0,u2.*omega2) + Fn(psi_h3,h2,0,u3.*omega3) + Fn(2*pi-psi_h4,h2,0,u4.*omega4)];
AFX7fr = [Fn(psi_h1,h2,0,u1.*omega1) + Fn(2*pi-psi_h2,h2,0,u2.*omega2),...
        Fn(psi_h3,h2,0,u3.*omega3) + Fn(2*pi-psi_h4,h2,0,u4.*omega4)];
AFY8 = [ones(length(du),1).*(v1.*omega1+v2.*omega2+v3.*omega3+v4.*omega4)];
AFY9 = [Fn(psi_h1,h2,0,v1.*omega1) + Fn(2*pi-psi_h2,h2,0,v2.*omega2) + Fn(psi_h3,h2,0,v3.*omega3) + Fn(2*pi-psi_h4,h2,0,v4.*omega4)];
AFY9fr = [Fn(psi_h1,h2,0,v1.*omega1) + Fn(2*pi-psi_h2,h2,0,v2.*omega2),...
        Fn(psi_h3,h2,0,v3.*omega3) + Fn(2*pi-psi_h4,h2,0,v4.*omega4)];
AN6  = [ones(length(du),1).*(b*SL(1)*u1.*omega1+SL(2)*b*u2.*omega2+SL(3)*b*u3.*omega3+SL(4)*b*u4.*omega4)];
AN7  = [b*SL(1)*Fn(psi_h1,h2,0,u1.*omega1) + b*SL(2)*Fn(2*pi-psi_h2,h2,0,u2.*omega2) + b*SL(3)*Fn(psi_h3,h2,0,u3.*omega3) + b*SL(4)*Fn(2*pi-psi_h4,h2,0,u4.*omega4)];
AN7fr  = [b*SL(1)*Fn(psi_h1,h2,0,u1.*omega1) + b*SL(2)*Fn(2*pi-psi_h2,h2,0,u2.*omega2),...
          b*SL(3)*Fn(psi_h3,h2,0,u3.*omega3) + b*SL(4)*Fn(2*pi-psi_h4,h2,0,u4.*omega4)];
AN8  = [ones(length(du),1).*(l*SM(1)*v1.*omega1+l*SM(2)*v2.*omega2+l*SM(3)*v3.*omega3+l*SM(4)*v4.*omega4)];
AN9  = [l*SM(1)*Fn(psi_h1,h2,0,v1.*omega1) + l*SM(2)*Fn(2*pi-psi_h2,h2,0,v2.*omega2)+l*SM(3)*Fn(psi_h3,h2,0,v3.*omega3) + l*SM(4)*Fn(2*pi-psi_h4,h2,0,v4.*omega4)];
AN9fr  = [l*SM(1)*Fn(psi_h1,h2,0,v1.*omega1) + l*SM(2)*Fn(2*pi-psi_h2,h2,0,v2.*omega2),...
          l*SM(3)*Fn(psi_h3,h2,0,v3.*omega3) + l*SM(4)*Fn(2*pi-psi_h4,h2,0,v4.*omega4)];
AN10  = [SN(1)*P33(u1,v1,w1,omega1) + SN(2)*P33(u2,v2,w2,omega2) + SN(3)*P33(u3,v3,w3,omega3) + SN(4)*P33(u4,v4,w4,omega4)];
AN11  = [ones(length(du),1).*(b*SL(1)*SN(1)*v1.*omega1+SL(2)*b*SN(2)*v2.*omega2+SL(3)*b*SN(3)*v3.*omega3+SL(4)*b*SN(4)*v4.*omega4)];
AN12  = [ones(length(du),1).*(l*SM(1)*SN(1)*u1.*omega1+l*SM(2)*SN(2)*u2.*omega2+l*SM(3)*SN(3)*u3.*omega3+l*SM(4)*SN(4)*u4.*omega4)];
AN13 = [SN(1)*Fn(psi_h1,h,0,dynhead1*Area,mu1,lc1) + SN(2)*Fn(2*pi-psi_h2,h,0,dynhead2*Area,mu2,lc2) ...
         + SN(3)*Fn(psi_h3,h,0,dynhead3*Area,mu3,lc3) + SN(4)*Fn(2*pi-psi_h4,h,0,dynhead4*Area,mu4,lc4)];
AFX_coup = [ones(length(du),1).*(SN(1)*v1.*omega1+SN(2)*v2.*omega2+SN(3)*v3.*omega3+SN(4)*v4.*omega4)];
AFY_coup = [ones(length(du),1).*(SN(1)*u1.*omega1+SN(2)*u2.*omega2+SN(3)*u3.*omega3+SN(4)*u4.*omega4)];
% AL4 = [P32(alpha,sign(v).*beta_abs,va.^2*rho*S)];
% AM5 = [P32(sign(u).*alpha,beta_abs,va.^2*rho*S)];
AL4 = P33(u,abs(v),w,sign(v));
AM5 = P33(u,abs(v),w);

%%
% % % frame T + M ; rotor dT dM
% AA =    [{AT0} {AT1} {0}   {0}   {0}   {0}     {0}            {0}           {0}    {0}           {0}; 
%          {0}   {0}   {0}   {0}   {0}   {0}     {[AFX6 AFX7]}  {0}           {0}    {0}           {0};   
%          {0}   {0}   {0}   {0}   {0}   {0}     {0}            {[AFY8 AFY9]} {0}    {0}           {0};
%          {0}   {AL1} {AL2} {AL3} {AL4} {0}     {0}            {0}           {0}    {0}           {AFY8};
%          {0}   {AM1} {AM2} {AM3} {0}   {AM5}   {0}            {0}           {0}    {AFX6}        {0};
%          {0}   {0}   {0}   {0}   {0}   {0}     {[AN6 AN7]}    {[AN8 AN9]}   {AN10} {0}           {0} ]; 
% % % frame T + M ; rotor dT dM
AA =    [{AT0}      {AT1}	{0}        {0}          {0}             {0}             {0}      {0}        {0}     {0}         {0}         {0}; 
         {0}        {0}     {[AFX6]}   {[AFX_coup]} {0}             {0}             {0}      {0}        {0}     {0}         {0}         {0};   
         {0}        {0}     {0}        {0}          {[AFY8]}        {[AFY_coup]}    {0}      {0}        {0}     {0}         {0}         {0};
         {0}        {AL1}	{0}        {0}          {0}             {0}             {AFY8}   {AFY_coup} {0}     {0}         {0}         {0};
         {0}        {AM1}	{0}        {0}          {0}             {0}             {0}      {0}        {AFX6}  {AFX_coup}  {0}         {0};
         {0}        {0}     {[AN6]}    {AN11}       {[AN8]}         {AN12}          {0}      {0}        {0}     {0}         {AN13}      {0}]; 

%%
columneliminate = [];
for j = 1:size(AA,2)
    for i = 1:size(AA,1)
        if AA{i,j} ~=0 
           size_cell = size(AA{i,j}); 
           break;
        end
        if i == size(AA,1)
            columneliminate = [columneliminate j];
        end
    end
    for i = 1:size(AA,1)
        if AA{i,j} ==0 
            AA{i,j} = zeros(size_cell);
        end
    end
end
AA(:,columneliminate) = [];

A_est = cell2mat(AA);

Z_est = [-Fz(du)+fz(du);
         Fx(du);
         Fy(du)
         Mx(du)-mx(du);
         My(du)-my(du);
         Mz(du)-mz(du)];

sigma = ones(size(Z_est));
sigma(1:3*length(du)) = 5e-3;
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

R2f = find_R2(Y(1:3*length(du)),Z_est(1:3*length(du)));
RMSf = find_RMS(Y(1:3*length(du)),Z_est(1:3*length(du)));

R2m = find_R2(Y(3*length(du)+1:end),Z_est(3*length(du)+1:end));
RMSm = find_RMS(Y(3*length(du)+1:end),Z_est(3*length(du)+1:end));

fprintf('R2f:\t%f\nRMSf:\t%f\n\n',R2f,RMSf);
fprintf('R2m:\t%f\nRMSm:\t%f\n\n',R2m,RMSm);
CC = corr(A_est);
[row,col] = find(abs(CC)>0.95);
% 
figure
plot(row,col,'.');

figure
plot(Z_est); hold on; 
plot(Y);

Y_fz = 0; Y_fx = 0; Y_fy = 0; Y_mx = 0; Y_my = 0; Y_mz = 0;
for i = 1:Nk
    Y_fz = Y_fz + AA{1,i}*k_model{i};
    Y_fx = Y_fx + AA{2,i}*k_model{i};
    Y_fy = Y_fy + AA{3,i}*k_model{i};
    Y_mx = Y_mx + AA{4,i}*k_model{i};
    Y_my = Y_my + AA{5,i}*k_model{i};
    Y_mz = Y_mz + AA{6,i}*k_model{i};
end
figure

subplot(3,1,1)
plot(Fx(du)); hold on;
plot(Y_fx); 
% R2 = find_R2(Fx(du),Y_fx);
% RMS = find_RMS(Fx(du),Y_fx);
% legend(num2str(R2),num2str(RMS))
plot(fx(du),'linestyle','-.');
subplot(3,1,2)
plot(Fy(du)); hold on;
plot(Y_fy); 
plot(fy(du),'linestyle','-.');
subplot(3,1,3)
plot(-Fz(du)); hold on;
plot(Y_fz-fz(du)); 
plot(-fz(du),'linestyle','-.');

figure
subplot(3,1,1)
plot(Mx(du)); hold on;
plot(mx(du)+Y_mx); plot(mx(du),'linestyle','-.'); 
subplot(3,1,2)
plot(My(du)); hold on;
plot(my(du)+Y_my); plot(my(du),'linestyle','-.'); 
subplot(3,1,3)
plot(Mz(du)); hold on;
plot(mz(du)+Y_mz); plot(mz(du),'linestyle','-.'); 

load('E:\system identification\thrust_model\Ct_model_BB2.mat');
load('E:\system identification\thrust_model\Cq_model_BB2.mat');

% save('E:\system identification\quadrotor identification\models\model_individual_damage_BB2.mat','Cq_model_BB2_surf','Ct_model_BB2_surf','k_model','h')
save('E:\damage controller\bebop_v3.0_validation\data\model_individual_damage_BB2.mat','Cq_model_BB2_surf','Ct_model_BB2_surf','k_model','h')
% return;
%% submodel check same characteristics

pp = [0 127/57.3 pi/2 3*pi/4]';
mmu = linspace(0,0.25,30)'; llc = linspace(-0.1,0.1,30)';

[pp_grid,llc_grid,mmu_grid] = meshgrid(pp,llc,mmu);

dCt = zeros(size(pp_grid));
dCq = zeros(size(pp_grid));
Ct = zeros(size(pp_grid));
Cq = zeros(size(pp_grid));
for i = 1:length(pp)
    for j = 1:length(llc)
        for k = 1:length(mmu)
            aalpha = atan(llc(j)./mmu(k))*57.3;
            vva = sqrt(llc(j).^2+mmu(k).^2);
            dCt(j,i,k) = Fn(pp(i),h,0,1,mmu(k),llc(j))*k_model{2};
            if dCt(j,i,k)>= 0.007
                dCt(j,i,k) = 0.007;
            elseif dCt(j,i,k)< -0.007
                dCt(j,i,k) =  -0.007;
            end
            dCq(j,i,k) = Fn(pp(i),h,0,1,mmu(k),llc(j))*k_model{11};
            if dCq(j,i,k)>= 0.0004
                dCq(j,i,k) = 0.0004;
            elseif dCq(j,i,k)< -0.0007
                dCq(j,i,k) =  -0.0007;
            end
            Ct(j,i,k) = Ct_model_BB2_surf(aalpha,vva);
            Cq(j,i,k) = Cq_model_BB2_surf(aalpha,vva);            
        end
    end
end
%%
figure
subplot(2,2,1)
a = pcolor(mmu,llc,squeeze(dCt(:,1,:))); hold on; colorbar; title('\beta_i = 0')
xlabel('\mu'); ylabel('\lambda_c')
subplot(2,2,2)
a = pcolor(mmu,llc,squeeze(dCt(:,2,:))); hold on; colorbar; title('\beta_i = 90')
xlabel('\mu'); ylabel('\lambda_c')
subplot(2,2,3)
a = pcolor(mmu,llc,squeeze(dCt(:,3,:))); hold on; colorbar; title('\beta_i = 180')
xlabel('\mu'); ylabel('\lambda_c')
subplot(2,2,4)
a = pcolor(mmu,llc,squeeze(dCt(:,4,:))); hold on; colorbar; title('\beta_i = 270')
xlabel('\mu'); ylabel('\lambda_c')

figure
subplot(2,2,1)
a = pcolor(mmu,llc,squeeze(dCq(:,1,:))); hold on; colorbar; title('\beta_i = 0')
xlabel('\mu'); ylabel('\lambda_c')
subplot(2,2,2)
a = pcolor(mmu,llc,squeeze(dCq(:,2,:))); hold on; colorbar; title('\beta_i = 90')
xlabel('\mu'); ylabel('\lambda_c')
subplot(2,2,3)
a = pcolor(mmu,llc,squeeze(dCq(:,3,:))); hold on; colorbar; title('\beta_i = 180')
xlabel('\mu'); ylabel('\lambda_c')
subplot(2,2,4)
a = pcolor(mmu,llc,squeeze(dCq(:,4,:))); hold on; colorbar; title('\beta_i = 270')
xlabel('\mu'); ylabel('\lambda_c')
%%
figure
subplot(2,2,1)
a = surf(mmu,llc,squeeze(Ct(:,1,:))); hold on; set(a,'FaceAlpha',0.5,'LineStyle','none','FaceColor','r');
b = surf(mmu,llc,squeeze(dCt(:,1,:)+Ct(:,1,:))); title('\beta_i = 0');set(b,'FaceAlpha',0.5);
xlabel('\mu'); ylabel('\lambda_c');zlabel('C_t')
% plot(mu1,lc1,':'); plot(mu2,lc2,':'); plot(mu3,lc3,':'); plot(mu4,lc4,':');
subplot(2,2,2)
a = surf(mmu,llc,squeeze(Ct(:,2,:))); hold on;  set(a,'FaceAlpha',0.5,'LineStyle','none','FaceColor','r');
b = surf(mmu,llc,squeeze(dCt(:,2,:)+Ct(:,1,:))); title('\beta_i = 90');set(b,'FaceAlpha',0.5);
xlabel('\mu'); ylabel('\lambda_c');zlabel('C_t')
% plot(mu1,lc1,':'); plot(mu2,lc2,':'); plot(mu3,lc3,':'); plot(mu4,lc4,':');
subplot(2,2,3)
a = surf(mmu,llc,squeeze(Ct(:,3,:))); hold on; set(a,'FaceAlpha',0.5,'LineStyle','none','FaceColor','r');
b = surf(mmu,llc,squeeze(dCt(:,3,:)+Ct(:,1,:))); title('\beta_i = 180');set(b,'FaceAlpha',0.5);
xlabel('\mu'); ylabel('\lambda_c');zlabel('C_t')
% plot(mu1,lc1,':'); plot(mu2,lc2,':'); plot(mu3,lc3,':'); plot(mu4,lc4,':');
subplot(2,2,4)
a = surf(mmu,llc,squeeze(Ct(:,4,:))); hold on; set(a,'FaceAlpha',0.5,'LineStyle','none','FaceColor','r');
b = surf(mmu,llc,squeeze(dCt(:,4,:)+Ct(:,1,:))); title('\beta_i = 270');set(b,'FaceAlpha',0.5);
xlabel('\mu'); ylabel('\lambda_c');zlabel('C_t')
% plot(mu1,lc1,':'); plot(mu2,lc2,':'); plot(mu3,lc3,':'); plot(mu4,lc4,':');

figure
subplot(2,2,1)
a = surf(mmu,llc,squeeze(Cq(:,1,:))); hold on; set(a,'FaceAlpha',0.5,'LineStyle','none','FaceColor','r');
b = surf(mmu,llc,squeeze(dCq(:,1,:)+Cq(:,1,:))); title('\beta_i = 0');set(b,'FaceAlpha',0.5);
xlabel('\mu'); ylabel('\lambda_c');zlabel('C_q')
% plot(mu1,lc1,':'); plot(mu2,lc2,':'); plot(mu3,lc3,':'); plot(mu4,lc4,':');
subplot(2,2,2)
a = surf(mmu,llc,squeeze(Cq(:,2,:))); hold on; set(a,'FaceAlpha',0.5,'LineStyle','none','FaceColor','r');
b = surf(mmu,llc,squeeze(dCq(:,2,:)+Cq(:,1,:))); title('\beta_i = 90');set(b,'FaceAlpha',0.5);
xlabel('\mu'); ylabel('\lambda_c');zlabel('C_q')
% plot(mu1,lc1,':'); plot(mu2,lc2,':'); plot(mu3,lc3,':'); plot(mu4,lc4,':');
subplot(2,2,3)
a = surf(mmu,llc,squeeze(Cq(:,3,:))); hold on; set(a,'FaceAlpha',0.5,'LineStyle','none','FaceColor','r');
b = surf(mmu,llc,squeeze(dCq(:,3,:)+Cq(:,1,:))); title('\beta_i = 180');set(b,'FaceAlpha',0.5);
xlabel('\mu'); ylabel('\lambda_c');zlabel('C_q')
% plot(mu1,lc1,':'); plot(mu2,lc2,':'); plot(mu3,lc3,':'); plot(mu4,lc4,':');
subplot(2,2,4)
a = surf(mmu,llc,squeeze(Cq(:,4,:))); hold on; set(a,'FaceAlpha',0.5,'LineStyle','none','FaceColor','r');
b = surf(mmu,llc,squeeze(dCq(:,4,:)+Cq(:,1,:))); title('\beta_i = 270');set(b,'FaceAlpha',0.5);
xlabel('\mu'); ylabel('\lambda_c');zlabel('C_q')
% plot(mu1,lc1,':'); plot(mu2,lc2,':'); plot(mu3,lc3,':'); plot(mu4,lc4,':');

return;
%% slice curves
%Ct
xslice = [0 pi/4 pi/2 3*pi/4];
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

% Cq
xslice = [];
yslice = [0];
zslice = [0];
figure;
slc = slice(pp_grid, llc_grid, mmu_grid, dCq, xslice, yslice, zslice);
for i = 1:size(slc)
    set(slc(i),'EdgeAlpha',0)
end
xlabel('psi'); ylabel('lambda_c');zlabel('mu');
cb = colorbar;
cb.Label.String = 'dCq';



%% plot revised Ct and Cq in each heading angle

load('E:\system identification\thrust_model\Ct_model_BB2.mat');
load('E:\system identification\thrust_model\Cq_model_BB2.mat');


for i = 1:length(pp)
    for j = 1:length(llc)
        for k = 1:length(mmu)
            aalpha = atan(llc(j)./mmu(k));
            vva = sqrt(llc(j).^2+mmu(k).^2);
            Ct(j,i,k) = Ct_model_BB2_surf(aalpha,vva);
            Cq(j,i,k) = Cq_model_BB2_surf(aalpha,vva);
        end
    end
end

figure
surf(mmu,llc,squeeze(Ct(:,1,:)));
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