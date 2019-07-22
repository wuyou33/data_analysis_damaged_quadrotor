%%
clear all;

addpath('D:\damage controller\data_analysis\data')
addpath(genpath('D:\Data\_code_import_files'));

%DRF r<0
DATA_raw{1} = load('9_9877-91680.mat');
DATA_raw{2} = load('17_8445-86200.mat');
% DATA_raw{2} = load('18_6865-88380.mat');
% %DRF r>0
% DATA_raw{4} = load('7_23650-96460.mat');
DATA_raw{3} = load('19_12350-88220.mat');
DATA_raw{4} = load('20_7174-87430.mat');
% %SRF
DATA_raw{5} = load('32_29110-84740.mat'); %r<0 lf
DATA_raw{6} = load('23_30700-100200.mat'); %r>0 lb 
% DATA_raw{7} = load('25_31500-92750.mat'); %r<0 rb
DATA_raw{7} = load('35_22600-53420.mat'); %r>0 rf

%Norminal
DATA_raw{8} = load('92_24380-117810.mat'); % psi = 0;
DATA_raw{9} = load('93_1-63210.mat'); % psi = -90;
% DATA_raw{1} = load('92_24380-90710.mat'); % psi = 0;
% DATA_raw{2} = load('93_1-63210.mat'); % psi = -90;
DATA_raw{10} = load('122_1-161100.mat'); % up and down
DATA_raw{11} = load('102_1-238300.mat'); % % psi = 0; heavy

% DATA_raw{1} = load('20_79690-87430.mat');
% DATA_raw{1} = load('18_81820-88380.mat');
% DATA_raw{1} = load('75_114700-119800.mat');
% DATA_raw{1} = load('18_62130-69970.mat');
Nc = 5;

DATA = DATA_raw;
names = fieldnames(DATA_raw{1});
for i = 1:length(DATA_raw)
    for j = 1:length(names)
        if ~strcmp(names{j},'DU') && ~strcmp(names{j},'rc')
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
alpha = [];  psi_mod = []; theta = []; phi = []; psi = [];
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
   psi = [psi ;DATA{i}.psi];   
   phi = [phi; DATA{i}.phi];
   theta = [theta; DATA{i}.theta];  
end
du = 1:length(fx);

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

if va>0.001
    u_bar = u./va;
    v_bar = v./va;
    w_bar = w./va;
else
    u_bar = 0; v_bar = 0; w_bar = 0;
end

% AT0 = [P32(alpha,abs(beta),va.^2*rho*S)];
% AFX0 = [P32(sign(u_bar).*abs(u_bar),w_bar,va.^2*S*1.225/2.*sign(u_bar))];
% AFY0 = [P32(sign(v_bar).*abs(v_bar),w_bar,va.^2*S*1.225/2.*sign(v_bar))];

AM0     = [P32(abs(u_bar),w_bar,va.^2*S*1.225/2.*sign(u_bar))];
AL0     = [P32(abs(v_bar),w_bar,va.^2*S*1.225/2.*sign(v_bar))];
AN0     = [P32(abs(v_bar),abs(u_bar),va.^2*S*1.225/2.*sign(v_bar).*sign(u_bar))];
% AFZ0    = [P32(abs(w_bar),abs(v_bar),va.^2*S*1.225/2.*sign(w_bar))];
AFZ0    = [P1n(abs(w_bar),2,0,va.^2*S*1.225/2.*sign(w_bar))];

AFX0 = [P1n(abs(u_bar),2,0,va.^2*S*1.225/2.*sign(u_bar))];
AFY0 = [P1n(abs(v_bar),2,0,va.^2*S*1.225/2.*sign(v_bar))];
% AFZ0  = [P1n(abs(w_bar),2,0,va.^2*S*1.225/2.*sign(w_bar))];
% AM0 = [P1n(abs(u_bar.*w_bar),2,0,b*va.^2*S*1.225/2.*sign(u_bar.*w_bar))];
% AL0 = [P1n(abs(v_bar.*w_bar),2,0,b*va.^2*S*1.225/2.*sign(v_bar.*w_bar))];
% AN0 = [P1n(abs(v_bar),2,0,va.^2*S*1.225/2.*sign(v_bar).*sign(u_bar))];

% AFX0 = [];
% AFY0 = [];
% AT0  = [];
% AM0 = [];
% AL0 = [];
% AN0 = [];


AFX6 = [ones(length(du),1).*(u1.*omega1+u2.*omega2+u3.*omega3+u4.*omega4)];
AFY8 = [ones(length(du),1).*(v1.*omega1+v2.*omega2+v3.*omega3+v4.*omega4)];

% AM6  = [ones(length(du),1).*(u1.*omega1+u2.*omega2+u3.*omega3+u4.*omega4), ones(length(du),1).*(u1.^2.*omega1+u2.^2.*omega2+u3.^2.*omega3+u4.^2.*omega4)]*1e-4;
% AL8  = [ones(length(du),1).*(v1.*omega1+v2.*omega2+v3.*omega3+v4.*omega4), ones(length(du),1).*(v1.^2.*omega1+v2.^2.*omega2+v3.^2.*omega3+v4.^2.*omega4)]*1e-4;

AN6  = [ones(length(du),1).*(b*SL(1)*u1.*omega1+SL(2)*b*u2.*omega2+SL(3)*b*u3.*omega3+SL(4)*b*u4.*omega4)];
AN8  = [ones(length(du),1).*(l*SM(1)*v1.*omega1+l*SM(2)*v2.*omega2+l*SM(3)*v3.*omega3+l*SM(4)*v4.*omega4)];
AN11  = [ones(length(du),1).*(b*SL(1)*SN(1)*v1.*omega1+SL(2)*b*SN(2)*v2.*omega2+SL(3)*b*SN(3)*v3.*omega3+SL(4)*b*SN(4)*v4.*omega4)];
AN12  = [ones(length(du),1).*(l*SM(1)*SN(1)*u1.*omega1+l*SM(2)*SN(2)*u2.*omega2+l*SM(3)*SN(3)*u3.*omega3+l*SM(4)*SN(4)*u4.*omega4)];
AFX_coup = [ones(length(du),1).*(SN(1)*v1.*omega1+SN(2)*v2.*omega2+SN(3)*v3.*omega3+SN(4)*v4.*omega4)];
AFY_coup = [ones(length(du),1).*(SN(1)*u1.*omega1+SN(2)*u2.*omega2+SN(3)*u3.*omega3+SN(4)*u4.*omega4)];

%%

% AA = [{0} {0};
%       {[u,u.^2]} {0};
%       {0} {[v,v.^2]};
%       {0} {0};
%       {0} {0};
%       {0} {0}];
% % % frame T; rotor dT dM, no heading effect.

% AA =    [ {AT0}    {0}      {0}     {0}             {0}                 {0}      {0}       ; 
%           {0}      {AFX0}   {0}     {[AFX6]}        {[AFX_coup]}        {0}      {0}       ;   
%           {0}      {0}      {AFY0}  {[AFY8]}        {-[AFY_coup]}       {0}      {0}       ;
%           {0}      {0}      {0}     {0}             {0}                 {-AFY8}  {AFY_coup};
%           {0}      {0}      {0}     {0}             {0}                 {AFX6}   {AFX_coup};
%           {0}      {0}      {0}     {[AN6]+[AN8]}   {AN11-AN12}         {0}      {0}       ]; 
% %       
AA =    [ {AFZ0}    {0}      {0}     {0}             {0}                 {0}      {0}        {0}     {0}     {0}; 
          {0}      {AFX0}   {0}     {[AFX6]}        {[AFX_coup]}        {0}      {0}        {0}     {0}     {0};   
          {0}      {0}      {AFY0}  {[AFY8]}        {-[AFY_coup]}       {0}      {0}        {0}     {0}     {0};
          {0}      {0}      {0}     {0}             {0}                 {-AFY8}  {AFY_coup} {AL0}   {0}     {0};
          {0}      {0}      {0}     {0}             {0}                 {AFX6}   {AFX_coup} {0}     {AM0}   {0};
          {0}      {0}      {0}     {[AN6]+[AN8]}   {AN11-AN12}         {0}      {0}        {0}     {0}     {AN0}]; 

% TODO : check correctness of AM0 !!!

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

Z_est = [Fz(du)-fz(du);
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
% 

%% import the simple model
load('D:\surfdrive\DATA\AeroModels\Ct_model_BB2_surf.mat');
load('D:\surfdrive\DATA\AeroModels\Cq_model_BB2_surf.mat');
Ct0 = Ct_model_BB2_surf(0,0);
Cq0 = Cq_model_BB2_surf(0,0);

R = 0.075;
Kt0 = Ct0*(pi*R^2)*R^2*1.225;
Kq0 = Cq0*(pi*R^2)*R^2*1.225;
T_compare = Kt0 * ( omega1.^2 + omega2.^2 + omega3.^2 + omega4.^2);
Fx_compare = -0.0965 * u - 0.0014*u.^2;
Fy_compare = -0.1060 * v - 0.0044*v.^2;
My_compare = Kt0 * l * ( omega1.^2 + omega2.^2 - omega3.^2 - omega4.^2);
Mx_compare = Kt0 * b * ( omega1.^2 - omega2.^2 - omega3.^2 + omega4.^2);
Mz_compare = Kq0 * ( omega1.^2 - omega2.^2 + omega3.^2 - omega4.^2);

%%

% 
% figure
% plot(Z_est); hold on; 
% plot(Y);

Y_fz = 0; Y_fx = 0; Y_fy = 0; Y_mx = 0; Y_my = 0; Y_mz = 0;
for i = 1:Nk
    Y_fz = Y_fz + AA{1,i}*k_model{i};
    Y_fx = Y_fx + AA{2,i}*k_model{i};
    Y_fy = Y_fy + AA{3,i}*k_model{i};
    Y_mx = Y_mx + AA{4,i}*k_model{i};
    Y_my = Y_my + AA{5,i}*k_model{i};
    Y_mz = Y_mz + AA{6,i}*k_model{i};
end


R2fx = find_R2(Y_fx,Fx(du));
RMSfx = find_RMS(Y_fx,Fx(du));

R2fy = find_R2(Y_fy,Fy(du));
RMSfy = find_RMS(Y_fy,Fy(du));

R2fz = find_R2(Y_fz+fz(du),Fz(du));
RMSfz = find_RMS(Y_fz+fz(du),Fz(du));

R2mx = find_R2(mx(du)+Y_mx,Mx(du));
RMSmx = find_RMS(mx(du)+Y_mx,Mx(du));

R2my = find_R2(my(du)+Y_my,My(du));
RMSmy = find_RMS(my(du)+Y_my,My(du));

R2mz = find_R2(mz(du)+Y_mz,Mz(du));
RMSmz = find_RMS(mz(du)+Y_mz,Mz(du));

fprintf('R2f:\t%f\t%f\t%f\nRMSf:\t%f\t%f\t%f\n\n',R2fx,R2fy,R2fz,RMSfx,RMSfy,RMSfz);
fprintf('R2m:\t%f\t%f\t%f\nRMSm:\t%f\t%f\t%f\n\n',R2mx,R2my,R2mz,RMSmx,RMSmy,RMSmz);
CC = corr(A_est);
[row,col] = find(abs(CC)>0.99);

figure
plot(row,col,'.');

figure
subplot(3,1,1)
plot(Fx(du)); hold on;
plot(Y_fx); 
% R2 = find_R2(Fx(du),Y_fx);
% RMS = find_RMS(Fx(du),Y_fx);
% legend(num2str(R2),num2str(RMS))
plot(Fx_compare(du),'linestyle','-.');
subplot(3,1,2)
plot(Fy(du)); hold on;
plot(Y_fy); 
plot(Fy_compare(du),'linestyle','-.');
subplot(3,1,3)
plot(Fz(du)); hold on;
plot((fz(du)+Y_fz)); 
% plot(-fz0(du),'linestyle','-.');
plot(-T_compare(du),'linestyle','-.');

figure
subplot(3,1,1)
plot(Mx(du)); hold on;
plot(mx(du)+Y_mx); plot(Mx_compare(du),'linestyle','-.');  
subplot(3,1,2)
plot(My(du)); hold on;
plot(my(du)+Y_my); plot(My_compare(du),'linestyle','-.'); 
subplot(3,1,3)
plot(Mz(du)); hold on;
plot(mz(du)+Y_mz); plot(Mz_compare(du),'linestyle','-.'); 

load('D:\surfdrive\DATA\AeroModels\Ct_model_BB2_v2.mat');
load('D:\surfdrive\DATA\AeroModels\Cq_model_BB2_v2.mat');

% save('D:\system identification\quadrotor identification\models\model_individual_damage_BB2.mat','Cq_model_BB2_surf','Ct_model_BB2_surf','k_model','h')
save('D:\damage controller\data_analysis\system_identification\model_individual_damage_BB2_simple_10th_Feb.mat','k_Cq0','k_Ct0','k_model')
return;


%% check P33 model
clear Cx Cy Cm Cl Cn Cz;
% 
% AFX0 = [P32(sign(u_bar).*abs(u_bar),w_bar,va.^2*S*1.225/2.*sign(u_bar))];
% AFY0 = [P32(sign(v_bar).*abs(v_bar),w_bar,va.^2*S*1.225/2.*sign(v_bar))];
% AT0 = [P33(abs(u_bar),abs(v_bar),w_bar,va.^2*S*1.225/2)];
% AM0 = [P32(sign(u_bar).*abs(u_bar),w_bar,va.^2*S*1.225/2.*sign(u_bar))];
% AL0 = [P32(sign(v_bar).*abs(v_bar),w_bar,va.^2*S*1.225/2.*sign(v_bar))];
% AN0 = [P32(sign(u_bar).*sign(v_bar).*abs(v_bar),w_bar,va.^2*S*1.225/2.*sign(v_bar).*sign(u_bar))];

u_bar_vector = linspace(-1,1,21)';
v_bar_vector = linspace(-1,1,21)';
w_bar_vector = linspace(-1,1,21)';

for i = 1: length(u_bar_vector)
    for j = 1 :length(w_bar_vector)
        for k = 1:length(v_bar_vector)
%             Cx(i,j) = P32(abs(u_bar_vector(i)),w_bar_vector(j),sign(u_bar_vector(i))) * k_model{2}; 
%             Cy(k,j) = P32(abs(v_bar_vector(k)),w_bar_vector(j),sign(v_bar_vector(k))) * k_model{3};
            Cx(i,j) = P1n(abs(u_bar_vector(i)),2,0,sign(u_bar_vector(i))) * k_model{2}; 
            Cy(k,j) = P1n(abs(v_bar_vector(k)),2,0,sign(v_bar_vector(k))) * k_model{3};
%             Cz(k,j) = P32(abs(w_bar_vector(j)),abs(v_bar_vector(k)),sign(w_bar_vector(j))) * k_model{1};
            Cz(k,j) = P1n(abs(w_bar_vector(j)),2,0,sign(w_bar_vector(j))) * k_model{1}; 
            Cl(k,j) = P32(abs(v_bar_vector(k)),w_bar_vector(j),sign(v_bar_vector(k))) * k_model{8}/b;
            Cm(i,j) = P32(abs(u_bar_vector(i)),w_bar_vector(j),sign(u_bar_vector(i))) * k_model{9}/b;             
%             Cn(k,i) = P1n(abs(v_bar_vector(k)),2,0,sign(v_bar_vector(k))*sign(u_bar_vector(i)))*  k_model{10};
            Cn(k,i) = P32(abs(v_bar_vector(k)),abs(u_bar_vector(i)),sign(v_bar_vector(k)).*sign(u_bar_vector(i)))*  k_model{10}/b;
            
            if w_bar_vector(j)^2+v_bar_vector(k)^2>1
               Cz(k,j) = nan;
               Cl(k,j) = nan;
            end
            if u_bar_vector(i)^2+v_bar_vector(k)^2>1
               Cn(k,i) = nan;
            end
            if u_bar_vector(i)^2+w_bar_vector(j)^2>1
               Cm(i,j) = nan;
            end
        end
    end
end

[u_bar_grid,w_bar_grid] = meshgrid(u_bar_vector,w_bar_vector);
[v_bar_grid,~] = meshgrid(v_bar_vector,w_bar_vector);
[v_bar_grid2,u_bar_grid2] = meshgrid(v_bar_vector,u_bar_vector);
figure
subplot(3,2,1)
surf(u_bar_grid,w_bar_grid,Cx'); xlabel('u'); ylabel('w'); zlabel('Cx');
hold on;
subplot(3,2,2)
surf(v_bar_grid,w_bar_grid,Cy'); xlabel('v'); ylabel('w'); zlabel('Cy');
hold on;
subplot(3,2,3)
surf(v_bar_grid,w_bar_grid,Cl'); xlabel('v'); ylabel('w'); zlabel('Cl');
hold on;
subplot(3,2,4)
surf(u_bar_grid,w_bar_grid,Cm'); xlabel('u'); ylabel('w'); zlabel('Cm');
hold on;
subplot(3,2,5)
surf(v_bar_grid2,u_bar_grid2,Cn'); xlabel('v'); ylabel('u'); zlabel('Cn');
hold on;
subplot(3,2,6)
surf(v_bar_grid,w_bar_grid,Cz'); xlabel('v'); ylabel('w'); zlabel('Cz');
hold on;

figure
fsub = subplot(2,3,1)
plot(u_bar_vector,Cx(:,1)); xlabel('$\bar{u}$'); ylabel('$C_x$');
hold on;
pos = get(fsub,'position');
set(fsub,'position',pos+[-.05,0,0,0]);
fsub = subplot(2,3,2)
plot(v_bar_vector,Cy(:,1)); xlabel('$\bar{v}$'); ylabel('$C_y$');
hold on;
pos = get(fsub,'position');
set(fsub,'position',pos+[-.05,0,0,0]);
fsub = subplot(2,3,3)
plot(w_bar_vector,Cz(11,:)); xlabel('$\bar{w}$'); ylabel('$C_z$'); hold on;
% plot(w_bar_vector,Cz(16,:));
% [leg3,] =legend('$\bar{v}$=0','$\bar{v}$=0.5');
% set(leg3,'Interpreter','latex');
pos = get(fsub,'position');
set(fsub,'position',pos+[-.05,0,0,0]);
fsub = subplot(2,3,4)
plot(v_bar_vector,Cl(:,16)); xlabel('$\bar{v}$'); ylabel('$C_l$'); hold on;
plot(v_bar_vector,Cl(:,11));
plot(v_bar_vector,Cl(:,5)); 
% plot(v_bar_vector,Cl(:,1)); 
leg4 = legend('$\bar{w}$=0.5','$\bar{w}$=0.0','$\bar{w}$=-0.5'); 
set(leg4,'Interpreter','latex');
pos = get(fsub,'position');
set(fsub,'position',pos+[-.05,0,0,0]);
fsub = subplot(2,3,5)
plot(u_bar_vector,Cm(:,16)); xlabel('$\bar{u}$'); ylabel('$C_m$'); hold on;
plot(u_bar_vector,Cm(:,11));
plot(u_bar_vector,Cm(:,5)); 
leg5 = legend('$\bar{w}$=0.5','$\bar{w}$=0.0','$\bar{w}$=-0.5');
set(leg5,'Interpreter','latex');
pos = get(fsub,'position');
set(fsub,'position',pos+[-.05,0,0,0]);
fsub = subplot(2,3,6)
% plot(v_bar_vector,Cn(:,16)); xlabel('v'); ylabel('Cn'); hold on;
% plot(v_bar_vector,Cn(:,11));
% plot(v_bar_vector,Cn(:,5)); 
% legend('u=0.5','u=0.0','\bar{u}=-0.5');
surf(v_bar_grid2,u_bar_grid2,Cn'); xlabel('$\bar{v}$'); ylabel('$\bar{u}$'); zlabel('$C_n$');
hold on; view([0,0,1]);
pos = get(fsub,'position');
set(fsub,'position',pos+[-.05,0,0,0]);
% leg = get(gca,'legend'); set(leg,'Interpreter','latex');
c=colorbar;
c.Label.String = 'C_n'