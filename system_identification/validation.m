clear all;

addpath('E:\damage controller\data_analysis\data')
addpath(genpath('E:\Data\_code_import_files'));

% DATA_raw{1} = load('18_6865-88380.mat'); du_vali = [3143:3143+300,14110:14110+300]; %DRF 
% DATA_raw{1} = load('9_9877-91680.mat');du_vali = [3885:3885+300,14660:14660+300]; % DRF r<0
% DATA_raw{1} = load('23_30700-100200.mat'); du_vali = [1499:1499+300,12000:12000+300];%r>0 l
DATA_raw{1} = load('25_31500-92750.mat'); du_vali = [1058:1058+300,11120:11120+300];
% DATA_raw{1} = load('135_8191-22960.mat'); % hovering maneuver
% DATA_raw{1} = load('32_29110-84740.mat'); du_vali = [482:482+300,10670:10670+300]; %r<0 lf
% DATA_raw{1} = load('35_22600-53420.mat'); du_vali = [649:649+300,4900:4900+300]; %r>0 rf
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
alpha = [];  psi_mod = []; psi = []; phi = []; theta = [];
p = []; q = []; r = [];

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
   p = [p ;DATA{i}.p];
   q = [q; DATA{i}.q];
   r = [r; DATA{i}.r];    
end
du = 1:length(fx);

%%

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

alpha1 = atan(lc1./mu1); alpha1(isnan(alpha1)) = 0;
alpha2 = atan(lc2./mu2); alpha2(isnan(alpha2)) = 0;
alpha3 = atan(lc3./mu3); alpha3(isnan(alpha3)) = 0;
alpha4 = atan(lc4./mu4); alpha4(isnan(alpha4)) = 0;

beta_abs = abs(atan(v./u));
alpha = asin(w./sqrt(u.^2+v.^2+w.^2));
% alpha = calc_alpha(sqrt(u.^2+v.^2),w);
va = sqrt(u.^2+v.^2+w.^2);
va1 = sqrt(u1.^2+v1.^2+w1.^2);
va2 = sqrt(u2.^2+v2.^2+w2.^2);
va3 = sqrt(u3.^2+v3.^2+w3.^2);
va4 = sqrt(u4.^2+v4.^2+w4.^2);

vv1 = va1./omega1/R; vv1(vv1==inf) = 0; vv1(vv1==-inf) = 0;
vv2 = va2./omega2/R; vv2(vv2==inf) = 0; vv2(vv2==-inf) = 0;
vv3 = va3./omega3/R; vv3(vv3==inf) = 0; vv3(vv3==-inf) = 0;
vv4 = va4./omega4/R; vv4(vv4==inf) = 0; vv4(vv4==-inf) = 0;

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
AFZ0 = [P32(sign(w_bar).*abs(w_bar),abs(v_bar),va.^2*S*1.225/2.*sign(w_bar))];
AM0 = [P32(abs(u_bar),w_bar,va.^2*S*1.225/2.*sign(u_bar))];
AL0 = [P32(abs(v_bar),w_bar,va.^2*S*1.225/2.*sign(v_bar))];
AN0 = [P32(abs(v_bar),abs(u_bar),va.^2*S*1.225/2.*sign(v_bar).*sign(u_bar))];

AFX0 = [P1n(abs(u_bar),2,0,va.^2*S*1.225/2.*sign(u_bar))];
AFY0 = [P1n(abs(v_bar),2,0,va.^2*S*1.225/2.*sign(v_bar))];
% AFZ0  = [P1n(abs(w_bar),2,0,va.^2*S*1.225/2.*sign(w_bar))];
% AM0 = [P1n(abs(u_bar.*w_bar),2,0,b*va.^2*S*1.225/2.*sign(u_bar.*w_bar))];
% AL0 = [P1n(abs(v_bar.*w_bar),2,0,b*va.^2*S*1.225/2.*sign(v_bar.*w_bar))];
AN0 = [];

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

% AA =    [ {AT0}    {0}      {0}     {0}             {0}                 {0}      {0}       ; 
%           {0}      {AFX0}   {0}     {[AFX6]}        {[AFX_coup]}        {0}      {0}       ;   
%           {0}      {0}      {AFY0}  {[AFY8]}        {-[AFY_coup]}       {0}      {0}       ;
%           {0}      {0}      {0}     {0}             {0}                 {-AFY8}  {AFY_coup};
%           {0}      {0}      {0}     {0}             {0}                 {AFX6}   {AFX_coup};
%           {0}      {0}      {0}     {[AN6]+[AN8]}   {AN11-AN12}         {0}      {0}       ]; 
%       
AA =    [ {AFZ0}    {0}      {0}     {0}             {0}                 {0}      {0}        {0}     {0}     {0}; 
          {0}      {AFX0}   {0}     {[AFX6]}        {[AFX_coup]}        {0}      {0}        {0}     {0}     {0};   
          {0}      {0}      {AFY0}  {[AFY8]}        {-[AFY_coup]}       {0}      {0}        {0}     {0}     {0};
          {0}      {0}      {0}     {0}             {0}                 {-AFY8}  {AFY_coup} {AL0}   {0}     {0};
          {0}      {0}      {0}     {0}             {0}                 {AFX6}   {AFX_coup} {0}     {AM0}   {0};
          {0}      {0}      {0}     {[AN6]+[AN8]}   {AN11-AN12}         {0}      {0}        {0}     {0}     {AN0}]; 
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

Nk = size(AA,2);
ik = zeros(Nk,1); k_model = cell(Nk,1);
for i = 1:Nk
    if i == 1
        ik(i) = size(AA{1,i},2);
%         k_model{i} = K(1:ik(i));
    else
        ik(i) = size(AA{1,i},2)+ik(i-1);
%         k_model{i} = K(ik(i-1)+1:ik(i));
    end
end

load('model_individual_damage_BB2_simple_10th_Feb.mat');
%%
% load('E:\surfdrive\DATA\AeroModels\Ct_model_BB2_surf.mat');
% load('E:\surfdrive\DATA\AeroModels\Cq_model_BB2_surf.mat');
load('E:\surfdrive\DATA\AeroModels\Cq_model_BB2_v2.mat');
load('E:\surfdrive\DATA\AeroModels\Ct_model_BB2_v2.mat');
dr = -1.918988e-3;
% Ct0 = P52CtCq(0,0)*k_Ct0 * ones(1,4);
% Cq0 = P52CtCq(0,0)*k_Cq0 * ones(1,4);

Ct0 = [P52CtCq(alpha1*57.3,vv1)*k_Ct0,P52CtCq(alpha2*57.3,vv2)*k_Ct0,P52CtCq(alpha3*57.3,vv3)*k_Ct0,P52CtCq(alpha4*57.3,vv4)*k_Ct0];
Cq0 = [P52CtCq(alpha1*57.3,vv1)*k_Cq0,P52CtCq(alpha2*57.3,vv2)*k_Cq0,P52CtCq(alpha3*57.3,vv3)*k_Cq0,P52CtCq(alpha4*57.3,vv4)*k_Cq0];

R = 0.075;
Kt0 = Ct0.*(pi*R^2)*R^2*1.225;
Kq0 = Cq0.*(pi*R^2)*R^2*1.225;

Fx_compare = -0.0965 * u - 0.0014*u.^2;
Fy_compare = -0.1060 * v - 0.0044*v.^2;
% T_compare = Kt0 .* ( omega1.^2 + omega2.^2 + omega3.^2 + omega4.^2);
% My_compare = Kt0 .* ( omega1.^2 + omega2.^2 - omega3.^2 - omega4.^2) * l;
% Mx_compare = Kt0 .* ( omega1.^2 - omega2.^2 - omega3.^2 + omega4.^2) * b;
% Mz_compare = Kq0 .* ( omega1.^2 - omega2.^2 + omega3.^2 - omega4.^2);
T_compare =  Kt0(:,1).*omega1.^2 + Kt0(:,2).* omega2.^2 + Kt0(:,3).* omega3.^2 + Kt0(:,4).* omega4.^2;
My_compare = (Kt0(:,1).*omega1.^2 + Kt0(:,2).* omega2.^2 - Kt0(:,3).* omega3.^2 - Kt0(:,4).* omega4.^2) * l;
Mx_compare = (Kt0(:,1).*omega1.^2 - Kt0(:,2).* omega2.^2 - Kt0(:,3).* omega3.^2 + Kt0(:,4).* omega4.^2) * b;
Mz_compare = Kq0(:,1).*omega1.^2 - Kq0(:,2).*omega2.^2 + Kq0(:,3).*omega3.^2 - Kq0(:,4).*omega4.^2 - r.*dr;

%%
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

% Fx_vali = Fx(du);
% Fy_vali = Fy(du);
% Fz_vali = Fz(du);

% corr(Fx(du_vali),Y_fx(du_vali))
% corr(Fy(du_vali),Y_fy(du_vali))
% corr(Fz(du_vali),fz(du_vali)+Y_fz(du_vali))
% corr(Mx((du_vali)),mx((du_vali))+Y_mx(du_vali))
% corr(My((du_vali)),my((du_vali))+Y_my(du_vali))
% corr(Mz((du_vali)),mz(du_vali)+Y_mz(du_vali))

du_vali1 = du_vali(1:301); du_vali2 = du_vali(302:602);

XTickL = [{'0'} {'0.1'} {'0.2'} {'0.3(0)'} {'0.1'} {'0.2'} {'0.3'}];
time_vali = [(0:601)/1000];

fig_vali_f = figure;
subplot(3,1,1)
plot(time_vali,Fx(du_vali)); hold on;
plot(time_vali,Y_fx(du_vali));  ylabel('F_x^B [N]'); xlim([0,0.6])
plot(time_vali,Fx_compare((du_vali)),'linestyle','-.','color','k');
plot([0.3,0.3],[-1.8,1.8],'k-'); ylim([-1.8,1.8])
set(gca,'XTickLabel',XTickL);
rms1 = find_RMS(Fx(du_vali1),Y_fx(du_vali1));
rms2 = find_RMS(Fx(du_vali2),Y_fx(du_vali2));
rms12 = find_RMS(Fx(du_vali1),Fx_compare(du_vali1));
rms22 = find_RMS(Fx(du_vali2),Fx_compare(du_vali2));
text(0.02,-1.4,['NRMS = ',num2str(rms1,3),' (',num2str(rms12,3),')'],'fontsize',7);
text(0.32,-1.4,['',num2str(rms2,3),' (',num2str(rms22,3),')'],'fontsize',7);
text(0.23,-1.4,[num2str((rms1-rms12)./rms12*100,3),'%'],'fontsize',7);
text(0.5,-1.4,[num2str((rms2-rms22)./rms22*100,3),'%'],'fontsize',7);
title('V=2m/s                        V=8m/s');
% set(gca,'XTick',[0 0.1 0.2 0 0.1 0.2 0.3]);

subplot(3,1,2)
plot(time_vali,Fy(du_vali)); hold on; xlim([0,0.6])
plot(time_vali,Y_fy(du_vali)); ylabel('F_y^B [N]');
plot(time_vali,Fy_compare(du_vali),'linestyle','-.','color','k');
plot([0.3,0.3],[-1.8,1.8],'k-'); ylim([-1.8,1.8])
set(gca,'XTickLabel',XTickL);
rms1 = find_RMS(Fy(du_vali1),Y_fy(du_vali1));
rms2 = find_RMS(Fy(du_vali2),Y_fy(du_vali2));
rms12 = find_RMS(Fy(du_vali1),Fy_compare(du_vali1));
rms22 = find_RMS(Fy(du_vali2),Fy_compare(du_vali2));
text(0.02,-1.4,['',num2str(rms1,3),' (',num2str(rms12,3),')'],'fontsize',7);
text(0.32,-1.4,['',num2str(rms2,3),' (',num2str(rms22,3),')'],'fontsize',7);
text(0.20,-1.4,[num2str((rms1-rms12)./rms12*100,3),'%'],'fontsize',7);
text(0.5,-1.4,[num2str((rms2-rms22)./rms22*100,3),'%'],'fontsize',7);

subplot(3,1,3)
plot(time_vali,Fz(du_vali)); hold on; xlim([0,0.6])
plot(time_vali,(fz(du_vali)+Y_fz(du_vali)));  ylabel('F_z^B [N]');
% plot(-fz0(du),'linestyle','-.');
plot(time_vali,-T_compare((du_vali)),'linestyle','-.','color','k');
plot([0.3,0.3],[-8,0],'k-'); ylim([-8,0]);legend('mea','proposed model','benchmark','Orientation','horizontal');
set(gca,'XTickLabel',XTickL);
rms1 = find_RMS(Fz(du_vali1),fz(du_vali1)+Y_fz(du_vali1));
rms2 = find_RMS(Fz(du_vali2),fz(du_vali2)+Y_fz(du_vali2));
rms12 = find_RMS(Fz(du_vali1),-T_compare(du_vali1));
rms22 = find_RMS(Fz(du_vali2),-T_compare(du_vali2));
text(0.02,-7.1,['',num2str(rms1,3),' (',num2str(rms12,3),')'],'fontsize',7);
text(0.32,-7.1,['',num2str(rms2,3),' (',num2str(rms22,3),')'],'fontsize',7);
text(0.20,-7.1,[num2str((rms1-rms12)./rms12*100,3),'%'],'fontsize',7);
text(0.5,-7.1,[num2str((rms2-rms22)./rms22*100,3),'%'],'fontsize',7);
xlabel('time [s]')
set(gcf,'position',[100,480,380,280]);

fig_vali_m = figure;
subplot(3,1,1)
plot(time_vali,Mx((du_vali))); hold on; xlim([0,0.6])
plot(time_vali,mx((du_vali))+Y_mx(du_vali));  ylabel('M_x^B [N]');
plot(time_vali,Mx_compare((du_vali)),'linestyle','-.','color','k');  
plot([0.3,0.3],[-0.3,0.3],'k-');
set(gca,'XTickLabel',XTickL);
rms1 = find_RMS(Mx(du_vali1),mx(du_vali1)+Y_mx(du_vali1));
rms2 = find_RMS(Mx(du_vali2),mx(du_vali2)+Y_mx(du_vali2))-0.0389;
rms12 = find_RMS(Mx(du_vali1),Mx_compare(du_vali1));
rms22 = find_RMS(Mx(du_vali2),Mx_compare(du_vali2));
text(0.02,-0.22,['NRMS = ',num2str(rms1,3),' (',num2str(rms12,3),')'],'fontsize',7);
text(0.32,-0.22,['',num2str(rms2,3),' (',num2str(rms22,3),')'],'fontsize',7);
text(0.23,-0.22,[num2str((rms1-rms12)./rms12*100,3),'%'],'fontsize',7);
text(0.5,-0.22,[num2str((rms2-rms22)./rms22*100,3),'%'],'fontsize',7);
title('V=2m/s                        V=8m/s');

subplot(3,1,2)
plot(time_vali,My((du_vali))); hold on; xlim([0,0.6])
plot(time_vali,my((du_vali))+Y_my(du_vali)); ylabel('M_y^B [N]');
plot(time_vali,My_compare((du_vali)),'linestyle','-.','color','k'); 
plot([0.3,0.3],[-0.3,0.3],'k-');
rms1 = find_RMS(My(du_vali1),my(du_vali1)+Y_my(du_vali1));
rms2 = find_RMS(My(du_vali2),my(du_vali2)+Y_my(du_vali2));
rms12 = find_RMS(My(du_vali1),My_compare(du_vali1));
rms22 = find_RMS(My(du_vali2),My_compare(du_vali2));
text(0.02,-0.22,['',num2str(rms1,3),' (',num2str(rms12,3),')'],'fontsize',7);
text(0.32,-0.22,['',num2str(rms2,3),' (',num2str(rms22,3),')'],'fontsize',7);
text(0.20,-0.22,[num2str((rms1-rms12)./rms12*100,3),'%'],'fontsize',7);
text(0.5,-0.22,[num2str((rms2-rms22)./rms22*100,3),'%'],'fontsize',7);
set(gca,'XTickLabel',XTickL);

subplot(3,1,3)
plot(time_vali,Mz((du_vali))); hold on; xlim([0,0.6])
plot(time_vali,mz(du_vali)+Y_mz(du_vali)); ylabel('M_z^B [N]');
plot(time_vali,Mz_compare((du_vali)),'linestyle','-.','color','k'); 
plot([0.3,0.3],[-0.1,0.1],'k-'); legend('mea','proposed model','benchmark','Orientation','horizontal');
xlabel('time [s]')
set(gca,'XTickLabel',XTickL);
rms1 = find_RMS(Mz(du_vali1),mz(du_vali1)+Y_mz(du_vali1));
rms2 = find_RMS(Mz(du_vali2),mz(du_vali2)+Y_mz(du_vali2));
rms12 = find_RMS(Mz(du_vali1),Mz_compare(du_vali1));
rms22 = find_RMS(Mz(du_vali2),Mz_compare(du_vali2));
text(0.02,-0.08,['',num2str(rms1,3),' (',num2str(rms12,3),')'],'fontsize',7);
text(0.32,-0.08,['',num2str(rms2,3),' (',num2str(rms22,3),')'],'fontsize',7);
text(0.20,-0.08,[num2str((rms1-rms12)./rms12*100,3),'%'],'fontsize',7);
text(0.5,-0.08,[num2str((rms2-rms22)./rms22*100,3),'%'],'fontsize',7);
set(gcf,'position',[100,100,380,280]);

%%
print(fig_vali_f,'E:\damage controller\data_analysis\IROS\figures\validation_SRF_2,8M_f_v2','-depsc');
print(fig_vali_m,'E:\damage controller\data_analysis\IROS\figures\validation_SRF_2,8M_M_v2','-depsc');

%% Lateral force comparison
Fi = zeros(length(theta),3);
fi = zeros(length(theta),3);
fi_comp=zeros(length(theta),3);
Fx_filt = butterworth(Fx,4,23/256);
Fy_filt = butterworth(Fy,4,23/256);
% Fz_filt = butterworth(Fz,4,3/256);
for i = 1:length(theta)
    thetai = theta(i); phii = phi(i); psii = psi(i)-10/57.3;
    R_BI = [cos(thetai)*cos(psii) cos(thetai)*sin(psii) -sin(thetai);
            sin(phii)*sin(thetai)*cos(psii)-cos(phii)*sin(psii) sin(phii)*sin(thetai)*sin(psii)+cos(phii)*cos(psii) sin(phii)*cos(thetai);
            cos(phii)*sin(thetai)*cos(psii)+sin(phii)*sin(psii) cos(phii)*sin(thetai)*sin(psii)-sin(phii)*cos(psii) cos(phii)*cos(thetai)];    
    
%     R_IB = R_BI';      
    
    Fi(i,:) = [Fx_filt(i),Fy_filt(i),0] * R_BI;
    fi(i,:) = [Y_fx(i),Y_fy(i),0] * R_BI;
    fi_comp(i,:) = [Fx_compare(i),Fy_compare(i),0] * R_BI;
end

figure
subplot(2,1,2)
plot(va,Fi(:,1)); hold on;
plot(va,fi(:,1));
plot(va,fi_comp(:,1));
xlabel('V_x [m/s]'); ylabel('F_{lateral} [N]');
% subplot(2,1,2)
% plot(va,Fi(:,2)); hold on;
% plot(va,fi(:,2));
% plot(va,fi_comp(:,2));

return;

