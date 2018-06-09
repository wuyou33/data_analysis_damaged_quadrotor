function [OT_a,OB_a,Delay] = align_signals(OT,OB,method)

switch method
    case 'phi'
        [s1,s2,Delay] = alignsignals(OB.phi_ot*57.3,OT.PHI);
    case 'theta'
        [s1,s2,Delay] = alignsignals(OB.theta*57.3,OT.THETA);          
    case 'Z'
        [s1,s2,Delay] = alignsignals(OB.PosNED(:,3),-OT.posCO_G(:,2)-(-OT.posCO_G(1,2)-OB.PosNED(1,3))); %using height
    case 'Z_top' %Use the upper 20% height to align signals
        hthr = min(OB.PosNED(:,3))*0.8;
        h1 = find(OB.PosNED(:,3)<=hthr); h2 = find(-OT.posCO_G(:,2)-(-OT.posCO_G(1,2)-OB.PosNED(1,3))<=hthr);
        hh1 = -0*ones(size(OB.PosNED(:,3))); hh2 = -0*ones(size(-OT.posCO_G(:,2)));
        hh1(h1) = OB.PosNED(h1,3); hh2(h2) = -OT.posCO_G(h2,2)-(-OT.posCO_G(1,2)-OB.PosNED(1,3));
        [s1,s2,Delay] = alignsignals(-hh1,-hh2);
    case 'VZ'
        [s1,s2,Delay] = alignsignals(OB.VelNED(:,3),OT.vel_E(:,3)); %using Vz
    case 'PQR'        
        pqr_filt_ot = [butterworth(OT.P,4,5/512);butterworth(OT.Q,4,5/512);butterworth(OT.R,4,5/512)];
        pqr_filt_ob = [butterworth(OB.p,4,5/512);butterworth(OB.q,4,5/512);butterworth(OB.r,4,5/512)];
        [s1,s2,Delay] = alignsignals(pqr_filt_ob,pqr_filt_ot);
    case 'P'        
        pqr_filt_ot = [butterworth(OT.P,4,5/512);];
        pqr_filt_ob = [butterworth(OB.p,4,5/512);];
        [s1,s2,Delay] = alignsignals(pqr_filt_ob,pqr_filt_ot);
    case 'Q'        
        pqr_filt_ot = [butterworth(OT.Q,4,5/512);];
        pqr_filt_ob = [butterworth(OB.q,4,5/512);];
        [s1,s2,Delay] = alignsignals(pqr_filt_ob,pqr_filt_ot);          
    case 'VZPQR'
        pqr_filt_ot = [butterworth(OT.P,4,5/512);butterworth(OT.Q,4,5/512);butterworth(OT.R,4,5/512)];
        pqr_filt_ob = [butterworth(OB.P,4,5/512);butterworth(OB.Q,4,5/512);butterworth(OB.R,4,5/512)];
        [s1,s2,Delay] = alignsignals([pqr_filt_ob;OB.VelNED(:,3)],...
                                     [pqr_filt_ot*57.3;OT.vel_E(:,3)]);
end
% %%
%%
L_OB = length(OB.TIME);
L_OT  = length(OT.TIME);

fields_OT = fieldnames(OT);
fields_OB = fieldnames(OB);

OT_a = struct; OB_a = struct;
if Delay >= 0
    L_SYN = min(L_OT-Delay,L_OB);
    for i = 1:length(fields_OT)
        OT_a.(fields_OT{i}) = OT.(fields_OT{i})(Delay+1:Delay+L_SYN,:);
    end
    for i = 1:length(fields_OB)
        OB_a.(fields_OB{i}) = OB.(fields_OB{i})(1:L_SYN,:);
    end
else
    L_SYN = min(L_OT,L_OB+Delay);
    for i = 1:length(fields_OT)
        OT_a.(fields_OT{i}) = OT.(fields_OT{i})(1:L_SYN,:);
    end
    for i = 1:length(fields_OB)
        OB_a.(fields_OB{i}) = OB.(fields_OB{i})(1-Delay:L_SYN-Delay,:);
    end
end
%%
% Check align result
figure('position',[0 0,700,200])
% plot(OB_a.PosNED(:,3)); hold on
% plot(-OT_a.posCO_G(:,2)-(-OT_a.posCO_G(1,2)-OB_a.PosNED(1,3))); ylabel(method);title('Check alignment');
plot(s1);hold on;
plot(s2);
end