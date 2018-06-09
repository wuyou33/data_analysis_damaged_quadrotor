function [OB_ekf,OT_ekf,ekf_interval] = data_for_ekf(OB,OT)


i_start = find(OT.Z<-0.2,1,'first');
i_end   = find(OT.Z<-0.2,1,'last');
ekf_interval = i_start:i_end;

OB_ekf.TIME = OB.TIME(ekf_interval);
OB_ekf.P = OB.P(ekf_interval);
OB_ekf.Q = OB.Q(ekf_interval);
OB_ekf.R = OB.R(ekf_interval);
OB_ekf.AX = OB.AX(ekf_interval);
OB_ekf.AY = OB.AY(ekf_interval);
OB_ekf.AZ = OB.AZ(ekf_interval);

OT_ekf.PHI = OT.PHI(ekf_interval);
OT_ekf.THETA = OT.THETA(ekf_interval);
OT_ekf.PSI = OT.PSI(ekf_interval);
OT_ekf.velCG_E = OT.velCG_E(ekf_interval,:);
end