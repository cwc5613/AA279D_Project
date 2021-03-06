%% PROBLEM SET TWO
% DESCRIPTIVE TEXT

clear;
clc;
close all;

TLE_Reader;

mu = 398600;
J2 = 1.081874*10^-3;
Re = 6378137*10^-3;
M0 = M;
tol = 10^-5;
D2R = pi/180;
R2D = 1/D2R;

a2 = a;
e2 = e;
inc2 = inc+1;
RAAN2 = RAAN;
w2 = w+1;
anom2 = anom+1;

[r_ECI_init, v_ECI_init]=...
    OE2ECI(a, e, inc, RAAN, w, anom, mu);
[r2_ECI_init, v2_ECI_init]=...
    OE2ECI(a2, e2, inc2, RAAN2, w2, anom2, mu);

th0 = D2R*w + D2R*anom;
r0 = norm(r_ECI_init);
h = norm(cross(r_ECI_init,v_ECI_init));
thd0 = h/(r0^2);
ENERGY = norm(v_ECI_init)^2/2-mu/r0;

T_ECI2RTN = ECI2RTN_rel(r_ECI_init,v_ECI_init);

rho_RTN_init = T_ECI2RTN*(r2_ECI_init-r_ECI_init);
rhod_RTN_init = T_ECI2RTN*(v2_ECI_init-v_ECI_init)-...
    cross([0;0;thd0],rho_RTN_init);

%rd0 = dot(v_ECI_init,r_ECI_init/norm(r_ECI_init));
rd0 = sqrt(2*(ENERGY+mu/r0)-thd0*h);

rel_pos0 = [rho_RTN_init;r0;th0];
rel_vel0 = [rhod_RTN_init;rd0;thd0];

options = simset('SrcWorkspace','current');
orbitCount = 10;
stopTime = orbitCount/revs_per_day*24*60*60;
stepSize = 5;
impulseTime = 0;
dV = [0;0;0];
sim('relative_orbit_prop',[],options);

rel_pos_RTN = squeeze(rel_pos_RTN);
rel_vel_RTN = squeeze(rel_vel_RTN);

figure
plot3(rel_pos_RTN(1,:),rel_pos_RTN(2,:),rel_pos_RTN(3,:));...
    xlabel('\delta r_R (km)');ylabel('\delta r_T (km)');zlabel('\delta r_N (km)');
figure
plot3(rel_vel_RTN(1,:),rel_vel_RTN(2,:),rel_vel_RTN(3,:));...
    xlabel('\delta v_R (km/s)');ylabel('\delta v_T (km/s)');zlabel('\delta v_N (km/s)');

% figure
% plot3(rel_vel_RTN(1,:),rel_vel_RTN(2,:),rel_vel_RTN(3,:),'LineWidth',2)

figure
subplot(3,1,1);plot(rel_pos_RTN(1,:));xlabel('Time (sec)');ylabel('\delta r_R (km)');...
        title('Relative RTN Position')
subplot(3,1,2);plot(rel_pos_RTN(2,:));xlabel('Time (sec)');ylabel('\delta r_T (km)');
subplot(3,1,3);plot(rel_pos_RTN(3,:));xlabel('Time (sec)');ylabel('\delta r_N (km)');

figure
subplot(3,1,1);plot(rel_vel_RTN(1,:));xlabel('Time (sec)');ylabel('\delta v_R (km/s)');...
        title('Relative RTN Velocity')
subplot(3,1,2);plot(rel_vel_RTN(2,:));xlabel('Time (sec)');ylabel('\delta v_T (km/s)');
subplot(3,1,3);plot(rel_vel_RTN(3,:));xlabel('Time (sec)');ylabel('\delta v_N (km/s)');

f0 = 0;
v0 = v_ECI_init;
r0 = r_ECI_init;

sim('NUM_PROP',[],options);

ch_r_ECI = reshape(NUM_r,[size(NUM_r,1),size(NUM_r,3)]);
ch_v_ECI = reshape(NUM_v,[size(NUM_v,1),size(NUM_v,3)]);

v0 = v2_ECI_init;
r0 = r2_ECI_init;

sim('NUM_PROP',[],options);

dep_r_ECI = reshape(NUM_r,[size(NUM_r,1),size(NUM_r,3)]);
dep_v_ECI = reshape(NUM_v,[size(NUM_v,1),size(NUM_v,3)]);

for i = 1:size(ch_r_ECI,2)
    %     ch_r = ch_r_ECI(:,i)/norm(ch_r_ECI(:,i));
    %     ch_h = cross(ch_r_ECI(:,i),ch_v_ECI(:,i))/...
    %         norm(cross(ch_r_ECI(:,i),ch_v_ECI(:,i)));
    ch_T_ECI2RTN = ECI2RTN_rel(ch_r_ECI(:,i),ch_v_ECI(:,i));
    ch_r_rtn(:,i) = ch_T_ECI2RTN*(dep_r_ECI(:,i)-ch_r_ECI(:,i));
    ch_v_rtn(:,i) = ch_T_ECI2RTN*(dep_v_ECI(:,i)-ch_v_ECI(:,i))-...
        cross([0;0;norm(cross(ch_r_ECI(:,i),ch_v_ECI(:,i)))/...
        (norm(ch_r_ECI(:,i))^2)], ch_r_rtn(:,i));
end

% figure
% plot3(ch_r_rtn(1,:),ch_r_rtn(2,:),ch_r_rtn(3,:),'r','LineWidth',2)
% 
% figure
% hold on
% plot3(ch_r_rtn(1,:),ch_r_rtn(2,:),ch_r_rtn(3,:),'r--','LineWidth',3)
% plot3(rel_pos_RTN(1,:),rel_pos_RTN(2,:),rel_pos_RTN(3,:),'LineWidth',1.5)
% hold off
% 
% figure
% hold on
% plot3(ch_v_rtn(1,:),ch_v_rtn(2,:),ch_v_rtn(3,:),'r--','LineWidth',3)
% plot3(rel_vel_RTN(1,:),rel_vel_RTN(2,:),rel_vel_RTN(3,:),'LineWidth',1.5)
% hold off

figure
subplot(3,1,1);...
    hold on;...
    plot(rel_pos_RTN(1,:),'LineWidth',1.5);...
    plot(ch_r_rtn(1,:),'r--','LineWidth',3);...
    hold off;...
    xlabel('Time (sec)');...
    ylabel('\delta r_R (km)');...
    title('Relative RTN Position Vectors');
subplot(3,1,2);...
    hold on;...
    plot(rel_pos_RTN(2,:),'LineWidth',1.5);...
    plot(ch_r_rtn(2,:),'r--','LineWidth',3);...
    hold off;...
    xlabel('Time (sec)');...
    ylabel('\delta r_T (km)');
subplot(3,1,3);...
    hold on;...
    plot(rel_pos_RTN(3,:),'LineWidth',1.5);...
    plot(ch_r_rtn(3,:),'r--','LineWidth',3);...
    hold off;...
    xlabel('Time (sec)');...
    ylabel('\delta r_N (km)');

figure
subplot(3,1,1);...
    hold on;...
    plot(rel_vel_RTN(1,:),'LineWidth',1.5);...
    plot(ch_v_rtn(1,:),'r--','LineWidth',3);...
    hold off;...
    xlabel('Time (sec)');...
    ylabel('\delta v_R (km)');...
    title('Relative RTN Velocity Vectors');
subplot(3,1,2);...
    hold on;...
    plot(rel_vel_RTN(2,:),'LineWidth',1.5);...
    plot(ch_v_rtn(2,:),'r--','LineWidth',3);...
    hold off;...
    xlabel('Time (sec)');...
    ylabel('\delta v_T (km)');
subplot(3,1,3);...
    hold on;...
    plot(rel_vel_RTN(3,:),'LineWidth',1.5);...
    plot(ch_v_rtn(3,:),'r--','LineWidth',3);...
    hold off;...
    xlabel('Time (sec)');...
    ylabel('\delta v_N (km)');

figure
subplot(3,1,1);...
    plot(rel_pos_RTN(1,:)-ch_r_rtn(1,:),'LineWidth',1.5);...
    xlabel('Time (sec)');...
    ylabel('\Delta\delta r_R (km)');...
    title('Relative RTN Position Error');
subplot(3,1,2);...
    plot(rel_pos_RTN(2,:)-ch_r_rtn(2,:),'LineWidth',1.5);...
    xlabel('Time (sec)');...
    ylabel('\Delta\delta r_T (km)');...
    subplot(3,1,3);...
    plot(rel_pos_RTN(3,:)-ch_r_rtn(3,:),'LineWidth',1.5);...
    xlabel('Time (sec)');...
    ylabel('\Delta\delta r_N (km)');

figure
subplot(3,1,1);...
    plot(rel_vel_RTN(1,:)-ch_v_rtn(1,:),'LineWidth',1.5);...
    xlabel('Time (sec)');...
    ylabel('\Delta\delta v_R (km)');...
    title('Relative RTN Velocity Error');
subplot(3,1,2);...
    plot(rel_vel_RTN(2,:)-ch_v_rtn(2,:),'LineWidth',1.5);...
    xlabel('Time (sec)');...
    ylabel('\Delta\delta v_T (km)');...
    subplot(3,1,3);...
    plot(rel_vel_RTN(3,:)-ch_v_rtn(3,:),'LineWidth',1.5);...
    xlabel('Time (sec)');...
    ylabel('\Delta\delta v_N (km)');