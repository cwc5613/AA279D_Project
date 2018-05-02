%% PROBLEM SET THREE
% DESCRIPTIVE TEXT

clear;
clc;
close all;

TLE_Reader;

e = 0;

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
RAAN2 = RAAN+1;
w2 = w+1;
anom2 = anom;

%% Orbital Element Differences

del_a = a2-a;
del_e = e2-e;
del_inc = inc2-inc;
del_RAAN = RAAN2-RAAN;
del_w = w2-w;
del_anom2 = anom2-anom;

%% Position and Velocity in ECI

[r_ECI_init, v_ECI_init]=...
    OE2ECI(a, e, D2R*inc, D2R*RAAN, D2R*w, D2R*anom, mu);
[r2_ECI_init, v2_ECI_init]=...
    OE2ECI(a2, e2, D2R*inc2, D2R*RAAN2, D2R*w2, D2R*anom2, mu);

th0 = D2R*w + D2R*anom;
r0 = norm(r_ECI_init);
h = norm(cross(r_ECI_init,v_ECI_init));
thd0 = h/(r0^2);
ENERGY = norm(v_ECI_init)^2/2-mu/r0;

%% Relative Position and Velocity in RTN

T_ECI2RTN = ECI2RTN_rel(r_ECI_init,v_ECI_init);

rho_RTN_init = T_ECI2RTN*(r2_ECI_init-r_ECI_init);
rhod_RTN_init = T_ECI2RTN*(v2_ECI_init-v_ECI_init)-...
    cross([0;0;thd0],rho_RTN_init);

del_r = abs(norm(rho_RTN_init))*1000 <= Re

options = simset('SrcWorkspace','current');
orbitCount = 15;
stopTime = orbitCount/revs_per_day*24*60*60;
stepSize = 5;
%sim('HWC',[],options);
sim('HWC_accel',[],options);

pos_RTN_HCW = squeeze(pos_RTN_HCW);
vel_RTN_HCW = squeeze(vel_RTN_HCW);

pos_RTN_HCW = pos_RTN_HCW';
vel_RTN_HCW = vel_RTN_HCW';

figure
plot3(pos_RTN_HCW(:,1),pos_RTN_HCW(:,2),pos_RTN_HCW(:,3),...
    'LineWidth',2)
xlabel('r_R (km)')
ylabel('r_T (km)')
zlabel('r_N (km)')
title('Relative RTN Position from HCW Equations')

figure
plot3(vel_RTN_HCW(:,1),vel_RTN_HCW(:,2),vel_RTN_HCW(:,3),...
    'LineWidth',2)
xlabel('v_R (km/s)')
ylabel('v_T (km/s)')
zlabel('v_N (km/s)')
title('Relative RTN Velocity from HCW Equations')

figure
subplot(3,1,1)
plot(pos_RTN_HCW(:,2),pos_RTN_HCW(:,1),...
    'LineWidth',2)
xlabel('r_T (km)')
ylabel('r_R (km)')
title('Relative TR Position from HCW Equations')
subplot(3,1,2)
plot(pos_RTN_HCW(:,3),pos_RTN_HCW(:,1),...
    'LineWidth',2)
xlabel('r_N (km)')
ylabel('r_R (km)')
title('Relative NR Position from HCW Equations')
subplot(3,1,3)
plot(pos_RTN_HCW(:,2),pos_RTN_HCW(:,3),...
    'LineWidth',2)
xlabel('r_T (km)')
ylabel('r_N (km)')
title('Relative TN Position from HCW Equations')

figure
subplot(3,1,1)
plot(vel_RTN_HCW(:,2),vel_RTN_HCW(:,1),...
    'LineWidth',2)
xlabel('v_T (km/s)')
ylabel('v_R (km/s)')
title('Relative TR Velocity from HCW Equations')
subplot(3,1,2)
plot(vel_RTN_HCW(:,3),vel_RTN_HCW(:,1),...
    'LineWidth',2)
xlabel('v_N (km/s)')
ylabel('v_R (km/s)')
title('Relative NR Velocity from HCW Equations')
subplot(3,1,3);
plot(vel_RTN_HCW(:,2),vel_RTN_HCW(:,3),...
    'LineWidth',2)
xlabel('v_T (km/s)')
ylabel('v_N (km/s)')
title('Relative TN Velocity from HCW Equations')

%% Curvilinear Coordinates

[curve_r, curve_v] = rect2curv(rho_RTN_init, rhod_RTN_init, r0);

options = simset('SrcWorkspace','current');
orbitCount = 15;
stopTime = orbitCount/revs_per_day*24*60*60;
stepSize = 5;
sim('HWC_accel_curve',[],options);

curve_r_RTN_HCW = squeeze(curve_r_RTN_HCW);
curve_v_RTN_HCW = squeeze(curve_v_RTN_HCW);

[pos_RTN_HCW2, vel_RTN_HCW2] =...
    curv2rect(curve_r_RTN_HCW, curve_v_RTN_HCW, r0*ones(1,size(curve_r_RTN_HCW,2)));

pos_RTN_HCW2 = pos_RTN_HCW2';
vel_RTN_HCW2 = vel_RTN_HCW2';

figure
subplot(3,1,1)
hold on;
plot(pos_RTN_HCW(:,2),pos_RTN_HCW(:,1),...
    'LineWidth',2)
plot(pos_RTN_HCW2(:,2),pos_RTN_HCW2(:,1),...
    'LineWidth',2)
hold off;
xlabel('r_T (km)')
ylabel('r_R (km)')
title('Relative TR Position from HCW Equations')
subplot(3,1,2)
hold on;
plot(pos_RTN_HCW(:,3),pos_RTN_HCW(:,1),...
    'LineWidth',2)
plot(pos_RTN_HCW2(:,3),pos_RTN_HCW2(:,1),...
    'LineWidth',2)
hold off;
xlabel('r_N (km)')
ylabel('r_R (km)')
title('Relative NR Position from HCW Equations')
subplot(3,1,3)
hold on;
plot(pos_RTN_HCW(:,2),pos_RTN_HCW(:,3),...
    'LineWidth',2)
plot(pos_RTN_HCW2(:,2),pos_RTN_HCW2(:,3),...
    'LineWidth',2)
xlabel('r_T (km/s)')
ylabel('r_N (km/s)')
title('Relative TN Position from HCW Equations')

figure
subplot(3,1,1)
hold on;
plot(vel_RTN_HCW(:,2),vel_RTN_HCW(:,1),...
    'LineWidth',2)
plot(vel_RTN_HCW2(:,2),vel_RTN_HCW2(:,1),...
    'LineWidth',2)
hold off;
xlabel('v_T (km/s)')
ylabel('v_R (km/s)')
title('Relative TR Velocity from HCW Equations')
subplot(3,1,2)
hold on;
plot(vel_RTN_HCW(:,3),vel_RTN_HCW(:,1),...
    'LineWidth',2)
plot(vel_RTN_HCW2(:,3),vel_RTN_HCW2(:,1),...
    'LineWidth',2)
xlabel('v_N (km/s)')
ylabel('v_R (km/s)')
title('Relative NR Velocity from HCW Equations')
subplot(3,1,3);
plot(vel_RTN_HCW(:,3),vel_RTN_HCW(:,1),...
    'LineWidth',2)
plot(vel_RTN_HCW2(:,2),vel_RTN_HCW2(:,3),...
    'LineWidth',2)
xlabel('v_T (km/s)')
ylabel('v_N (km/s)')
title('Relative TN Velocity from HCW Equations')

figure
hold on
plot3(pos_RTN_HCW(:,1),pos_RTN_HCW(:,2),pos_RTN_HCW(:,3),...
    'LineWidth',2)
plot3(pos_RTN_HCW2(:,1),pos_RTN_HCW2(:,2),pos_RTN_HCW2(:,3),...
    'LineWidth',2)
xlabel('r_R (km)')
ylabel('r_T (km)')
zlabel('r_N (km)')
title('Relative RTN Position from HCW Equations')
hold off;

figure
hold on;
plot3(vel_RTN_HCW(:,1),vel_RTN_HCW(:,2),vel_RTN_HCW(:,3),...
    'LineWidth',2)
plot3(vel_RTN_HCW2(:,1),vel_RTN_HCW2(:,2),vel_RTN_HCW2(:,3),...
    'LineWidth',2)
xlabel('v_R (km/s)')
ylabel('v_T (km/s)')
zlabel('v_N (km/s)')
title('Relative RTN Velocity from HCW Equations')
hold off;

%% Integration Constants

% Rectilinear -----------------------------

c1 = 4*rho_RTN_init(1) + 2*rhod_RTN_init(2)/n
c2 = rhod_RTN_init(1)/n
c3 = -3*rho_RTN_init(1) - 2*rhod_RTN_init(2)/n
c4 = rho_RTN_init(2) - 2*rhod_RTN_init(1)/n
c5 = rhod_RTN_init(3)/n
c6 = rho_RTN_init(3)

% Curvilinear -----------------------------

cc1 = 4*curve_r(1) + 2*a*curve_v(2)/n
cc2 = curve_v(1)/n
cc3 = -3*curve_r(1) - 2*a*curve_v(2)/n
cc4 = a*curve_r(2) - 2*curve_v(1)/n
cc5 = a*curve_v(3)/n
cc6 = a*curve_r(3)
