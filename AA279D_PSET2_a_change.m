% ----------------------------------------------------

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

delta_a = 2;

a2 = a+delta_a;
e2 = e;
inc2 = inc+1;
RAAN2 = RAAN+1;
w2 = w+1;
anom2 = anom;

[r_ECI_init, v_ECI_init]=...
    OE2ECI(a, e, D2R*inc, D2R*RAAN, D2R*w, D2R*anom, mu);
[r2_ECI_init, v2_ECI_init]=...
    OE2ECI(a2, e2, D2R*inc2, D2R*RAAN2, D2R*w2, D2R*anom2, mu);

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
    min_finder(:,i) = norm(rel_pos_RTN(:,i));
end

figure
title('Osculating Orbital Elements with J2 Effects')
subplot(3,1,1);hold on;plot(rel_pos_RTN(1,:),'LineWidth',1.5);plot(ch_r_rtn(1,:),'r--','LineWidth',3);hold off;xlabel('Time (sec)');ylabel('\delta r_R (km)');
subplot(3,1,2);hold on;plot(rel_pos_RTN(2,:),'LineWidth',1.5);plot(ch_r_rtn(2,:),'r--','LineWidth',3);hold off;xlabel('Time (sec)');ylabel('\delta r_T (km)');
subplot(3,1,3);hold on;plot(rel_pos_RTN(3,:),'LineWidth',1.5);plot(ch_r_rtn(3,:),'r--','LineWidth',3);hold off;xlabel('Time (sec)');ylabel('\delta r_N (km)');

figure
title('Osculating Orbital Elements with J2 Effects')
subplot(3,1,1);hold on;plot(rel_vel_RTN(1,:),'LineWidth',1.5);plot(ch_v_rtn(1,:),'r--','LineWidth',3);hold off;xlabel('Time (sec)');ylabel('\delta v_R (km/s)');
subplot(3,1,2);hold on;plot(rel_vel_RTN(2,:),'LineWidth',1.5);plot(ch_v_rtn(2,:),'r--','LineWidth',3);hold off;xlabel('Time (sec)');ylabel('\delta v_T (km/s)');
subplot(3,1,3);hold on;plot(rel_vel_RTN(3,:),'LineWidth',1.5);plot(ch_v_rtn(3,:),'r--','LineWidth',3);hold off;xlabel('Time (sec)');ylabel('\delta v_N (km/s)');

% DELTA_V -------------------------------------------------

tout = 1:stepSize:stopTime;

x = rel_pos_RTN(1,:); dx = rel_vel_RTN(1,:);
y = rel_pos_RTN(2,:); dy = rel_vel_RTN(2,:);
z = rel_pos_RTN(3,:); dz = rel_vel_RTN(3,:);
r = rel_pos_RTN(4,:); dr = rel_vel_RTN(4,:);
th = rel_pos_RTN(5,:); dth = rel_vel_RTN(5,:);

Tc = 1/revs_per_day*24*60*60;
T1 = find(tout>0.99*Tc,1);

for t = 1:T1
    Vxm = dx(t)-dth(t)*y(t)+dr(t);
    Vym = dy(t)+dth(t)*(x(t)+r(t));
    Vzm = dz(t);
    r1 = sqrt((r(t)+x(t))^2+y(t)^2+z(t)^2);
    v1 = sqrt(Vxm^2+Vym^2+Vzm^2);
    term = -1+1/v1*sqrt(mu*(2*a-r1)/(a*r1));
    DVx = Vxm*term;
    DVy = Vym*term;
    DVz = Vzm*term;
    DVtot(t) = sqrt(DVx^2+DVy^2+DVz^2);
end

impulseTime = tout(DVtot == min(DVtot));
loc = find(DVtot == min(DVtot));

Vxm = dx(loc)-dth(loc)*y(loc)+dr(loc);
Vym = dy(loc)+dth(loc)*(x(loc)+r(loc));
Vzm = dz(loc);
r1 = sqrt((r(loc)+x(loc))^2+y(loc)^2+z(loc)^2);
v1 = sqrt(Vxm^2+Vym^2+Vzm^2);
term = -1+1/v1*sqrt(mu*(2*a-r1)/(a*r1));
DVx = Vxm*term;
DVy = Vym*term;
DVz = Vzm*term;

dV = [DVx;DVy;DVz]/stepSize

sim('relative_orbit_prop',[],options);

rel_pos_RTN = squeeze(rel_pos_RTN);
rel_vel_RTN = squeeze(rel_vel_RTN);

figure
plot(tout(1:T1)/Tc,DVtot)

figure
title('Osculating Orbital Elements with J2 Effects')
subplot(3,1,1);hold on;plot(tout,rel_pos_RTN(1,:),'LineWidth',1.5);plot(tout(loc),rel_pos_RTN(1,loc),'ro','LineWidth',3);hold off;xlabel('Time (sec)');ylabel('\delta r_R (km)');
subplot(3,1,2);hold on;plot(tout,rel_pos_RTN(2,:),'LineWidth',1.5);plot(tout(loc),rel_pos_RTN(2,loc),'ro','LineWidth',3);hold off;xlabel('Time (sec)');ylabel('\delta r_T (km)');
subplot(3,1,3);hold on;plot(tout,rel_pos_RTN(3,:),'LineWidth',1.5);plot(tout(loc),rel_pos_RTN(3,loc),'ro','LineWidth',3);hold off;xlabel('Time (sec)');ylabel('\delta r_N (km)');

figure
plot3(rel_pos_RTN(1,:),rel_pos_RTN(2,:),rel_pos_RTN(3,:),'LineWidth',2)



