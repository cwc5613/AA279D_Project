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

delta_a = 4;

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
orbitCount = 4;
stopTime = orbitCount/revs_per_day*24*60*60;
stepSize = 1;
impulseTime = 0;
deltaV_theta = 0;
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
    apoapsis_finder(:,i) = norm(dep_r_ECI(:,i));
end

figure
hold on
plot3(ch_r_rtn(1,:),ch_r_rtn(2,:),ch_r_rtn(3,:),'r--','LineWidth',3)
plot3(rel_pos_RTN(1,:),rel_pos_RTN(2,:),rel_pos_RTN(3,:),'LineWidth',1.5)
hold off

dep_h = norm(cross(r2_ECI_init,v2_ECI_init));
dep_r = norm(r2_ECI_init);
dep_p = a2*(1-e2^2);

deltaV_theta = (delta_a*dep_r*dep_h)/(2*a2^2*dep_p)
[apoapsis,impulseTime] = min(apoapsis_finder)

sim('relative_orbit_prop',[],options);

rel_pos_RTN = squeeze(rel_pos_RTN);
rel_vel_RTN = squeeze(rel_vel_RTN);

ENERGY1 = norm(v_ECI_init)^2/2-mu/norm(r_ECI_init);
ENERGY2 = norm(rel_vel_RTN(1:3,end))^2/2-...
    mu/ norm(rel_pos_RTN(1:3,end));

figure
plot3(rel_pos_RTN(1,:),rel_pos_RTN(2,:),rel_pos_RTN(3,:),'LineWidth',2)



