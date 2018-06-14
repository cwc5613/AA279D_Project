%% PROBLEM SET THREE
% DESCRIPTIVE TEXT

clear;
clc;
close all;

TLE_Reader;

mu = 398600;
J2 = 1.081874*10^-3;
Re = 6378137*10^-3;
M0 = M*pi/180;
tol = 10^-5;
D2R = pi/180;
R2D = 1/D2R;

anom = 0;
e = 0.1;

a2 = a;
e2 = e;
inc2 = inc-0.001;
RAAN2 = RAAN;
w2 = w+0.001;
anom2 = anom+0.001;
M2 = E2M(anom2E(D2R*anom2,e2),e2);
n2 = sqrt(mu/a2^3);

%% Position and Velocity in ECI

[r_ECI_init, v_ECI_init]=...
    OE2ECI(a, e, inc, RAAN, w, anom, mu);
[r2_ECI_init, v2_ECI_init]=...
    OE2ECI(a2, e2, inc2, RAAN2, w2, anom2, mu);

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

%% Integration Constants (YA)

f = anom*pi/180;
k = 1+e*cos(f);

t = 0;
I = mu^2/h^3*t;

s = k*sin(f);
c = k*cos(f);
dsdf = cos(f) - e + 2*e*cos(f)^2;
dcdf = -sin(f)*(2*e*cos(f) + 1);

phif = [s  c  2-3*e*s*I 0 0 0;...
    dsdf  dcdf  -3*e*(dsdf*I+s/k^2) 0 0 0;...
    c*(1+1/k) -s*(1+1/k) -3*k^2*I 1 0 0;...
    -2*s e-2*c -3*(1-2*e*s*I) 0 0 0;...
    0 0 0 0 cos(f) sin(f);
    0 0 0 0 -sin(f) cos(f)];

fdot = sqrt(mu/(a^3*(1-e^2)^3))*(1+e*cos(f))^2;

xb = rho_RTN_init(1)/r0;
dxb = rhod_RTN_init(1)/r0/fdot;
yb = rho_RTN_init(2)/r0;
dyb = rhod_RTN_init(2)/r0/fdot;
zb = rho_RTN_init(3)/r0;
dzb = rhod_RTN_init(3)/r0/fdot;

consts = phif\[xb;dxb;yb;dyb;zb;dzb];
consts'

%% YA Propagation

orbitCount = 15;
stopTime = orbitCount/revs_per_day*24*60*60;
stepSize = 0.01;

eta = sqrt(1-e^2);
phi_inv_f0 = 1/eta^2*[-3*s*(k+e^2)/k^2, c-2*e, 0, -s*(k+1)/k, 0 , 0;
    -3*(e+c/k), -s, 0, -(c*(k+1)/k+e), 0, 0;
    3*k-eta^2, e*s, 0, k^2, 0, 0;
    -3*e*s*(k+1)/k^2, -2+e*c, eta^2, -e*s*(k+1)/k, 0, 0;
    0, 0, 0, 0, eta^2*cos(f), -eta^2*sin(f);
    0, 0, 0, 0, eta^2*sin(f), eta^2*cos(f)];

x_vec = [];
f0 = f;
for f = f0+stepSize:stepSize:30*pi+stepSize
    func = @(f) 1./((1+e*cos(f)).^2);
    I = integral(func, f0, f);
    
    k = 1+e*cos(f);
    
    s = k*sin(f);
    c = k*cos(f);
    dsdf = cos(f) - e + 2*e*cos(f)^2;
    dcdf = -sin(f)*(2*e*cos(f) + 1);
    
    phif = [s  c  2-3*e*s*I 0 0 0;...
        dsdf  dcdf  -3*e*(dsdf*I+s/k^2) 0 0 0;...
        c*(1+1/k) -s*(1+1/k) -3*k^2*I 1 0 0;...
        -2*s e-2*c -3*(1-2*e*s*I) 0 0 0;...
        0 0 0 0 cos(f) sin(f);
        0 0 0 0 -sin(f) cos(f)];
    
    xb_tot = phif*phi_inv_f0*[xb;dxb;yb;dyb;zb;dzb];
    x_vec = [x_vec xb_tot];
end

x_vec = x_vec'*r0;
x_vec(:,2) = x_vec(:,2)*fdot;
x_vec(:,4) = x_vec(:,4)*fdot;
x_vec(:,6) = x_vec(:,6)*fdot;

fvec = f0+stepSize:stepSize:30*pi+stepSize;
[X] = YApropagation(fvec,consts,f0,e);

%% YA PLOT

% figure
% plot3(x_vec(:,1),x_vec(:,3),x_vec(:,5),...
%     'LineWidth',2)
% xlabel('r_R (km)')
% ylabel('r_T (km)')
% zlabel('r_N (km)')
% title('Relative RTN Position from HCW Equations')
% 
% figure
% plot3(x_vec(:,2),x_vec(:,4),x_vec(:,6),...
%     'LineWidth',2)
% xlabel('v_R (km/s)')
% ylabel('v_T (km/s)')
% zlabel('v_N (km/s)')
% title('Relative RTN Velocity from HCW Equations')
% 
% figure
% subplot(3,1,1)
% plot(x_vec(:,3),x_vec(:,1),...
%     'LineWidth',2)
% xlabel('r_T (km)')
% ylabel('r_R (km)')
% title('Relative TR Position from HCW Equations')
% subplot(3,1,2)
% plot(x_vec(:,5),x_vec(:,1),...
%     'LineWidth',2)
% xlabel('r_N (km)')
% ylabel('r_R (km)')
% title('Relative NR Position from HCW Equations')
% subplot(3,1,3)
% plot(x_vec(:,3),x_vec(:,5),...
%     'LineWidth',2)
% xlabel('r_T (km)')
% ylabel('r_N (km)')
% title('Relative TN Position from HCW Equations')
% 
% figure
% subplot(3,1,1)
% plot(x_vec(:,4),x_vec(:,2),...
%     'LineWidth',2)
% xlabel('v_T (km/s)')
% ylabel('v_R (km/s)')
% title('Relative TR Velocity from HCW Equations')
% subplot(3,1,2)
% plot(x_vec(:,6),x_vec(:,2),...
%     'LineWidth',2)
% xlabel('v_N (km/s)')
% ylabel('v_R (km/s)')
% title('Relative NR Velocity from HCW Equations')
% subplot(3,1,3);
% plot(x_vec(:,4),x_vec(:,6),...
%     'LineWidth',2)
% xlabel('v_T (km/s)')
% ylabel('v_N (km/s)')
% title('Relative TN Velocity from HCW Equations')

%% Orbital Element Differences

del_a = a2-a;
del_inc = (inc2-inc)*pi/180;
del_RAAN = (RAAN2-RAAN)*pi/180;
del_anom2 = (anom2-anom)*pi/180;
del_n = n2-n;
del_M0 = E2M(anom2E(anom2*D2R,e2),e2) - E2M(anom2E(anom*D2R,e),e);

q1_c = e*cos(w*D2R); q1_d = e2*cos(w2*D2R); del_q1 = q1_d-q1_c;
q2_c = e*sin(w*D2R); q2_d = e2*sin(w2*D2R); del_q2 = q2_d-q2_c;

del_e = 1/sqrt(q1_c^2+q2_c^2)*(q1_c*del_q1 + q2_c*del_q2);
del_w = 1/(q1_c^2+q2_c^2)*(q1_c*del_q2 - q2_c*del_q1);

%% Schaub Propagation for Arbitrary Eccentricity

xSh_vec = [];
for f = f0+stepSize:stepSize:30*pi+stepSize
    
    k = 1+e*cos(f);
    
    numOrb = floor((f+0.0000001)/(2*pi));
    M = E2M(anom2E(f,e),e)+numOrb*2*pi;
    del_M = del_M0 - 3/2*del_a/a*(M-M0);
    fdot = sqrt(mu/(a^3*(1-e^2)^3))*(1+e*cos(f))^2;
    
    x = del_a/a + k*e*del_M/eta^3*sin(f)-k*del_e/eta^2*cos(f);
    y = k^2*del_M/eta^3 + del_w + (2+e*cos(f))*del_e/eta^2*sin(f)+...
        cos(inc*pi/180)*del_RAAN;
    z = sin(f+w*pi/180)*del_inc - cos(f+w*pi/180)*sin(inc*pi/180)*del_RAAN;
    
    xdot = (e*cos(f)*...
        (del_M0+(3*del_a*(M0-acos((e + cos(f))/(e*cos(f)+ 1))+...
        e*(1 -(e + cos(f))^2/(e*cos(f)+ 1)^2)^(1/2)))/(2*a))*...
        (e*cos(f)+ 1))/(1 - e^2)^(3/2)-(e^2*sin(f)^2*...
        (del_M0 +(3*del_a*(M0 - acos((e + cos(f))/(e*cos(f)+ 1))+...
        e*(1 -(e + cos(f))^2/(e*cos(f)+ 1)^2)^(1/2)))/(2*a)))/...
        (1 - e^2)^(3/2)-(del_e*e*cos(f)*sin(f))/(e^2 - 1)-(del_e*sin(f)...
        *(e*cos(f)+ 1))/(e^2 - 1)-(3*del_a*e*sin(f)*...
        ((sin(f)/(e*cos(f)+ 1)-(e*sin(f)*(e + cos(f)))/...
        (e*cos(f)+ 1)^2)/(1 -(e + cos(f))^2/(e*cos(f)+ 1)^2)^(1/2)-...
        (e*((2*sin(f)*(e + cos(f)))/(e*cos(f)+ 1)^2 -(2*e*sin(f)*...
        (e + cos(f))^2)/(e*cos(f)+ 1)^3))/(2*(1 -(e + cos(f))^2/...
        (e*cos(f)+ 1)^2)^(1/2)))*(e*cos(f)+ 1))/(2*a*(1 - e^2)^(3/2));

    ydot = (del_e*e*sin(f)^2)/(e^2 - 1)-(del_e*cos(f)*(e*cos(f) + 2))/...
        (e^2 - 1) - (2*e*sin(f)*(del_M0 + (3*del_a*(M0 - acos((e +...
        cos(f))/(e*cos(f) + 1)) + e*(1 - (e + cos(f))^2/(e*cos(f) +...
        1)^2)^(1/2)))/(2*a))*(e*cos(f) + 1))/(1 - e^2)^(3/2) -...
        (3*del_a*((sin(f)/(e*cos(f) + 1) - (e*sin(f)*(e + cos(f)))/...
        (e*cos(f) + 1)^2)/(1 - (e + cos(f))^2/(e*cos(f) + 1)^2)^(1/2) -...
        (e*((2*sin(f)*(e + cos(f)))/(e*cos(f) + 1)^2 - (2*e*sin(f)*...
        (e + cos(f))^2)/(e*cos(f) + 1)^3))/(2*(1 - (e + cos(f))^2/...
        (e*cos(f) + 1)^2)^(1/2)))*(e*cos(f) + 1)^2)/(2*a*(1 - e^2)^(3/2));
    
    zdot = del_inc*cos(f + (pi*w)/180) + del_RAAN*sin(f +...
        (pi*w)/180)*sin((pi*inc)/180);


%     xdot = -e^2*del_M/eta^3*sin(f)+k*e/eta^3*(del_n)/fdot*sin(f)+...
%         k*e*del_M/eta^3*cos(f)+k*del_e/eta^2*sin(f)-...
%         e*del_e/eta^2*sin(f)*cos(f);
%     ydot = -2*k*e*sin(f)*del_M/eta^3 + k^2/eta^3*(del_n)/fdot +...
%         2*del_e/eta^2*cos(f)+ e*del_e/eta^2*(-sin(f)^2+cos(f)^2);
%     zdot = cos(f)*cos(w*pi/180)*del_inc - sin(f)*sin(w*pi/180)*del_inc +...
%         sin(f)*cos(w*pi/180)*sin(inc*pi/180)*del_RAAN +...
%         cos(f)*sin(w*pi/180)*sin(inc*pi/180)*del_RAAN;
    
    xSh_vec = [xSh_vec [x;xdot*fdot;y;ydot*fdot;z;zdot*fdot]];
    
end

xSh_vec = xSh_vec';%*r0;

oe = [a,e,inc,RAAN,w,E2M(anom2E(anom,e),e)];
del_oe = [del_a,del_e,del_inc,del_RAAN,del_w,del_M0];
[X2] = SCHAUBprop(fvec,del_oe,oe,mu);

%% SCHAUB PLOT

% figure
% plot3(xSh_vec(:,1),xSh_vec(:,3),xSh_vec(:,5),...
%     'LineWidth',2)
% xlabel('r_R (km)')
% ylabel('r_T (km)')
% zlabel('r_N (km)')
% title('Relative RTN Position from HCW Equations')
% 
% figure
% plot3(xSh_vec(:,2),xSh_vec(:,4),xSh_vec(:,6),...
%     'LineWidth',2)
% xlabel('v_R (km/s)')
% ylabel('v_T (km/s)')
% zlabel('v_N (km/s)')
% title('Relative RTN Velocity from HCW Equations')
% 
% figure
% subplot(3,1,1)
% plot(xSh_vec(:,3),xSh_vec(:,1),...
%     'LineWidth',2)
% xlabel('r_T (km)')
% ylabel('r_R (km)')
% title('Relative TR Position from HCW Equations')
% subplot(3,1,2)
% plot(xSh_vec(:,5),xSh_vec(:,1),...
%     'LineWidth',2)
% xlabel('r_N (km)')
% ylabel('r_R (km)')
% title('Relative NR Position from HCW Equations')
% subplot(3,1,3)
% plot(xSh_vec(:,3),xSh_vec(:,5),...
%     'LineWidth',2)
% xlabel('r_T (km)')
% ylabel('r_N (km)')
% title('Relative TN Position from HCW Equations')
% 
% figure
% subplot(3,1,1)
% plot(xSh_vec(:,4),xSh_vec(:,2),...
%     'LineWidth',2)
% xlabel('v_T (km/s)')
% ylabel('v_R (km/s)')
% title('Relative TR Velocity from HCW Equations')
% subplot(3,1,2)
% plot(xSh_vec(:,6),xSh_vec(:,2),...
%     'LineWidth',2)
% xlabel('v_N (km/s)')
% ylabel('v_R (km/s)')
% title('Relative NR Velocity from HCW Equations')
% subplot(3,1,3);
% plot(xSh_vec(:,4),xSh_vec(:,6),...
%     'LineWidth',2)
% xlabel('v_T (km/s)')
% ylabel('v_N (km/s)')
% title('Relative TN Velocity from HCW Equations')

%% Orbit Prop

rd0 = sqrt(2*(ENERGY+mu/r0)-thd0*h);

rel_pos0 = [rho_RTN_init;r0;th0];
rel_vel0 = [rhod_RTN_init;rd0;thd0];

options = simset('SrcWorkspace','current');
orbitCount = 15;
stopTime = orbitCount/revs_per_day*24*60*60;
stepSize = 5;
impulseTime = 0;
dV = [0;0;0];
sim('relative_orbit_prop',[],options);

rel_pos_RTN = squeeze(rel_pos_RTN);
rel_vel_RTN = squeeze(rel_vel_RTN);

%% Plotting

X = X'*r0;
X(:,2) = X(:,2)*fdot;
X(:,4) = X(:,4)*fdot;
X(:,6) = X(:,6)*fdot;

X2 = X2'*r0;
X2(:,2) = X2(:,2)*fdot;
X2(:,4) = X2(:,4)*fdot;
X2(:,6) = X2(:,6)*fdot;

% figure
% hold on
% plot3(x_vec(:,1),x_vec(:,3),x_vec(:,5),...
%     'LineWidth',2)
% plot3(xSh_vec(:,1),xSh_vec(:,3),xSh_vec(:,5),...
%     'LineWidth',2)
% plot3(rel_pos_RTN(1,:),rel_pos_RTN(2,:),rel_pos_RTN(3,:),...
%     'LineWidth',2)
% xlabel('r_R (km)')
% ylabel('r_T (km)')
% zlabel('r_N (km)')
% title('Relative RTN Position from HCW Equations')
% hold off
% legend('YA','Schaub','Prop')
% 
% figure
% hold on
% plot3(x_vec(:,2),x_vec(:,4),x_vec(:,6),...
%     'LineWidth',2)
% plot3(xSh_vec(:,2),xSh_vec(:,4),xSh_vec(:,6),...
%     'LineWidth',2)
% plot3(rel_vel_RTN(1,:),rel_vel_RTN(2,:),rel_vel_RTN(3,:),...
%     'LineWidth',2)
% xlabel('v_R (km/s)')
% ylabel('v_T (km/s)')
% zlabel('v_N (km/s)')
% title('Relative RTN Velocity from HCW Equations')
% hold off
% legend('YA','Schaub','Prop')
% 
% figure
% subplot(3,1,1)
% hold on
% plot(x_vec(:,3),x_vec(:,1),...
%     'LineWidth',2)
% plot(xSh_vec(:,3),xSh_vec(:,1),...
%     'LineWidth',2)
% plot(rel_pos_RTN(2,:),rel_pos_RTN(1,:),...
%     'LineWidth',2)
% hold off
% xlabel('r_T (km)')
% ylabel('r_R (km)')
% title('Relative TR Position from HCW Equations')
% subplot(3,1,2)
% hold on
% plot(x_vec(:,5),x_vec(:,1),...
%     'LineWidth',2)
% plot(xSh_vec(:,5),xSh_vec(:,1),...
%     'LineWidth',2)
% plot(rel_pos_RTN(3,:),rel_pos_RTN(1,:),...
%     'LineWidth',2)
% hold off
% xlabel('r_N (km)')
% ylabel('r_R (km)')
% title('Relative NR Position from HCW Equations')
% subplot(3,1,3)
% hold on
% plot(x_vec(:,3),x_vec(:,5),...
%     'LineWidth',2)
% plot(xSh_vec(:,3),xSh_vec(:,5),...
%     'LineWidth',2)
% plot(rel_pos_RTN(2,:),rel_pos_RTN(3,:),...
%     'LineWidth',2)
% hold off
% xlabel('r_T (km)')
% ylabel('r_N (km)')
% title('Relative TN Position from HCW Equations')
% legend('YA','Schaub','Prop')

% figure
% subplot(3,1,1)
% hold on
% plot(x_vec(:,4),x_vec(:,2),...
%     'LineWidth',2)
% plot(xSh_vec(:,4),xSh_vec(:,2),...
%     'LineWidth',2)
% plot(rel_vel_RTN(2,:),rel_vel_RTN(1,:),...
%     'LineWidth',2)
% hold off
% xlabel('v_T (km/s)')
% ylabel('v_R (km/s)')
% title('Relative TR Velocity from HCW Equations')
% subplot(3,1,2)
% hold on
% plot(x_vec(:,6),x_vec(:,2),...
%     'LineWidth',2)
% plot(xSh_vec(:,6),xSh_vec(:,2),...
%     'LineWidth',2)
% plot(rel_vel_RTN(3,:),rel_vel_RTN(1,:),...
%     'LineWidth',2)
% hold off
% xlabel('v_N (km/s)')
% ylabel('v_R (km/s)')
% title('Relative NR Velocity from HCW Equations')
% subplot(3,1,3);
% hold on
% plot(x_vec(:,4),x_vec(:,6),...
%     'LineWidth',2)
% plot(xSh_vec(:,4),xSh_vec(:,6),...
%     'LineWidth',2)
% plot(rel_vel_RTN(2,:),rel_vel_RTN(3,:),...
%     'LineWidth',2)
% hold off
% xlabel('v_T (km/s)')
% ylabel('v_N (km/s)')
% title('Relative TN Velocity from HCW Equations')
% legend('YA','Schaub','Prop')

figure
hold on
plot3(X(:,1),X(:,3),X(:,5),...
    'LineWidth',2)
plot3(X2(:,1),X2(:,3),X2(:,5),...
    'LineWidth',2)
plot3(rel_pos_RTN(1,:),rel_pos_RTN(2,:),rel_pos_RTN(3,:),...
    'LineWidth',2)
xlabel('r_R (km)')
ylabel('r_T (km)')
zlabel('r_N (km)')
title('Relative RTN Position from HCW Equations')
hold off
legend('YA','Schaub','Prop')

figure
hold on
plot3(X(:,2),X(:,4),X(:,6),...
    'LineWidth',2)
plot3(X2(:,2),X2(:,4),X2(:,6),...
    'LineWidth',2)
plot3(rel_vel_RTN(1,:),rel_vel_RTN(2,:),rel_vel_RTN(3,:),...
    'LineWidth',2)
xlabel('v_R (km/s)')
ylabel('v_T (km/s)')
zlabel('v_N (km/s)')
title('Relative RTN Velocity from HCW Equations')
hold off
legend('YA','Schaub','Prop')

figure
subplot(3,1,1)
hold on
plot(X(:,3),X(:,1),...
    'LineWidth',2)
plot(X2(:,3),X2(:,1),...
    'LineWidth',2)
plot(rel_pos_RTN(2,:),rel_pos_RTN(1,:),...
    'LineWidth',2)
hold off
xlabel('r_T (km)')
ylabel('r_R (km)')
title('Relative TR Position from HCW Equations')
subplot(3,1,2)
hold on
plot(X(:,5),X(:,1),...
    'LineWidth',2)
plot(X2(:,5),X2(:,1),...
    'LineWidth',2)
plot(rel_pos_RTN(3,:),rel_pos_RTN(1,:),...
    'LineWidth',2)
hold off
xlabel('r_N (km)')
ylabel('r_R (km)')
title('Relative NR Position from HCW Equations')
subplot(3,1,3)
hold on
plot(X(:,3),X(:,5),...
    'LineWidth',2)
plot(X2(:,3),X2(:,5),...
    'LineWidth',2)
plot(rel_pos_RTN(2,:),rel_pos_RTN(3,:),...
    'LineWidth',2)
hold off
xlabel('r_T (km)')
ylabel('r_N (km)')
title('Relative TN Position from HCW Equations')
legend('YA','Schaub','Prop')

figure
subplot(3,1,1)
hold on
plot(X(:,4),X(:,2),...
    'LineWidth',2)
plot(X2(:,4),X2(:,2),...
    'LineWidth',2)
plot(rel_vel_RTN(2,:),rel_vel_RTN(1,:),...
    'LineWidth',2)
hold off
xlabel('v_T (km/s)')
ylabel('v_R (km/s)')
title('Relative TR Velocity from HCW Equations')
subplot(3,1,2)
hold on
plot(X(:,6),X(:,2),...
    'LineWidth',2)
plot(X2(:,6),X2(:,2),...
    'LineWidth',2)
plot(rel_vel_RTN(3,:),rel_vel_RTN(1,:),...
    'LineWidth',2)
hold off
xlabel('v_N (km/s)')
ylabel('v_R (km/s)')
title('Relative NR Velocity from HCW Equations')
subplot(3,1,3);
hold on
plot(X(:,4),X(:,6),...
    'LineWidth',2)
plot(X2(:,4),X2(:,6),...
    'LineWidth',2)
plot(rel_vel_RTN(2,:),rel_vel_RTN(3,:),...
    'LineWidth',2)
hold off
xlabel('v_T (km/s)')
ylabel('v_N (km/s)')
title('Relative TN Velocity from HCW Equations')
legend('YA','Schaub','Prop')