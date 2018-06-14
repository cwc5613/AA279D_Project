% AA 279D  HW 4

clear;
clc;
close all;

% A
mu = 3.986e5;       % Grav. Parameter   [km^3/s^2]

% Chief
a = 25000;          % Semi-major axis   [km]
e = 0.7;            % Eccentricity      []
inc = 0.68;         % Inclination       [rad]
RAAN = 6.251;       % RAAN              [rad]
w = 6.261;        % Arg. of Periapsis [rad]
M = 0;              % Mean Anomaly      [rad]


% NEAR ECCENTRIC
a = 7500; e = 0.1; inc = deg2rad(38.961); RAAN = deg2rad(358.1559); w = deg2rad(358.729); anom = 0;
a2 = 7500; e2 = 0.1; inc2 = deg2rad(38.960); RAAN2 = deg2rad(358.1559); w2 = deg2rad(358.73); anom2 = 0.001;

% SEMI MAJOR AXIS CHANGE
% a = 7500; e = 0.1; inc = deg2rad(38.961); RAAN = deg2rad(358.1559); w = deg2rad(358.729); anom = 0;
% a2 = 7501; e2 = 0.1; inc2 = deg2rad(38.960); RAAN2 = deg2rad(358.1559); w2 = deg2rad(358.73); anom2 = 0.001;

% FULL ORBIT
% a = 25000; e = 0.7; inc = deg2rad(38.961); RAAN = deg2rad(358.1559); w = deg2rad(358.729); anom = 0;
% a2 = 24998; e2 = 0.7; inc2 = deg2rad(38.960); RAAN2 = deg2rad(358.1559); w2 = deg2rad(358.73); anom2 = 0.001;

M = E2M(anom2E(anom,e),e);
% % a = 6783; e = 0.1; inc = deg2rad(51.6376); RAAN = deg2rad(291.1196); w = deg2rad(3.3044);
% % a = 7400; e = 0.1; inc = deg2rad(29); RAAN = deg2rad(20); w = deg2rad(50);
n = sqrt(mu/a^3);   % Mean motion       [rad/sec]
Tc = 2*pi/n;        % Orbit period      [sec]

[r_ECI_init,v_ECI_init] = OE2ECI(a, e, rad2deg(inc), rad2deg(RAAN), rad2deg(w), rad2deg(anom), mu);
%
% a2 = a;          % Semi-major axis   [km]
% e2 = e;            % Eccentricity      []
% inc2 = inc - deg2rad(.001);
% Om2 = RAAN;%- deg2rad(.01);
% w2 = w;% + .001; % UPDATE IN TABLE FOR REPORT...
% anom2 = deg2rad(.001);
M2 = E2M(anom2E(anom2,e2),e2);
n2 = sqrt(mu/a2^3);   % Mean motion       [rad/sec]
Td = 2*pi/n2;
[r2_ECI_init,v2_ECI_init] = OE2ECI(a2, e2, rad2deg(inc2), rad2deg(RAAN2), rad2deg(w2), rad2deg(anom2), mu);

%%

% Relative Position
% dcm = ECI2RTN(R, V);
rho0 = (r2_ECI_init - r_ECI_init);
drho0 = (v2_ECI_init - v_ECI_init);

rmag = norm(r_ECI_init);
vmag = norm(v_ECI_init);
hb = cross(r_ECI_init,v_ECI_init);
r_hat = r_ECI_init/norm(r_ECI_init);
n_hat = hb/norm(hb);
t_hat = cross(n_hat,r_hat)/norm(cross(n_hat,r_hat));

dtheta = norm(hb)/rmag^2;

Rot = [r_hat, t_hat, n_hat]';
rho_rtn = Rot*rho0;
drho_rtn = Rot*drho0 - cross([0;0;dtheta],rho_rtn);

linearityheck = norm(rho_rtn)/norm(r_ECI_init);


%% TRUE POSITIONS ...
% % Keplerian Chief
% a = a;e = e; inc = inc; RAAN = RAAN; w = w; n = n;
% tol = 1e-10;
% M0 = M;
% StartTime = 0; numOrb = 15;
% StopTime = numOrb*Tc;
% StepSize = Tc/1000;
% %COMM t_kepC = sim('stilltesting2.slx');
% 
% r_kepC = r_ECI_init;
% v_kepC = v_ECI_init;
% 
% % This creates a vector of mean anomaly at all simulated time steps that is
% % then turned into a true anomaly without wraping to 2pi
% %COMM Mvec = M0+n*t_kepC;
% % for ll = 1:length(Mvec)
% %     fvec_val = E2anom(M2E(mod(Mvec(ll),2*pi),e,1e-10),e);
% %     for nn = 1:numOrb
% %         
% %         if Mvec(ll) >= nn*2*pi-.0001
% %             fvec_val = fvec_val+2*pi;
% %         end
% %     end
% %     fvec(ll) = fvec_val;
% % end
% 
% % do the same but in reverse (used in schaub propagation)
% % numorbs = floor((fvec+.0001)/(2*pi));
% % for ll = 1:length(Mvec)
% %     M2(ll) = E2M(anom2E(mod(fvec(ll),2*pi),e),e) + numorbs(ll)*2*pi;
% % end
% 
% % figure; plot(Mvec);hold on; plot(fvec); hold on; plot(M2)
% 
% % Keplerian Deputy
% a = a2;e = e2; inc = inc2; RAAN = Om2; w = w2; n = n2;
% tol = 1e-10;
% M0 = M2;
% %COMM t_kepD = sim('stilltesting2.slx');
% 
% r_kepD = r2_ECI_init;
% v_kepD = v2_ECI_init;
% 
% dr_xyz = (r_kepD - r_kepC);
% dv_xyz = (v_kepD - v_kepC);
% 
% % Rotate to RTN
% for j = 1:length(r_kepC)
%     
%     rmag = norm(r_kepC(j,:));
%     vmag = norm(v_kepC(j,:));
%     h(j,:) = cross(r_kepC(j,:),v_kepC(j,:));
%     %     e_vec(j,:) = 1/mu*((vmag^2-mu/rmag)*r_kepC(j,:) - dot(r_kepC(j,:),v_kepC(j,:))*v_kepC(j,:));
%     %     spec_e_num(j) = vmag^2/2 - mu/rmag;
%     
%     r_hat = r_kepC(j,:)./norm(r_kepC(j,:));
%     n_hat = h(j,:)/norm(h(j,:));
%     t_hat = cross(n_hat,r_hat)/norm(cross(n_hat,r_hat));
%     
%     R_eci2rtn = [r_hat', t_hat', n_hat']';
%     dtheta = norm(h(j,:))/rmag^2;
%     dr_rtn(j,:) = R_eci2rtn*dr_xyz(j,:)';
%     dv_rtn(j,:) = R_eci2rtn*dv_xyz(j,:)' - cross([0;0;dtheta],dr_rtn(j,:)');
%     
%     fdot = sqrt(mu/(a^3*(1-e^2)^3))*(1+e*cos(fvec(j)))^2;
%     
%     dr_rtn_norm(j,:) = dr_rtn(j,:)/rmag;
%     dv_rtn_norm(j,:) = dv_rtn(j,:)/rmag/fdot;
%     %     dv_rtn_norm(j,:) = R_eci2rtn*dv_xyz(j,:)'/rmag/fdot - cross([0;0;dtheta],dr_rtn_norm(j,:)');
%     %     % Oscullating OEs
%     %     [a_num(j), e_numos(j), i_num(j), Om_num(j), w_num(j), anom_num(j)] = ...
%     %         ECI2OE(r_num(j,:), v_num(j,:), mu);
% end
% 
% %  figure;plot(dr_rtn(:,1));hold on;plot(dr_rtn(:,2));hold on; plot(dr_rtn(:,3))

[r_ECI_init, v_ECI_init]=...
    OE2ECI(a, e, inc, RAAN, w, anom, mu);
[r2_ECI_init, v2_ECI_init]=...
    OE2ECI(a2, e2, inc2, RAAN2, w2, anom2, mu);

[r_ECI_init,v_ECI_init] =...
    OE2ECI(a, e, rad2deg(inc), rad2deg(RAAN), rad2deg(w), rad2deg(anom), mu);
[r2_ECI_init,v2_ECI_init] =...
    OE2ECI(a2, e2, rad2deg(inc2), rad2deg(RAAN2), rad2deg(w2), rad2deg(anom2), mu);

th0 = w + anom;
r_ECI_init = norm(r_ECI_init);
h = norm(cross(r_ECI_init,v_ECI_init));
thd0 = h/(r_ECI_init^2);
ENERGY = norm(v_ECI_init)^2/2-mu/r_ECI_init;

rd0 = sqrt(2*(ENERGY+mu/r_ECI_init)-thd0*h);

T_ECI2RTN = ECI2RTN_rel(r_ECI_init,v_ECI_init);

rho_RTN_init = T_ECI2RTN*(r2_ECI_init-r_ECI_init);
rhod_RTN_init = T_ECI2RTN*(v2_ECI_init-v_ECI_init)-...
    cross([0;0;thd0],rho_RTN_init);

Re = 6378137*10^-3;
del_r = abs(norm(rho_RTN_init))*1000 <= Re

rel_pos0 = [rho_RTN_init;r_ECI_init;th0];
rel_vel0 = [rhod_RTN_init;rd0;thd0];

options = simset('SrcWorkspace','current');
orbitCount = 15;
stopTime = orbitCount/15.542654381087770*24*60*60;
stepSize = 5;
impulseTime = 0;
dV = [0;0;0];
sim('relative_orbit_prop',[],options);

f0 = anom;
fdot = sqrt(mu/(a^3*(1-e^2)^3))*(1+e*cos(f0))^2;

dr_rtn_norm = squeeze(rel_pos_RTN)'/r_ECI_init;
dv_rtn_norm = squeeze(rel_vel_RTN)'/r_ECI_init/fdot;

%%
% Integration Constants
f0 = anom;
fdot = sqrt(mu/(a^3*(1-e^2)^3))*(1+e*cos(f0))^2;
r_0 = norm(r_ECI_init);
Xbar0 = zeros(6,1);
Xbar0([1; 3; 5]) = rho_rtn/r_0;
Xbar0([2; 4; 6]) = drho_rtn/r_0/fdot;

k = 1+e*cos(f0);
c = k*cos(f0);
s = k*sin(f0);
eta = sqrt(1-e^2);
e = e;
Phi_inv = 1/eta^2*[-3*s*(k+e^2)/k^2, c-2*e, 0, -s*(k+1)/k, 0 , 0;
    -3*(e+c/k), -s, 0, -(c*(k+1)/k+e), 0, 0;
    3*k-eta^2, e*s, 0, k^2, 0, 0;
    -3*e*s*(k+1)/k^2, -2+e*c, eta^2, -e*s*(k+1)/k, 0, 0;
    0, 0, 0, 0, eta^2*cos(f0), -eta^2*sin(f0);
    0, 0, 0, 0, eta^2*sin(f0), eta^2*cos(f0)];

C = Phi_inv*Xbar0;

% Propagate YA equations for 15 orbits...
% fvec = 0:pi/50:30*pi;
% k1 = 1+e*cos(fvec);
% fvec = mod(fvec,2*pi); % why is this different
% k2 = 1+e*cos(fvec);
stepSize = 0.01;
fvec = f0+stepSize:stepSize:30*pi+stepSize;
[X] = YApropagation(fvec,C,f0,e);

% % SCHAUB OE ELEMENTS
del_a = a2 - a;
del_inc = inc2 - inc;
del_RAAN = RAAN2 - RAAN;

q1 = e*cos(w); q12 = e2*cos(w2); del_q1 = q12-q1;
q2 = e*sin(w); q22 = e2*sin(w2); del_q2 = q22-q2;

del_e = 1/sqrt(q1^2+q2^2)*(q1*del_q1 + q2*del_q2);
del_w = 1/(q1^2+q2^2)*(q1*del_q2 - q2*del_q1);

del_M0 = E2M(anom2E(anom2,e2),e2) - E2M(anom2E(anom,e),e);

oe = [a,e,inc,RAAN,w,E2M(anom2E(anom,e),e)];
del_oe = [del_a,del_e,del_inc,del_RAAN,del_w,del_M0];


[X2] = SCHAUBprop(fvec,del_oe,oe,mu);



%% RELATIVE POSITION FIGURES
figure;
plot3(X(1,:),X(3,:),X(5,:)); hold on
plot3(X2(1,:),X2(3,:),X2(5,:)); hold on
plot3(dr_rtn_norm(:,1),dr_rtn_norm(:,2),dr_rtn_norm(:,3))
xlabel('Radial Position ');
ylabel('Along Track Position ');
zlabel('Cross Track Position ');

figure;
subplot(3,1,1)
plot(X(3,:),X(1,:)); hold on
plot(X2(3,:),X2(1,:));
xlabel('Along Track Position');
ylabel('Radial Position ');

subplot(3,1,2)
plot(X(5,:),X(1,:)); hold on
plot(X2(5,:),X2(1,:));
xlabel('Cross Track Position');
ylabel('Radial Position ');

subplot(3,1,3)
plot(X(3,:),X(5,:)); hold on
plot(X2(3,:),X2(5,:));
xlabel('Along Track Position');
ylabel('Cross Track Position ');

% RELATIVE VELOCITY FIGURES
figure;
plot3(X(2,:),X(4,:),X(6,:)); hold on
plot3(X2(2,:),X2(4,:),X2(6,:)); hold on;
plot3(dv_rtn_norm(:,1),dv_rtn_norm(:,2),dv_rtn_norm(:,3))
xlabel('Radial Velocity ');
ylabel('Along Track Velocity ');
zlabel('Cross Track Velocity ');

figure;
subplot(3,1,1)
plot(X(4,:),X(2,:)); hold on
plot(X2(4,:),X2(2,:));
xlabel('Along Track Velocity ');
ylabel('Radial Velocity');

subplot(3,1,2)
plot(X(6,:),X(2,:)); hold on
plot(X2(6,:),X2(2,:));
xlabel('Cross Track Velocity ');
ylabel('Radial Velocity');

subplot(3,1,3)
plot(X(4,:),X(6,:)); hold on
plot(X2(4,:),X2(6,:));
xlabel('Along Track Velocity ');
ylabel('Cross Track Velocity ');


% ERRORS
figure;
subplot(3,1,1)
plot(X(3,:)-dr_rtn_norm(:,2)',X(1,:)-dr_rtn_norm(:,1)'); hold on
plot(X2(3,:)-dr_rtn_norm(:,2)',X2(1,:)-dr_rtn_norm(:,1)');
xlabel('Along Track Position');
ylabel('Radial Position');

subplot(3,1,2)
plot(X(5,:)-dr_rtn_norm(:,3)',X(1,:)-dr_rtn_norm(:,1)'); hold on
plot(X2(5,:)-dr_rtn_norm(:,3)',X2(1,:)-dr_rtn_norm(:,1)');
xlabel('Cross Track Position');
ylabel('Radial Position');

subplot(3,1,3)
plot(X(3,:)-dr_rtn_norm(:,2)',X(5,:)-dr_rtn_norm(:,3)'); hold on
plot(X2(3,:)-dr_rtn_norm(:,2)',X2(5,:)-dr_rtn_norm(:,3)');
xlabel('Along Track Position ');
ylabel('Cross Track Position ');

% ERRORS
figure;
subplot(3,1,1)
plot(X(4,:)-dv_rtn_norm(:,2)',X(2,:)-dv_rtn_norm(:,1)'); hold on
plot(X2(4,:)-dv_rtn_norm(:,2)',X2(2,:)-dv_rtn_norm(:,1)');
xlabel('Along Track Velocity');
ylabel('Radial Velocity');

subplot(3,1,2)
plot(X(6,:)-dv_rtn_norm(:,3)',X(2,:)-dv_rtn_norm(:,1)'); hold on
plot(X2(6,:)-dv_rtn_norm(:,3)',X2(2,:)-dv_rtn_norm(:,1)');
xlabel('Cross Track Velocity');
ylabel('Radial Velocity');

subplot(3,1,3)
plot(X(4,:)-dv_rtn_norm(:,2)',X(6,:)-dv_rtn_norm(:,3)'); hold on
plot(X2(4,:)-dv_rtn_norm(:,2)',X2(6,:)-dv_rtn_norm(:,3)');
xlabel('Along Track Velocity ');
ylabel('Cross Track Velocity ');