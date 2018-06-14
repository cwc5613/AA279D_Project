% AA 279D  HW 4

clear;
clc;
close all;

% A
mu = 3.986e5;       % Grav. Parameter   [km^3/s^2]

% Chief
a_c = 25000;          % Semi-major axis   [km]
e_c = 0.7;            % Eccentricity      []
i_c = 0.68;         % Inclination       [rad]
Om_c = 6.251;       % RAAN              [rad]
w_c = 6.261;        % Arg. of Periapsis [rad]
M_c = 0;              % Mean Anomaly      [rad]

type = 3;
plotsol = 0;
plot_rtnerr = 0;
plotcomp = 0;

if type == 1
    % NEAR ECCENTRIC
    a_c = 6782.423; e_c = 0.15; i_c = deg2rad(51.6438); Om_c = deg2rad(331.1221); w_c = deg2rad(355.8915); anom_c = 0;
    a_d = 6782.423; e_d = 0.15; i_d = deg2rad(51.6428); Om_d = deg2rad(331.1221); w_d = deg2rad(355.8925); anom_d = 0.001;
elseif type == 2
    % SEMI MAJOR AXIS CHANGE
    a_c = 6782.423; e_c = 0.1; i_c = deg2rad(51.6438); Om_c = deg2rad(331.1221); w_c = deg2rad(355.8915); anom_c = 0;
    a_d = 6783.423; e_d = 0.1; i_d = deg2rad(51.6428); Om_d = deg2rad(331.1221); w_d = deg2rad(355.8925); anom_d = 0.001;
elseif type == 3
    % FULL ORBIT
    a_c = 6782.423; e_c = 0.7; i_c = deg2rad(51.6438); Om_c = deg2rad(331.1221); w_c = deg2rad(355.8915); anom_c = 0;
    a_d = 6782.423; e_d = 0.7; i_d = deg2rad(51.6428); Om_d = deg2rad(331.1221); w_d = deg2rad(355.8925); anom_d = 0.001;
end

M_c = E2M(anom2E(anom_c,e_c),e_c);
% % a_c = 6783; e_c = 0.1; i_c = deg2rad(51.6376); Om_c = deg2rad(291.1196); w_c = deg2rad(3.3044);
% % a_c = 7400; e_c = 0.1; i_c = deg2rad(29); Om_c = deg2rad(20); w_c = deg2rad(50);
n_C = sqrt(mu/a_c^3);   % Mean motion       [rad/sec]
Tc = 2*pi/n_C;        % Orbit period      [sec]

[r0_c,v0_c] = OE2ECI(a_c, e_c, rad2deg(i_c), rad2deg(Om_c), rad2deg(w_c), rad2deg(anom_c), mu);
%
% a_d = a_c;          % Semi-major axis   [km]
% e_d = e_c;            % Eccentricity      []
% i_d = i_c - deg2rad(.001);
% Om_d = Om_c;%- deg2rad(.01);
% w_d = w_c;% + .001; % UPDATE IN TABLE FOR REPORT...
% anom_d = deg2rad(.001);
M_d = E2M(anom2E(anom_d,e_d),e_d);
n_D = sqrt(mu/a_d^3);   % Mean motion       [rad/sec]
Td = 2*pi/n_D;
[r0_d,v0_d] = OE2ECI(a_d, e_d, rad2deg(i_d), rad2deg(Om_d), rad2deg(w_d), rad2deg(anom_d), mu);

% Relative Position
% dcm = ECI2RTN(R, V);
rho0 = (r0_d - r0_c);
drho0 = (v0_d - v0_c);

rmag = norm(r0_c);
vmag = norm(v0_c);
hb = cross(r0_c,v0_c);
r_hat = r0_c/norm(r0_c);
n_hat = hb/norm(hb);
t_hat = cross(n_hat,r_hat)/norm(cross(n_hat,r_hat));

dtheta = norm(hb)/rmag^2;

Rot = [r_hat, t_hat, n_hat]';
rho_rtn = Rot*rho0;
drho_rtn = Rot*drho0 - cross([0;0;dtheta],rho_rtn);

linearity_check = norm(rho_rtn)/norm(r0_c);

Re = 6378137*10^-3;

del_r = abs(norm(rho_rtn))*1000 <= Re

%% TRUE POSITIONS ...
% Keplerian Chief
a = a_c;e = e_c; i = i_c; Om = Om_c; w = w_c; n = n_C;
tol = 1e-10;
M0 = M_c;
StartTime = 0; numOrb = 15;
StopTime = numOrb*Tc;
StepSize = Tc/1000;
t_kepC = sim('stilltesting2.slx');

r_kepC = r_eci;
v_kepC = v_eci;

% This creates a vector of mean anomaly at all simulated time steps that is
% then turned into a true anomaly without wraping to 2pi
Mvec = M0+n*t_kepC;
for ll = 1:length(Mvec)
    fvec_val = E2anom(M2E(mod(Mvec(ll),2*pi),e,1e-10),e);
    for nn = 1:numOrb
        
        if Mvec(ll) >= nn*2*pi-.0001
            fvec_val = fvec_val+2*pi;
        end
    end
    fvec(ll) = fvec_val;
end

% do the same but in reverse (used in schaub propagation)
numorbs = floor((fvec+.0001)/(2*pi));
for ll = 1:length(Mvec)
    M2(ll) = E2M(anom2E(mod(fvec(ll),2*pi),e),e) + numorbs(ll)*2*pi;
end

% figure; plot(Mvec);hold on; plot(fvec); hold on; plot(M2)

% Keplerian Deputy
a = a_d;e = e_d; i = i_d; Om = Om_d; w = w_d; n = n_D;
tol = 1e-10;
M0 = M_d;
t_kepD = sim('stilltesting2.slx');

r_kepD = r_eci;
v_kepD = v_eci;

dr_xyz = (r_kepD - r_kepC);
dv_xyz = (v_kepD - v_kepC);

% Rotate to RTN
for j = 1:length(r_kepC)
    
    rmag = norm(r_kepC(j,:));
    vmag = norm(v_kepC(j,:));
    h(j,:) = cross(r_kepC(j,:),v_kepC(j,:));
    %     e_vec(j,:) = 1/mu*((vmag^2-mu/rmag)*r_kepC(j,:) - dot(r_kepC(j,:),v_kepC(j,:))*v_kepC(j,:));
    %     spec_e_num(j) = vmag^2/2 - mu/rmag;
    
    r_hat = r_kepC(j,:)./norm(r_kepC(j,:));
    n_hat = h(j,:)/norm(h(j,:));
    t_hat = cross(n_hat,r_hat)/norm(cross(n_hat,r_hat));
    
    R_eci2rtn = [r_hat', t_hat', n_hat']';
    dtheta = norm(h(j,:))/rmag^2;
    dr_rtn(j,:) = R_eci2rtn*dr_xyz(j,:)';
    dv_rtn(j,:) = R_eci2rtn*dv_xyz(j,:)' - cross([0;0;dtheta],dr_rtn(j,:)');
    
    fdot = sqrt(mu/(a_c^3*(1-e_c^2)^3))*(1+e_c*cos(fvec(j)))^2;
    
    dr_rtn_norm(j,:) = dr_rtn(j,:)/rmag;
    dv_rtn_norm(j,:) = dv_rtn(j,:)/rmag/fdot;
    %     dv_rtn_norm(j,:) = R_eci2rtn*dv_xyz(j,:)'/rmag/fdot - cross([0;0;dtheta],dr_rtn_norm(j,:)');
    %     % Oscullating OEs
    %     [a_num(j), e_numos(j), i_num(j), Om_num(j), w_num(j), anom_num(j)] = ...
    %         ECI2OE(r_num(j,:), v_num(j,:), mu);
end

%  figure;plot(dr_rtn(:,1));hold on;plot(dr_rtn(:,2));hold on; plot(dr_rtn(:,3))
%%
% Integration Constants
f0 = anom_c;

fdot = sqrt(mu/(a_c^3*(1-e_c^2)^3))*(1+e_c*cos(f0))^2;
r_0 = norm(r0_c);
Xbar0 = zeros(6,1);
Xbar0([1; 3; 5]) = rho_rtn/r_0;
Xbar0([2; 4; 6]) = drho_rtn/r_0/fdot;

k = 1+e_c*cos(f0);
c = k*cos(f0);
s = k*sin(f0);
eta = sqrt(1-e_c^2);
e = e_c;
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
[X] = YApropagation(fvec,C,f0,e);

% % SCHAUB OE ELEMENTS
del_a = a_d - a_c;
del_i = i_d - i_c;
del_Om = Om_d - Om_c;

q1_c = e_c*cos(w_c); q1_d = e_d*cos(w_d); del_q1 = q1_d-q1_c;
q2_c = e_c*sin(w_c); q2_d = e_d*sin(w_d); del_q2 = q2_d-q2_c;

del_e = 1/sqrt(q1_c^2+q2_c^2)*(q1_c*del_q1 + q2_c*del_q2);
del_w = 1/(q1_c^2+q2_c^2)*(q1_c*del_q2 - q2_c*del_q1);

del_M0 = E2M(anom2E(anom_d,e_d),e_d) - E2M(anom2E(anom_c,e_c),e_c);

oe = [a_c,e_c,i_c,Om_c,w_c,E2M(anom2E(anom_c,e_c),e_c)];
del_oe = [del_a,del_e,del_i,del_Om,del_w,del_M0];


[X2] = SCHAUBprop(fvec,del_oe,oe,mu);

%% YA PLOTS

if plotsol == 1
    
    figure;
    plot3(X(1,:),X(3,:),X(5,:),'-','LineWidth',2);
    xlabel('Radial Position ');
    ylabel('Along Track Position ');
    zlabel('Cross Track Position ');
    
    figure;
    plot3(X(2,:),X(4,:),X(6,:),'-','LineWidth',2);
    xlabel('Radial Velocity ');
    ylabel('Along Track Velocity ');
    zlabel('Cross Track Velocity ');
    
    figure;
    subplot(3,1,1)
    plot(X(3,:),X(1,:),'-','LineWidth',2);
    xlabel('Along Track Position');
    ylabel('Radial Position ');
    subplot(3,1,2)
    plot(X(5,:),X(1,:),'-','LineWidth',2);
    xlabel('Cross Track Position');
    ylabel('Radial Position ');
    subplot(3,1,3)
    plot(X(3,:),X(5,:),'-','LineWidth',2);
    xlabel('Along Track Position');
    ylabel('Cross Track Position ');
    
    figure;
    subplot(3,1,1)
    plot(X(4,:),X(2,:),'-','LineWidth',2);
    xlabel('Along Track Velocity ');
    ylabel('Radial Velocity');
    subplot(3,1,2)
    plot(X(6,:),X(2,:),'-','LineWidth',2);
    xlabel('Cross Track Velocity ');
    ylabel('Radial Velocity');
    subplot(3,1,3)
    plot(X(4,:),X(6,:),'-','LineWidth',2);
    xlabel('Along Track Velocity ');
    ylabel('Cross Track Velocity ');
    
    %% SCHAUB PLOTS
    
    figure;
    plot3(X2(1,:),X2(3,:),X2(5,:),'-','LineWidth',2);
    xlabel('Radial Position ');
    ylabel('Along Track Position ');
    zlabel('Cross Track Position ');
    
    figure;
    plot3(X2(2,:),X2(4,:),X2(6,:),'-','LineWidth',2);
    xlabel('Radial Velocity ');
    ylabel('Along Track Velocity ');
    zlabel('Cross Track Velocity ');
    
    figure;
    subplot(3,1,1)
    plot(X2(3,:),X2(1,:),'-','LineWidth',2);
    xlabel('Along Track Position');
    ylabel('Radial Position ');
    subplot(3,1,2)
    plot(X2(5,:),X2(1,:),'-','LineWidth',2);
    xlabel('Cross Track Position');
    ylabel('Radial Position ');
    subplot(3,1,3)
    plot(X2(3,:),X2(5,:),'-','LineWidth',2);
    xlabel('Along Track Position');
    ylabel('Cross Track Position ');
    
    figure;
    subplot(3,1,1)
    plot(X2(4,:),X2(2,:),'-','LineWidth',2);
    xlabel('Along Track Velocity ');
    ylabel('Radial Velocity');
    subplot(3,1,2)
    plot(X2(6,:),X2(2,:),'-','LineWidth',2);
    xlabel('Cross Track Velocity ');
    ylabel('Radial Velocity');
    subplot(3,1,3)
    plot(X2(4,:),X2(6,:),'-','LineWidth',2);
    xlabel('Along Track Velocity ');
    ylabel('Cross Track Velocity ');
    
end

%% RELATIVE POSITION FIGURES

if plotcomp == 1
    
    figure;
    plot3(X(1,:),X(3,:),X(5,:),'-','LineWidth',2); hold on
    plot3(X2(1,:),X2(3,:),X2(5,:),'-','LineWidth',2); hold on
    plot3(dr_rtn_norm(:,1),dr_rtn_norm(:,2),dr_rtn_norm(:,3),'-','LineWidth',2)
    xlabel('Radial Position ');
    ylabel('Along Track Position ');
    zlabel('Cross Track Position ');
    legend('YA','Schaub','Prop')
    
    figure;
    subplot(3,1,1)
    plot(X(3,:),X(1,:),'-','LineWidth',2); hold on
    plot(X2(3,:),X2(1,:),'-','LineWidth',3); hold on
    plot(dr_rtn_norm(:,2),dr_rtn_norm(:,1),'-','LineWidth',1);
    xlabel('Along Track Position');
    ylabel('Radial Position ');
    subplot(3,1,2)
    plot(X(5,:),X(1,:),'-','LineWidth',2); hold on
    plot(X2(5,:),X2(1,:),'-','LineWidth',3); hold on
    plot(dr_rtn_norm(:,3),dr_rtn_norm(:,1),'-','LineWidth',1);
    xlabel('Cross Track Position');
    ylabel('Radial Position ');
    legend('YA','Schaub','Prop')
    subplot(3,1,3)
    plot(X(3,:),X(5,:),'-','LineWidth',2); hold on
    plot(X2(3,:),X2(5,:),'-','LineWidth',3); hold on
    plot(dr_rtn_norm(:,2),dr_rtn_norm(:,3),'-','LineWidth',1);
    xlabel('Along Track Position');
    ylabel('Cross Track Position ');
    
    %% RELATIVE VELOCITY FIGURES
    
    figure;
    plot3(X(2,:),X(4,:),X(6,:),'-','LineWidth',2); hold on
    plot3(X2(2,:),X2(4,:),X2(6,:),'-','LineWidth',2); hold on;
    plot3(dv_rtn_norm(:,1),dv_rtn_norm(:,2),dv_rtn_norm(:,3),'-','LineWidth',2);
    xlabel('Radial Velocity ');
    ylabel('Along Track Velocity ');
    zlabel('Cross Track Velocity ');
    legend('YA','Schaub','Prop')
    
    figure;
    subplot(3,1,1)
    plot(X(4,:),X(2,:),'-','LineWidth',2); hold on
    plot(X2(4,:),X2(2,:),'-','LineWidth',2); hold on
    plot(dv_rtn_norm(:,2),dv_rtn_norm(:,1),'-','LineWidth',2);
    xlabel('Along Track Velocity ');
    ylabel('Radial Velocity');
    subplot(3,1,2)
    plot(X(6,:),X(2,:),'-','LineWidth',2); hold on
    plot(X2(6,:),X2(2,:),'-','LineWidth',2); hold on
    plot(dv_rtn_norm(:,3),dv_rtn_norm(:,1),'-','LineWidth',2);
    xlabel('Cross Track Velocity ');
    ylabel('Radial Velocity');
    legend('YA','Schaub','Prop')
    subplot(3,1,3)
    plot(X(4,:),X(6,:),'-','LineWidth',2); hold on
    plot(X2(4,:),X2(6,:),'-','LineWidth',2); hold on
    plot(dv_rtn_norm(:,2),dv_rtn_norm(:,3),'-','LineWidth',2);
    xlabel('Along Track Velocity ');
    ylabel('Cross Track Velocity ');
    
end
%% Position Errors

if plot_rtnerr == 1
    
    figure;
    subplot(3,1,1)
    plot(X(3,:)-dr_rtn_norm(:,2)',X(1,:)-dr_rtn_norm(:,1)','-','LineWidth',2); hold on
    plot(X2(3,:)-dr_rtn_norm(:,2)',X2(1,:)-dr_rtn_norm(:,1)','-','LineWidth',2);
    xlabel('Along Track Position');
    ylabel('Radial Position');
    title('Relative Position Errors')
    subplot(3,1,2)
    plot(X(5,:)-dr_rtn_norm(:,3)',X(1,:)-dr_rtn_norm(:,1)','-','LineWidth',2); hold on
    plot(X2(5,:)-dr_rtn_norm(:,3)',X2(1,:)-dr_rtn_norm(:,1)','-','LineWidth',2);
    xlabel('Cross Track Position');
    ylabel('Radial Position');
    legend('YA','Schaub')
    subplot(3,1,3)
    plot(X(3,:)-dr_rtn_norm(:,2)',X(5,:)-dr_rtn_norm(:,3)','-','LineWidth',2); hold on
    plot(X2(3,:)-dr_rtn_norm(:,2)',X2(5,:)-dr_rtn_norm(:,3)','-','LineWidth',2);
    xlabel('Along Track Position ');
    ylabel('Cross Track Position ');
    
    %% Velocity Errors
    
    figure;
    subplot(3,1,1)
    plot(X(4,:)-dv_rtn_norm(:,2)',X(2,:)-dv_rtn_norm(:,1)','-','LineWidth',2); hold on
    plot(X2(4,:)-dv_rtn_norm(:,2)',X2(2,:)-dv_rtn_norm(:,1)','-','LineWidth',2);
    xlabel('Along Track Velocity');
    ylabel('Radial Velocity');
    title('Relative Velocity Errors')
    subplot(3,1,2)
    plot(X(6,:)-dv_rtn_norm(:,3)',X(2,:)-dv_rtn_norm(:,1)','-','LineWidth',2); hold on
    plot(X2(6,:)-dv_rtn_norm(:,3)',X2(2,:)-dv_rtn_norm(:,1)','-','LineWidth',2);
    xlabel('Cross Track Velocity');
    ylabel('Radial Velocity');
    legend('YA','Schaub')
    subplot(3,1,3)
    plot(X(4,:)-dv_rtn_norm(:,2)',X(6,:)-dv_rtn_norm(:,3)','-','LineWidth',2); hold on
    plot(X2(4,:)-dv_rtn_norm(:,2)',X2(6,:)-dv_rtn_norm(:,3)','-','LineWidth',2);
    xlabel('Along Track Velocity ');
    ylabel('Cross Track Velocity ');
    
end

%% TIME ERRORS

figure;
subplot(3,1,1)
plot(X(1,:)-dr_rtn_norm(:,1)','-','LineWidth',2); hold on
plot(X2(1,:)-dr_rtn_norm(:,1)','-','LineWidth',2);
ylabel('Radial Position');
xlabel('Time (sec)');
title('Relative Position Errors')
subplot(3,1,2)
plot(X(3,:)-dr_rtn_norm(:,2)','-','LineWidth',2); hold on
plot(X2(3,:)-dr_rtn_norm(:,2)','-','LineWidth',2);
xlabel('Time (sec)');
ylabel('Along Track Position');
legend('YA','Schaub')
subplot(3,1,3)
plot(X(5,:)-dr_rtn_norm(:,3)','-','LineWidth',2); hold on
plot(X2(5,:)-dr_rtn_norm(:,3)','-','LineWidth',2);
xlabel('Time (sec)');
ylabel('Cross Track Position ');

figure;
subplot(3,1,1)
plot(X(2,:)-dv_rtn_norm(:,1)','-','LineWidth',2); hold on
plot(X2(2,:)-dv_rtn_norm(:,1)','-','LineWidth',2);
xlabel('Time (sec)');
ylabel('Radial Velocity');
title('Relative Velocity Errors')
subplot(3,1,2)
plot(X(4,:)-dv_rtn_norm(:,2)','-','LineWidth',2); hold on
plot(X2(4,:)-dv_rtn_norm(:,2)','-','LineWidth',2);
xlabel('Time (sec)');
ylabel('Along Track Velocity');
legend('YA','Schaub')
subplot(3,1,3)
plot(X(6,:)-dv_rtn_norm(:,3)','-','LineWidth',2); hold on
plot(X2(6,:)-dv_rtn_norm(:,3)','-','LineWidth',2);
xlabel('Time (sec)');
ylabel('Cross Track Velocity ');