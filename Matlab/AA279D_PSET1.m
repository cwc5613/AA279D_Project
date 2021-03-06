clear;
clc;
close all;

TLE_Reader;

mu = 398600;
J2 = 1.081874*10^-3;
Re = 6378137*10^-3;
M0 = M;
tol = 10^-5;

[r_eci, v_eci] = OE2ECI(a, e, inc, RAAN, w, anom, mu);

initCond = [r_eci, v_eci];

f0 = 0;
v0 = v_eci;
r0 = r_eci;

options = simset('SrcWorkspace','current');
orbitCount = 10;
stopTime = orbitCount/revs_per_day*24*60*60;
%stepSize = 24*60*60/revs_per_day/100000;
stepSize = 1;
timeSpan = [0:stepSize:stopTime]';

sim('NUM_PROP',[],options);

out_r_NUM_PROP = reshape(NUM_r,[size(NUM_r,1),size(NUM_r,3)]);
out_v_NUM_PROP = reshape(NUM_v,[size(NUM_v,1),size(NUM_v,3)]);

% figure
% hold on
% plot3(out_r_NUM_PROP(1,:),out_r_NUM_PROP(2,:),out_r_NUM_PROP(3,:))
% earthPlot(1)
% hold off

sim('J2_Accel',[],options);

out_r_J2_Accel = reshape(J2_r,[size(J2_r,1),size(J2_r,3)]);
out_v_J2_Accel = reshape(J2_v,[size(J2_v,1),size(J2_v,3)]);

figure
hold on
plot3(out_r_J2_Accel(1,:),out_r_J2_Accel(2,:),out_r_J2_Accel(3,:))
earthPlot(1)
hold off

sim('Orbit_Sim_ECEF',[],options);

v_eci2(1,:) = v_eci;

r_eci2 = r_eci2';
v_eci2 = v_eci2';

err_r = [out_r_NUM_PROP(1,:)-r_eci2(1,:);...
    out_r_NUM_PROP(2,:)-r_eci2(2,:);...
    out_r_NUM_PROP(3,:)-r_eci2(3,:)];

err_v = [out_v_NUM_PROP(1,:)-v_eci2(1,:);...
    out_v_NUM_PROP(2,:)-v_eci2(2,:);...
    out_v_NUM_PROP(3,:)-v_eci2(3,:)];

for i = 1:size(err_r,2)
    [r_rtn1(:,i),v_rtn1(:,i)] = ECI2RTN(out_r_NUM_PROP(:,i), out_v_NUM_PROP(:,i));
    [r_rtn2(:,i),v_rtn2(:,i)] = ECI2RTN(r_eci2(:,i), v_eci2(:,i));
    osc_h(:,i) = norm(cross(out_r_NUM_PROP(:,i),out_v_NUM_PROP(:,i)));
    osc_h_J2(:,i) = norm(cross(out_r_J2_Accel(:,i),out_v_J2_Accel(:,i)));
end

% err_r_rtn = [r_rtn1(1,:)-r_rtn2(1,:);...
%     r_rtn1(2,:)-r_rtn2(2,:);...
%     r_rtn1(3,:)-r_rtn2(3,:)];
%
% err_v_rtn = [v_rtn1(1,:)-v_rtn2(1,:);...
%     v_rtn1(2,:)-v_rtn2(2,:);...
%     v_rtn1(3,:)-v_rtn2(3,:)];


% figure('color',[1 1 1]);
% subplot(3,1,1); plot(err_r(1,:)); title(''); ylabel('Pos Err (km)');
% subplot(3,1,2); plot(err_r(2,:)); title(''); ylabel('Pos Err (km)');
% subplot(3,1,3); plot(err_r(3,:)); title(''); ylabel('Pos Err (km)');
%
% figure('color',[1 1 1]);
% subplot(3,1,1); plot(err_v(1,:)); title(''); ylabel('Vel Err (km)');
% subplot(3,1,2); plot(err_v(2,:)); title(''); ylabel('Vel Err (km)');
% subplot(3,1,3); plot(err_v(3,:)); title(''); ylabel('Vel Err (km)');
%
% r_eci = r_eci';
% v_eci = v_eci';
% r_rtn = r_rtn';
% v_rtn = v_rtn';

% figure
% hold on
% plot3(out_r_NUM_PROP(1,:),out_r_NUM_PROP(2,:),out_r_NUM_PROP(3,:))
% plot3(r_eci2(1,:),r_eci2(2,:),r_eci2(3,:))
% earthPlot(1)
% hold off

% figure('color',[1 1 1]);
% subplot(3,1,1); plot(err_r_rtn(1,:)); title(''); ylabel('Radial Position Error (km)');
% subplot(3,1,2); plot(err_r_rtn(2,:)); title(''); ylabel('Tangential Position Error (km)');
% subplot(3,1,3); plot(err_r_rtn(3,:)); title(''); ylabel('Normal Position Error (km)');
%
% figure('color',[1 1 1]);
% subplot(3,1,1); plot(err_v_rtn(1,:)); title(''); ylabel('Radial Velocity Error (km/s)');
% subplot(3,1,2); plot(err_v_rtn(2,:)); title(''); ylabel('Tangential Velocity Error (km/s)');
% subplot(3,1,3); plot(err_v_rtn(3,:)); title(''); ylabel('Normal Velocity Error (km/s)');

%%

a_bar = mean(osc_a_J2);
e_bar = mean(osc_e_J2);
inc_bar = mean(osc_inc_J2);
RAAN_bar = mean(osc_RAAN_J2);
w_bar = mean(osc_w_J2);
anom_bar = mean(osc_anom_J2);

sim('MeanOE17',[],options);

R2D = 180/pi;

mean_a = a + mean_da.*timeSpan;
mean_e = e + mean_de.*timeSpan;
mean_inc = inc + R2D*mean_dinc.*timeSpan;
mean_RAAN = RAAN + R2D*mean_dRAAN.*timeSpan;
mean_w = w + R2D*mean_dw.*timeSpan;
mean_anom = anom + R2D*mean_danom.*timeSpan;

for i =1:size(err_r,2)
    osc_evec(:,i) = ((norm(out_v_NUM_PROP(:,i))^2/mu)-(1/norm(out_r_NUM_PROP(:,i))))*...
        out_r_NUM_PROP(:,i)-(dot(out_r_NUM_PROP(:,i),out_v_NUM_PROP(:,i)))*out_v_NUM_PROP(:,i)/mu;
    osc_evec_J2(:,i) = ((norm(out_v_J2_Accel(:,i))^2/mu)-(1/norm(out_r_J2_Accel(:,i))))*...
        out_r_J2_Accel(:,i)-(dot(out_r_J2_Accel(:,i),out_v_J2_Accel(:,i)))*out_v_J2_Accel(:,i)/mu;
    osc_hvec(:,i) = cross(out_r_NUM_PROP(:,i),out_v_NUM_PROP(:,i));
    osc_hvec_J2(:,i) = cross(out_r_J2_Accel(:,i),out_v_J2_Accel(:,i));
    
end

%%

figure
subplot(6,1,1);hold on; plot(osc_a_J2);plot(osc_a);xlabel('Time (sec)');ylabel('Semimajor Axis (km)');hold off;
subplot(6,1,2);hold on; plot(osc_e_J2);plot(osc_e);xlabel('Time (sec)');ylabel('Eccentricity');hold off;
subplot(6,1,3);hold on; plot(osc_inc_J2);plot(osc_inc);xlabel('Time (sec)');ylabel('Inclination (\circ)');hold off;
subplot(6,1,4);hold on; plot(osc_RAAN_J2);plot(osc_RAAN);xlabel('Time (sec)');ylabel('RAAN (\circ)');hold off;
subplot(6,1,5);hold on; plot(osc_w_J2);plot(osc_w);xlabel('Time (sec)');ylabel('Arg. of Periapsis (\circ)');hold off;
subplot(6,1,6);hold on; plot(osc_anom_J2);plot(osc_anom);xlabel('Time (sec)');ylabel('True Anomaly (\circ)');hold off;

figure
title('Osculating Orbital Elements with J2 Effects')
subplot(3,1,1);hold on; plot(osc_evec(1,:));plot(osc_evec_J2(1,:));xlabel('Time (sec)');ylabel('e_x');hold off;
subplot(3,1,2);hold on; plot(osc_evec(2,:));plot(osc_evec_J2(2,:));xlabel('Time (sec)');ylabel('e_y');hold off;
subplot(3,1,3);hold on; plot(osc_evec(3,:));plot(osc_evec_J2(3,:));xlabel('Time (sec)');ylabel('e_z');hold off;

figure
title('Osculating Orbital Elements with J2 Effects')
subplot(4,1,1);hold on; plot(osc_hvec(1,:));plot(osc_hvec_J2(1,:));xlabel('Time (sec)');ylabel('h_x (km^2/s)');hold off;
subplot(4,1,2);hold on; plot(osc_hvec(2,:));plot(osc_hvec_J2(2,:));xlabel('Time (sec)');ylabel('h_y (km^2/s)');hold off;
subplot(4,1,3);hold on; plot(osc_hvec(3,:));plot(osc_hvec_J2(3,:));xlabel('Time (sec)');ylabel('h_z (km^2/s)');hold off;
subplot(4,1,4);hold on; plot(osc_ME);plot(osc_ME_J2);xlabel('Time (sec)');ylabel('\epsilon (km^2/s^2)');hold off;

figure
subplot(6,1,1);hold on; plot(osc_a_J2);plot(mean_a);xlabel('Time (sec)');ylabel('Semimajor Axis (km)');hold off;
subplot(6,1,2);hold on; plot(osc_e_J2);plot(mean_e);xlabel('Time (sec)');ylabel('Eccentricity');hold off;
subplot(6,1,3);hold on; plot(osc_inc_J2);plot(mean_inc);xlabel('Time (sec)');ylabel('Inclination (\circ)');hold off;
subplot(6,1,4);hold on; plot(osc_RAAN_J2);plot(mean_RAAN);xlabel('Time (sec)');ylabel('RAAN (\circ)');hold off;
subplot(6,1,5);hold on; plot(osc_w_J2);plot(mean_w);xlabel('Time (sec)');ylabel('Arg. of Periapsis (\circ)');hold off;
subplot(6,1,6);hold on; plot(osc_anom_J2);plot(mean_anom);xlabel('Time (sec)');ylabel('True Anomaly (\circ)');hold off;

% figure
% title('Unperturbed Osculating Orbital Elements')
% subplot(4,2,1);plot(osc_a);xlabel('Time (sec)');ylabel('Semimajor Axis (km)');
% subplot(4,2,2);plot(osc_e);xlabel('Time (sec)');ylabel('Eccentricity');
% subplot(4,2,3);plot(osc_inc);xlabel('Time (sec)');ylabel('Inclination (\circ)');
% subplot(4,2,4);plot(osc_RAAN);xlabel('Time (sec)');ylabel('RAAN (\circ)');
% subplot(4,2,5);plot(osc_w);xlabel('Time (sec)');ylabel('Arg. of Periapsis (\circ)');
% subplot(4,2,6);plot(osc_anom);xlabel('Time (sec)');ylabel('True Anomaly (\circ)');
% subplot(4,2,7);plot(osc_h);xlabel('Time (sec)');ylabel('Angular Momentum (km^2/s)');
% subplot(4,2,8);plot(osc_ME);xlabel('Time (sec)');ylabel('Specific Mech. Energy (km^2/s^2)');
%
% figure
% title('Osculating Orbital Elements with J2 Effects')
% subplot(3,1,1);plot(osc_evec_J2(1,:));xlabel('Time (sec)');ylabel('Semimajor Axis (km)');
% subplot(3,1,2);plot(osc_evec_J2(2,:));xlabel('Time (sec)');ylabel('Semimajor Axis (km)');
% subplot(3,1,3);plot(osc_evec_J2(3,:));xlabel('Time (sec)');ylabel('Semimajor Axis (km)');
%
% figure
% title('Osculating Orbital Elements with J2 Effects')
% subplot(3,1,1);plot(osc_hvec_J2(1,2:end));xlabel('Time (sec)');ylabel('Semimajor Axis (km)');
% subplot(3,1,2);plot(osc_hvec_J2(2,2:end));xlabel('Time (sec)');ylabel('Semimajor Axis (km)');
% subplot(3,1,3);plot(osc_hvec_J2(3,2:end));xlabel('Time (sec)');ylabel('Semimajor Axis (km)');
%
% figure
% title('Unperturbed Osculating Orbital Elements')
% subplot(4,2,1);plot(osc_a);xlabel('Time (sec)');ylabel('Semimajor Axis (km)');
% subplot(4,2,2);plot(osc_e);xlabel('Time (sec)');ylabel('Eccentricity');
% subplot(4,2,3);plot(osc_inc);xlabel('Time (sec)');ylabel('Inclination (\circ)');
% subplot(4,2,4);plot(osc_RAAN);xlabel('Time (sec)');ylabel('RAAN (\circ)');
% subplot(4,2,5);plot(osc_w);xlabel('Time (sec)');ylabel('Arg. of Periapsis (\circ)');
% subplot(4,2,6);plot(osc_anom);xlabel('Time (sec)');ylabel('True Anomaly (\circ)');
% subplot(4,2,7);plot(osc_h);xlabel('Time (sec)');ylabel('Angular Momentum (km^2/s)');
% subplot(4,2,8);plot(osc_ME);xlabel('Time (sec)');ylabel('Specific Mech. Energy (km^2/s^2)');
%
% figure
% subplot(4,2,1);hold on; plot(osc_a_J2);plot(mean_a);xlabel('Time (sec)');ylabel('Semimajor Axis (km)');hold off;
% subplot(4,2,2);hold on; plot(osc_e_J2);plot(mean_e);xlabel('Time (sec)');ylabel('Eccentricity');hold off;
% subplot(4,2,3);hold on; plot(osc_inc_J2);plot(mean_inc);xlabel('Time (sec)');ylabel('Inclination (\circ)');hold off;
% subplot(4,2,4);hold on; plot(osc_RAAN_J2);plot(mean_RAAN);xlabel('Time (sec)');ylabel('RAAN (\circ)');hold off;
% subplot(4,2,5);hold on; plot(osc_w_J2);plot(mean_w);xlabel('Time (sec)');ylabel('Arg. of Periapsis (\circ)');hold off;
% subplot(4,2,6);hold on; plot(osc_anom_J2);plot(mean_anom);xlabel('Time (sec)');ylabel('True Anomaly (\circ)');hold off;
%
