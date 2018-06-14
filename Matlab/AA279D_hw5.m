%% PROBLEM SET FIVE
% DESCRIPTIVE TEXT

clear;
clc;
close all;

%% (a)

TLE_Reader;

mu = 398600;
J2 = 1.081874*10^-3;
Re = 6378137*10^-3;
M0 = M;
tol = 10^-5;
D2R = pi/180;
R2D = 1/D2R;

%% (b)

% a = 7500;
% e = 0.001;
% inc = 0.68*R2D;
% RAAN = 6.251*R2D;
% w = 6.261*R2D;
% M = 3.409*R2D;
% M0 = M;
% 
% anom = E2anom(M2E(M*D2R,e,10e-7),e);
% anom = anom*R2D;

OE = [a,e,inc,RAAN,w,anom];

roe = [0,100,50,100,30,200]/a*10e-4;
oe_d_init = roe2oe(roe, OE, 1);

a2 = oe_d_init(1);
e2 = oe_d_init(2);
inc2 = oe_d_init(3);
RAAN2 = oe_d_init(4);
w2 = oe_d_init(5);
anom2 = oe_d_init(6);

%% (c)

[r_ECI_init, v_ECI_init]=...
    OE2ECI(a, e, inc, RAAN, w, anom, mu);
[r2_ECI_init, v2_ECI_init]=...
    OE2ECI(a2, e2, inc2, RAAN2, w2, anom2, mu);

f0 = 0;
options = simset('SrcWorkspace','current');
orbitCount = 20;
stopTime = orbitCount/revs_per_day*24*60*60;
stepSize = 5;

%(C - i) Osculating Chief OE - No J2
v0 = v_ECI_init;
r0 = r_ECI_init;
sim('NUM_PROP',[],options);
oe_c = [osc_a, osc_e, osc_inc, osc_RAAN, osc_w, osc_anom, osc_ang];
r_c_eci = squeeze(NUM_r);
v_c_eci = squeeze(NUM_v);

%(C - i) Osculating Deputy OE - No J2
v0 = v2_ECI_init;
r0 = r2_ECI_init;
sim('NUM_PROP',[],options);
oe_d = [osc_a, osc_e, osc_inc, osc_RAAN, osc_w, osc_anom, osc_ang];
r_d_eci = squeeze(NUM_r);
v_d_eci = squeeze(NUM_v);

%(C - i) Osculating Chief OE - With J2
v0 = v_ECI_init;
r0 = r_ECI_init;
sim('J2_Accel',[],options);
oe_c_J2 = [osc_a_J2, osc_e_J2, osc_inc_J2, osc_RAAN_J2, osc_w_J2, osc_anom_J2, osc_ang_J2];
r_c_eci_J2 = squeeze(J2_r);
v_c_eci_J2 = squeeze(J2_v);

%(C - i) Osculating Deputy OE - With J2
v0 = v2_ECI_init;
r0 = r2_ECI_init;
sim('J2_Accel',[],options);
oe_d_J2 = [osc_a_J2, osc_e_J2, osc_inc_J2, osc_RAAN_J2, osc_w_J2, osc_anom_J2, osc_ang_J2];
r_d_eci_J2 = squeeze(J2_r);
v_d_eci_J2 = squeeze(J2_v);

for i = 1:size(oe_c,1)
    %(C - ii) Osculating ROE
    roe_no_J2(i,:) = oe2roe(oe_c(i,1:6), oe_d(i,1:6), 1);
    roe_J2(i,:) = oe2roe(oe_c_J2(i,1:6), oe_d_J2(i,1:6), 1);
    
    if roe_no_J2(i,2)<-350
        roe_no_J2(i,2)=roe_no_J2(i,2)+360;
    end
    if roe_no_J2(i,2)>350
        roe_no_J2(i,2)=roe_no_J2(i,2)-360;
    end
    
    if roe_J2(i,2)<-350
        roe_J2(i,2)=roe_J2(i,2)+360;
    end
    if roe_J2(i,2)>350
        roe_J2(i,2)=roe_J2(i,2)-360;
    end
    
    %(C - iii) Mean OE
    oe_c_mean(i,:) = osc2mean(oe_c(i,1:6), Re, 0);
    oe_d_mean(i,:) = osc2mean(oe_d(i,1:6), Re, 0);
    oe_c_mean_J2(i,:) = osc2mean(oe_c_J2(i,1:6), Re, J2);
    oe_d_mean_J2(i,:) = osc2mean(oe_d_J2(i,1:6), Re, J2);
    
    oe_quas_c_mean(i,:) = oeFormat(oe_c_mean(i,1:6));
    oe_quas_d_mean(i,:) = oeFormat(oe_d_mean(i,1:6));
    oe_quas_c_mean_J2(i,:) = oeFormat(oe_c_mean_J2(i,1:6));
    oe_quas_d_mean_J2(i,:) = oeFormat(oe_d_mean_J2(i,1:6));
    
    oe_quas_c(i,:) = oeFormat(oe_c(i,1:6));
    oe_quas_d(i,:) = oeFormat(oe_d(i,1:6));
    oe_quas_c_J2(i,:) = oeFormat(oe_c_J2(i,1:6));
    oe_quas_d_J2(i,:) = oeFormat(oe_d_J2(i,1:6));

    %(C - iv) Mean ROE
    roe_mean_no_J2(i,:) = oe2roe(oe_c_mean(i,1:6), oe_d_mean(i,1:6), 1);
    roe_mean_J2(i,:) = oe2roe(oe_c_mean_J2(i,1:6), oe_d_mean_J2(i,1:6), 1);
    
    if roe_mean_no_J2(i,2)<-350
        roe_mean_no_J2(i,2)=roe_mean_no_J2(i,2)+360;
    end
    if roe_mean_no_J2(i,2)>350
        roe_mean_no_J2(i,2)=roe_mean_no_J2(i,2)-360;
    end
    
    if roe_mean_J2(i,2)<-350
        roe_mean_J2(i,2)=roe_mean_J2(i,2)+360;
    end
    if roe_mean_J2(i,2)>350
        roe_mean_J2(i,2)=roe_mean_J2(i,2)-360;
    end

end

plotting = 1;

if plotting == 1
    
    figure
    subplot(3,2,1)
    plot(oe_quas_c(:,1)); hold on; plot(oe_quas_c_J2(:,1))
    xlabel('Time (sec)')
    ylabel('a (km)')
    subplot(3,2,2)
    plot(oe_quas_c(:,2)); hold on; plot(oe_quas_c_J2(:,2))
    xlabel('Time (sec)')
    ylabel('u (deg)')
    subplot(3,2,3)
    plot(oe_quas_c(:,3)); hold on; plot(oe_quas_c_J2(:,3))
    xlabel('Time (sec)')
    ylabel('e_x')
    subplot(3,2,4)
    plot(oe_quas_c(:,4)); hold on; plot(oe_quas_c_J2(:,4))
    xlabel('Time (sec)')
    ylabel('e_y')
    subplot(3,2,5)
    plot(oe_quas_c(:,5)); hold on; plot(oe_quas_c_J2(:,5))
    xlabel('Time (sec)')
    ylabel('i (deg)')
    subplot(3,2,6)
    plot(oe_quas_c(:,6)); hold on; plot(oe_quas_c_J2(:,6))
    xlabel('Time (sec)')
    ylabel('\Omega (deg)')
    legend('No J2', 'With J2')
    
    figure
    subplot(3,2,1)
    plot(roe_no_J2(:,1)); hold on; plot(roe_J2(:,1))
    xlabel('Time (sec)')
    ylabel('a (km)')
    subplot(3,2,2)
    plot(roe_no_J2(:,2)); hold on; plot(roe_J2(:,2))
    xlabel('Time (sec)')
    ylabel('u (deg)')
    subplot(3,2,3)
    plot(roe_no_J2(:,3)); hold on; plot(roe_J2(:,3))
    xlabel('Time (sec)')
    ylabel('e_x')
    subplot(3,2,4)
    plot(roe_no_J2(:,4)); hold on; plot(roe_J2(:,4))
    xlabel('Time (sec)')
    ylabel('e_y')
    subplot(3,2,5)
    plot(roe_no_J2(:,5)); hold on; plot(roe_J2(:,5))
    xlabel('Time (sec)')
    ylabel('i (deg)')
    subplot(3,2,6)
    plot(roe_no_J2(:,6)); hold on; plot(roe_J2(:,6))
    xlabel('Time (sec)')
    ylabel('\Omega (deg)')
    legend('No J2', 'With J2')
    
    figure
    subplot(3,2,1)
    plot(oe_quas_c_mean(:,1)); hold on; plot(oe_quas_c_mean_J2(:,1))
    xlabel('Time (sec)')
    ylabel('a (km)')
    subplot(3,2,2)
    plot(oe_quas_c_mean(:,2)); hold on; plot(oe_quas_c_mean_J2(:,2))
    xlabel('Time (sec)')
    ylabel('u (deg)')
    subplot(3,2,3)
    plot(oe_quas_c_mean(:,3)); hold on; plot(oe_quas_c_mean_J2(:,3))
    xlabel('Time (sec)')
    ylabel('e_x')
    subplot(3,2,4)
    plot(oe_quas_c_mean(:,4)); hold on; plot(oe_quas_c_mean_J2(:,4))
    xlabel('Time (sec)')
    ylabel('e_y')
    subplot(3,2,5)
    plot(oe_quas_c_mean(:,5)); hold on; plot(oe_quas_c_mean_J2(:,5))
    xlabel('Time (sec)')
    ylabel('i (deg)')
    subplot(3,2,6)
    plot(oe_quas_c_mean(:,6)); hold on; plot(oe_quas_c_mean_J2(:,6))
    xlabel('Time (sec)')
    ylabel('\Omega (deg)')
    legend('No J2', 'With J2')
    
    figure
    subplot(3,2,1)
    plot(roe_mean_no_J2(:,1)); hold on; plot(roe_mean_J2(:,1))
    xlabel('Time (sec)')
    ylabel('a (km)')
    subplot(3,2,2)
    plot(roe_mean_no_J2(:,2)); hold on; plot(roe_mean_J2(:,2))
    xlabel('Time (sec)')
    ylabel('u (deg)')
    subplot(3,2,3)
    plot(roe_mean_no_J2(:,3)); hold on; plot(roe_mean_J2(:,3))
    xlabel('Time (sec)')
    ylabel('e_x')
    subplot(3,2,4)
    plot(roe_mean_no_J2(:,4)); hold on; plot(roe_mean_J2(:,4))
    xlabel('Time (sec)')
    ylabel('e_y')
    subplot(3,2,5)
    plot(roe_mean_no_J2(:,5)); hold on; plot(roe_mean_J2(:,5))
    xlabel('Time (sec)')
    ylabel('i (deg)')
    subplot(3,2,6)
    plot(roe_mean_no_J2(:,6)); hold on; plot(roe_mean_J2(:,6))
    xlabel('Time (sec)')
    ylabel('\Omega (deg)')
    legend('No J2', 'With J2')
    
end

%% (d)

plotting = 1;

for i = 1:size(r_c_eci,2)
    T = ECI2RTN(r_c_eci(:,i), v_c_eci(:,i));
    T_J2 = ECI2RTN(r_c_eci_J2(:,i), v_c_eci_J2(:,i));
    r_c_rtn(:,i) = T*(r_c_eci(:,i)-r_d_eci(:,i));
    r_c_rtn_J2(:,i) = T_J2*(r_c_eci_J2(:,i)-r_d_eci_J2(:,i));
end

r_c_rtn = r_c_rtn.*1000;
r_c_rtn_J2 = r_c_rtn_J2.*1000;

if plotting == 1
    
    figure
    hold on;
    plot3(r_c_rtn(1,:),r_c_rtn(2,:),r_c_rtn(3,:));...
        xlabel('\delta r_R (km)');ylabel('\delta r_T (km)');zlabel('\delta r_N (km)');
    plot3(r_c_rtn_J2(1,:),r_c_rtn_J2(2,:),r_c_rtn_J2(3,:));...
        xlabel('\delta r_R (km)');ylabel('\delta r_T (km)');zlabel('\delta r_N (km)');
    hold off;
    legend('No J2','With J2')
    axis equal;

    
    figure;
    subplot(3,1,1)
    plot(r_c_rtn(2,:),r_c_rtn(1,:),'-','LineWidth',3); hold on
    plot(r_c_rtn_J2(2,:),r_c_rtn_J2(1,:),'-','LineWidth',0.5); hold on
    xlabel('Along Track Position');
    ylabel('Radial Position ');
    axis equal;
    subplot(3,1,2)
    plot(r_c_rtn(3,:),r_c_rtn(1,:),'-','LineWidth',3); hold on
    plot(r_c_rtn_J2(3,:),r_c_rtn_J2(1,:),'-','LineWidth',0.5); hold on
    xlabel('Cross Track Position');
    ylabel('Radial Position ');
    legend('No J2','J2')
    axis equal;
    subplot(3,1,3)
    plot(r_c_rtn(2,:),r_c_rtn(3,:),'-','LineWidth',3); hold on
    plot(r_c_rtn_J2(2,:),r_c_rtn_J2(3,:),'-','LineWidth',0.5); hold on
    xlabel('Along Track Position');
    ylabel('Cross Track Position ');
    axis equal;
    
end

%% (e)

% ROE = (da; dlam; dex; dey; dix; diy)

% roe_no_J2(:,2)=roe_no_J2(:,2)*D2R;roe_no_J2(:,5:6)=roe_no_J2(:,5:6)*D2R;
% roe_J2(:,2)=roe_J2(:,2)*D2R;roe_J2(:,5:6)=roe_J2(:,5:6)*D2R;
% roe_mean_no_J2(:,2)=roe_mean_no_J2(:,2)*D2R;roe_mean_no_J2(:,5:6)=roe_mean_no_J2(:,5:6)*D2R;
% roe_mean_J2(:,2)=roe_mean_J2(:,2)*D2R;roe_mean_J2(:,5:6)=roe_mean_J2(:,5:6)*D2R;

% roe_no_J2 = a*roe_no_J2*1000;
% roe_J2 = a*roe_J2/10e-4;
% roe_mean_no_J2 = a*roe_mean_no_J2/10e-4;
% roe_mean_J2 = a*roe_mean_J2/10e-4;

plotting = 1;

if plotting == 1
    figure
    hold on;
    plot(a*1000*roe_no_J2(:,3),a*1000*roe_no_J2(:,4))
    plot(a*1000*roe_mean_no_J2(:,3),a*1000*roe_mean_no_J2(:,4))
    hold off;
    xlabel('Relative Eccentricity X-Component')
    ylabel('Relative Eccentricity Y-Component')
    title('Unperturbed Relative Eccentricity Vector')
    legend('Osculating','Mean')
    axis equal;
    
    figure
    hold on;
    plot(a*1000*roe_J2(:,3),a*1000*roe_J2(:,4))
    plot(a*1000*roe_mean_J2(:,3),a*1000*roe_mean_J2(:,4))
    hold off;
    xlabel('Relative Eccentricity X-Component')
    ylabel('Relative Eccentricity Y-Component')
    title('J2 Relative Eccentricity Vector')
    legend('Osculating','Mean')
    axis equal;
    
    figure
    hold on;
    plot(a*1000*roe_no_J2(:,5),a*1000*roe_no_J2(:,6),'.')
    plot(a*1000*roe_mean_no_J2(:,5),a*1000*roe_mean_no_J2(:,6),'.')
    hold off;
    xlabel('Relative Inclination X-Component')
    ylabel('Relative Inclination Y-Component')
    title('Unperturbed Relative Inclination Vector')
    legend('Osculating','Mean')
    axis equal;
    
    figure
    hold on;
    plot(a*1000*roe_J2(:,5),a*1000*roe_J2(:,6))
    plot(a*1000*roe_mean_J2(:,5),a*1000*roe_mean_J2(:,6))
    hold off;
    xlabel('Relative Inclination X-Component')
    ylabel('Relative Inclination Y-Component')
    title('J2 Relative Inclination Vector')
    legend('Osculating','Mean')
    axis equal;
    
    figure
    hold on;
    plot(a*1000*roe_no_J2(:,2),a*1000*roe_no_J2(:,1),'o')
    plot(a*1000*roe_mean_no_J2(:,2),a*1000*roe_mean_no_J2(:,1))
    hold off;
    xlabel('a\delta \lambda')
    ylabel('a\delta a')
    title('Unperturbed Relative Mean Longitude and Semi-Major Axis')
    legend('Osculating','Mean')
    axis equal;
    
    figure
    hold on;
    plot(a*1000*roe_J2(:,2),a*1000*roe_J2(:,1),'o')
    plot(a*1000*roe_mean_J2(:,2),a*1000*roe_mean_J2(:,1))
    hold off;
    xlabel('a\delta \lambda')
    ylabel('a\delta a')
    title('J2 Relative Mean Longitude and Semi-Major Axis')
    legend('Osculating','Mean')
    axis equal;
end

