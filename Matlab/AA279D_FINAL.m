
clear;
clc;
close all;

TLE_Reader;

run_ukf = 1;

dv_Flag = 1;
stepSize = 5;

if dv_Flag == 0
    bool_dv = 0;
else
    bool_dv = 1/stepSize;
end

dV_Timer = 5;
Q_LQR = 1/(10)^2*eye(6);
R_LQR = 1/(0.1)^2*eye(3);
N = zeros(6,3);

[~,D] = eig([Q_LQR N;N' R_LQR]);

w = 0;
RAAN = 0;
anom = 0;

mu = 398600;
J2 = 1.081874*10^-3;
Re = 6378137*10^-3;
M0 = M;
tol = 10^-5;
D2R = pi/180;
R2D = 1/D2R;
options = simset('SrcWorkspace','current');

OE = [a,e,inc,RAAN,w,anom];
OE_d =  [a+1000/1000,e,inc+0.001,RAAN+0.001,w+0.01,anom-0.1];
OE_d2 = [a+350/1000,e,inc+0.0001,RAAN+0.0001,w+0.0001,anom-0.0001];
roe = oe2roe(OE, OE_d)
aroe = oe2roe(OE, OE_d)*a*1000;
roe_man = oe2roe(OE, OE_d2);

[r_eci_D, v_eci_D]=...
    OE2ECI(OE_d(1), OE_d(2), OE_d(3), OE_d(4), OE_d(5), OE_d(6), mu);
[r_eci_C, v_eci_C]=...
    OE2ECI(OE(1), OE(2), OE(3), OE(4), OE(5), OE(6), mu);

da = roe(1);
dlam = roe(2);
dex = roe(3);
dey = roe(4);
dix = roe(5);
diy = roe(6);

di = [dix,diy];
dix_man = roe_man(5);
diy_man = roe_man(6);
di_man = [dix_man,diy_man];

DVn = n*a*norm(di_man-di);
uM(1) = atan((diy_man-diy)/(dix_man-dix));

de = [dex,dey];
dex_man = roe_man(3);
dey_man = roe_man(4);
de_man = [dex_man,dey_man];

da_man = roe_man(1);

DVt1 = (n*a/4)*((da_man-da)+norm(de_man-de));
DVt2 = (n*a/4)*((da_man-da)-norm(de_man-de));
uM(2) = atan((dey_man-dey)/(dex_man-dex));
uM(3) = uM(2) + pi;

uM = mod(uM,2*pi);

dV = zeros(3);
dV(3,1) = DVn;
dV(2,2) = DVt1;
dV(2,3) = DVt2;

dv_tot = sum(vecnorm(dV))

startTime = 0.0;
stopTime = 2500*stepSize;
burnCount = 1;

sig0 = 2*eye(6);

sim('FINAL_MODEL',[],options);

rRTN1 = rRTN;
roe_out1 = a*roe_out*1000;
r1 = r;
rd1 = rd;
dV_LQR1 = dV_LQR;
dV1 = dV;

ukf_out1 = ukf_out;
for i = 1:size(sig_out,3)
    sig_out1(i) = norm(chol(sig_out(:,:,i)));
end
roe_man1 = a*1000*roe_man;

res_pre1 = res_pre;
res_post1 = res_post;

%%

burnCount = 2;
OE = OE_C(end,:);
if run_ukf == 1
    OE_d = roe2oe(roe_out(end,:), OE, 1);
else
    OE_d =  OE_D(end,:);
end
OE_d2 = [OE(1)+10/1000,OE(2),OE(3),OE(4),OE(5),OE(6)];
roe = oe2roe(OE, OE_d);
aroe = oe2roe(OE, OE_d)*a*1000;
roe_man = oe2roe(OE, OE_d2);

[r_eci_D, v_eci_D]=...
    OE2ECI(OE_d(1), OE_d(2), OE_d(3), OE_d(4), OE_d(5), OE_d(6), mu);
[r_eci_C, v_eci_C]=...
    OE2ECI(OE(1), OE(2), OE(3), OE(4), OE(5), OE(6), mu);

da = roe(1);
dlam = roe(2);
dex = roe(3);
dey = roe(4);
dix = roe(5);
diy = roe(6);

di = [dix,diy];
dix_man = roe_man(5);
diy_man = roe_man(6);
di_man = [dix_man,diy_man];

DVn = n*a*norm(di_man-di);
uM(1) = atan((diy_man-diy)/(dix_man-dix));

de = [dex,dey];
dex_man = roe_man(3);
dey_man = roe_man(4);
de_man = [dex_man,dey_man];

da_man = roe_man(1);

DVt1 = (n*a/4)*((da_man-da)+norm(de_man-de));
DVt2 = (n*a/4)*((da_man-da)-norm(de_man-de));
uM(2) = atan((dey_man-dey)/(dex_man-dex));
uM(3) = uM(2) + pi;

uM = mod(uM,2*pi);

dV = zeros(3);
dV(3,1) = DVn;
dV(2,2) = DVt1;
dV(2,3) = DVt2;

dv_tot2 = dv_tot + sum(vecnorm(dV))

startTime = stopTime+stepSize;
orbitCount = 15;
stopTime = startTime + 8000*stepSize;

sig0 = sig_out(:,:,end);
roe = roe_out(end,:);
sim('FINAL_MODEL',[],options);

% rRTN = [rRTN1;rRTN];
roe_out2 = [roe_out1;a*roe_out*1000];
roe_man2 = a*1000*roe_man;

iss = [vecnorm(rd1'),vecnorm(rd')];
dragon = [vecnorm(r1'),vecnorm(r')];

%%

close all;
clc;

t = (0:stepSize:(length(iss)-1)*stepSize)/60/60;
% t = 0:(length(iss)-1);
dist = (iss-dragon)*1000;

% burn = [116,830,1385,3513,4411,4966];
burn = [116,830,1385,3519,4130,4673];

ph1 = repmat(roe_man1,1,abs(burn(4)-burn(3)));
ph2 = repmat(roe_man2,1,abs(size(roe_out2,1)-burn(6)));

ukf = [ukf_out1;ukf_out];

for i = 1:size(sig_out,3)
    sig_out2(i) = norm(chol(sig_out(:,:,i)));
end

sig = [sig_out1,sig_out2];

pre = [res_pre1;res_pre]*a*1000;
post = [res_post1;res_post]*a*1000;

ukf = a*100*ukf;

dv = [vecnorm(dV_LQR1'),vecnorm(dV_LQR')];
dv(burn) = dv(burn) + [vecnorm(dV1),vecnorm(dV)];

dv_time = cumsum(dv);

figure
plot(t,dv_time*1000,'k','LineWidth',1)
xlabel('Time (hours)')
ylabel('deltaV [m/s]')

figure
hold on
plot(t,dist+20,'k','LineWidth',1)
plot(t(burn),dist(burn)+20,'ro','LineWidth',1)
hold off
xlabel('Time (hours)')
ylabel('Distance Apart [m]')

% figure
% hold on
% plot3(rd1(600:end,1),rd1(600:end,2),rd1(600:end,3),'c','LineWidth',1)
% plot3(rd(1:2500,1),rd(1:2500,2),rd(1:2500,3),'r','LineWidth',1)
% earthPlot(1)
% axis equal
% hold off
% xlabel('[km]')
% ylabel('[km]')
% zlabel('[km]')
%
figure
hold on
plot3(rRTN1(:,1)*1000,rRTN1(:,2)*1000,rRTN1(:,3)*1000,'r','LineWidth',1)
plot3(rRTN(:,1)*1000,rRTN(:,2)*1000,rRTN(:,3)*1000,'c','LineWidth',1)
hold off
xlabel('R [m]')
ylabel('T [m]')
zlabel('N [m]')

figure
subplot(3,2,1),hold on,plot(t,roe_out2(:,1)),hold off;...
    ylabel('a\delta a [m]'),xlabel('time [hr]')
subplot(3,2,2),hold on,plot(t,roe_out2(:,2)),hold off;...
    ylabel('a\delta \lambda [m]'),xlabel('time [hr]')
subplot(3,2,3),hold on,plot(t,roe_out2(:,3)),hold off;...
    ylabel('a\delta e_x [m]'),xlabel('time [hr]')
subplot(3,2,4),hold on,plot(t,roe_out2(:,4)),hold off;...
    ylabel('a\delta e_y [m]'),xlabel('time [hr]')
subplot(3,2,5),hold on,plot(t,roe_out2(:,5)),hold off;...
    ylabel('a\delta i_x [m]'),xlabel('time [hr]')
subplot(3,2,6),hold on,plot(t,roe_out2(:,6)),hold off;...
    ylabel('a\delta i_y [m]'),xlabel('time [hr]')

upb = roe_out2 + 3*sig';
lb = roe_out2 - 3*sig';

ukf_diff = roe_out2-ukf*10;

figure
subplot(3,2,1),hold on,plot(t,ukf(:,1)*10,'b'),...
    plot(t,upb(:,1),'r'),plot(t,lb(:,1),'r'),...
    hold off;...
    ylabel('a\delta a [m]'),xlabel('time [hr]')
subplot(3,2,2),hold on,plot(t,ukf(:,2)*10,'b'),...
    plot(t,upb(:,2),'r'),plot(t,lb(:,2),'r'),...
    hold off;...
    ylabel('a\delta \lambda [m]'),xlabel('time [hr]')
subplot(3,2,3),hold on,plot(t,ukf(:,3)*10,'b'),...
    plot(t,upb(:,3),'r'),plot(t,lb(:,3),'r'),...
    hold off;...    
    ylabel('a\delta e_x [m]'),xlabel('time [hr]')
subplot(3,2,4),hold on,plot(t,ukf(:,4)*10,'b'),...
    plot(t,upb(:,4),'r'),plot(t,lb(:,4),'r'),...
    hold off;...    
    ylabel('a\delta e_y [m]'),xlabel('time [hr]')
subplot(3,2,5),hold on,plot(t,ukf(:,5)*10,'b'),...
    plot(t,upb(:,5),'r'),plot(t,lb(:,5),'r'),...
    hold off;...    
    ylabel('a\delta i_x [m]'),xlabel('time [hr]')
subplot(3,2,6),hold on,plot(t,ukf(:,6)*10,'b'),...
    plot(t,upb(:,6),'r'),plot(t,lb(:,6),'r'),...
    hold off;...    
    ylabel('a\delta i_y [m]'),xlabel('time [hr]')

figure
subplot(3,2,1),hold on,...
    plot(t(burn(3)+2:burn(4)-1),roe_out2(burn(3)+2:burn(4)-1,1)),...
    plot(t(burn(3)+2:burn(4)-1),ph1(1,3:end),'r'),hold off;...
    ylabel('a\delta a [m]'),xlabel('time [hr]')
subplot(3,2,2),hold on,...
    plot(t(burn(3)+2:burn(4)-1),ph1(2,3:end),'r'),...
    plot(t(burn(3)+2:burn(4)-1),roe_out2(burn(3)+2:burn(4)-1,2)),...
    hold off;...
    ylabel('a\delta \lambda [m]'),xlabel('time [hr]')
subplot(3,2,3),hold on,...
    plot(t(burn(3)+2:burn(4)-1),ph1(3,3:end),'r'),
    plot(t(burn(3)+2:burn(4)-1),roe_out2(burn(3)+2:burn(4)-1,3)),...
    hold off;...
    ylabel('a\delta e_x [m]'),xlabel('time [hr]')
subplot(3,2,4),hold on,...
    plot(t(burn(3)+2:burn(4)-1),roe_out2(burn(3)+2:burn(4)-1,4)),...
    plot(t(burn(3)+2:burn(4)-1),ph1(4,3:end),'r'),hold off;...
    ylabel('a\delta e_y [m]'),xlabel('time [hr]')
subplot(3,2,5),hold on,...
    plot(t(burn(3)+2:burn(4)-1),roe_out2(burn(3)+2:burn(4)-1,5)),...
    plot(t(burn(3)+2:burn(4)-1),ph1(5,3:end),'r'),hold off;...
    ylabel('a\delta i_x [m]'),xlabel('time [hr]')
subplot(3,2,6),hold on,...
    plot(t(burn(3)+2:burn(4)-1),roe_out2(burn(3)+2:burn(4)-1,6)),...
    plot(t(burn(3)+2:burn(4)-1),ph1(6,3:end),'r'),hold off;...
    ylabel('a\delta i_y [m]'),xlabel('time [hr]')

figure
subplot(3,2,1),hold on,...
    plot(t(burn(6)+2:size(roe_out2,1)-1),roe_out2(burn(6)+2:size(roe_out2,1)-1,1)),...
    plot(t(burn(6)+2:size(roe_out2,1)-1),ph2(1,3:end),'r'),hold off;...
    ylabel('a\delta a [m]'),xlabel('time [hr]')
subplot(3,2,2),hold on,...
    plot(t(burn(6)+2:size(roe_out2,1)-1),roe_out2(burn(6)+2:size(roe_out2,1)-1,2)),...
    plot(t(burn(6)+2:size(roe_out2,1)-1),ph2(2,3:end),'r'),hold off;...
    ylabel('a\delta \lambda [m]'),xlabel('time [hr]')
subplot(3,2,3),hold on,...
    plot(t(burn(6)+2:size(roe_out2,1)-1),roe_out2(burn(6)+2:size(roe_out2,1)-1,3)),...
    plot(t(burn(6)+2:size(roe_out2,1)-1),ph2(3,3:end),'r'),hold off;...
    ylabel('a\delta e_x [m]'),xlabel('time [hr]')
subplot(3,2,4),hold on,...
    plot(t(burn(6)+2:size(roe_out2,1)-1),roe_out2(burn(6)+2:size(roe_out2,1)-1,4)),...
    plot(t(burn(6)+2:size(roe_out2,1)-1),ph2(4,3:end),'r'),hold off;...
    ylabel('a\delta e_y [m]'),xlabel('time [hr]')
subplot(3,2,5),hold on,...
    plot(t(burn(6)+2:size(roe_out2,1)-1),roe_out2(burn(6)+2:size(roe_out2,1)-1,5)),...
    plot(t(burn(6)+2:size(roe_out2,1)-1),ph2(5,3:end),'r'),hold off;...
    ylabel('a\delta i_x [m]'),xlabel('time [hr]')
subplot(3,2,6),hold on,...
    plot(t(burn(6)+2:size(roe_out2,1)-1),roe_out2(burn(6)+2:size(roe_out2,1)-1,6)),...
    plot(t(burn(6)+2:size(roe_out2,1)-1),ph2(6,3:end),'r'),hold off;...
    ylabel('a\delta i_y [m]'),xlabel('time [hr]')

figure
hold on
plot3(rRTN1(:,1)*1000,rRTN1(:,2)*1000,rRTN1(:,3)*1000,'r','LineWidth',1)
plot3(rRTN(:,1)*1000,rRTN(:,2)*1000,rRTN(:,3)*1000,'c','LineWidth',1)
hold off
xlabel('R [m]')
ylabel('T [m]')
zlabel('N [m]')

figure
subplot(3,2,1),hold on,plot(t,ukf(:,1)*10-0.35,'r'),plot(t,roe_out2(:,1)),...
    hold off;...
    ylabel('a\delta a [m]'),xlabel('time [hr]')
subplot(3,2,2),hold on,plot(t,ukf(:,2)*10-0.35,'r'),plot(t,roe_out2(:,2)),...
    hold off;...
    ylabel('a\delta \lambda [m]'),xlabel('time [hr]')
subplot(3,2,3),hold on,plot(t,ukf(:,3)*10-0.35,'r'),plot(t,roe_out2(:,3)),...
    hold off;...
    ylabel('a\delta e_x [m]'),xlabel('time [hr]')
subplot(3,2,4),hold on,plot(t,ukf(:,4)*10-0.35,'r'),plot(t,roe_out2(:,4)),...
    hold off;...
    ylabel('a\delta e_y [m]'),xlabel('time [hr]')
subplot(3,2,5),hold on,plot(t,ukf(:,5)*10-0.35,'r'),plot(t,roe_out2(:,5)),...
    hold off;...
    ylabel('a\delta i_x [m]'),xlabel('time [hr]')
subplot(3,2,6),hold on,plot(t,ukf(:,6)*10-0.35,'r'),plot(t,roe_out2(:,6)),...
    hold off;...
    ylabel('a\delta i_y [m]'),xlabel('time [hr]')

figure
subplot(3,2,1),hold on,plot(t,pre(:,1)),plot(t,post(:,1),'r'),hold off;...
    ylabel('a\delta a [m]'),xlabel('time [hr]')
subplot(3,2,2),hold on,plot(t,pre(:,2)),plot(t,post(:,2),'r'),hold off;...
    ylabel('a\delta \lambda [m]'),xlabel('time [hr]')
subplot(3,2,3),hold on,plot(t,pre(:,3)),plot(t,post(:,3),'r'),hold off;...
    ylabel('a\delta e_x [m]'),xlabel('time [hr]')
subplot(3,2,4),hold on,plot(t,pre(:,4)),plot(t,post(:,4),'r'),hold off;...
    ylabel('a\delta e_y [m]'),xlabel('time [hr]')
subplot(3,2,5),hold on,plot(t,pre(:,5)),plot(t,post(:,5),'r'),hold off;...
    ylabel('a\delta i_x [m]'),xlabel('time [hr]')
subplot(3,2,6),hold on,plot(t,pre(:,6)),plot(t,post(:,6),'r'),hold off;...
    ylabel('a\delta i_y [m]'),xlabel('time [hr]')

figure
subplot(3,2,1),hold on,plot(t,ukf_diff(:,1)),...
    hold off;...
    ylabel('a\delta a [m]'),xlabel('time [hr]')
subplot(3,2,2),hold on,plot(t,ukf_diff(:,2)),...
    hold off;...
    ylabel('a\delta \lambda [m]'),xlabel('time [hr]')
subplot(3,2,3),hold on,plot(t,ukf_diff(:,3)),...
    hold off;...    
    ylabel('a\delta e_x [m]'),xlabel('time [hr]')
subplot(3,2,4),hold on,plot(t,ukf_diff(:,4)),...
    hold off;...    
    ylabel('a\delta e_y [m]'),xlabel('time [hr]')
subplot(3,2,5),hold on,plot(t,ukf_diff(:,5)),...
    hold off;...    
    ylabel('a\delta i_x [m]'),xlabel('time [hr]')
subplot(3,2,6),hold on,plot(t,ukf_diff(:,6)),...
    hold off;...    
    ylabel('a\delta i_y [m]'),xlabel('time [hr]')

figure
plot(t,sig)
ylabel('Stanfard Deviation [m]'),xlabel('time [hr]')

