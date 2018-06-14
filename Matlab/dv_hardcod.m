function [cout,dV_vec] = impulse_man(n,cin,u,uM,dv,dcm,t)

i = cin+1;

if double(t) == 116 && n == 1
    dV_vec = dcm*[dv(1,1);dv(2,1);dv(3,1)];
%     fprintf('Completed Burn No.%0.0f at t = %0.2f\n',i,double(t));
    cout = cin+1;
elseif n == 2 && double(t) == 3519
    dV_vec = dcm*[dv(1,1);dv(2,1);dv(3,1)];
%     fprintf('Completed Burn No.%0.0f at t = %0.2f\n',i,double(t));
    cout = cin+1;
elseif double(t) == 830 && n == 1
    dV_vec = -dcm*[dv(1,2);dv(2,2);dv(3,2)];
%     fprintf('Completed Burn No.%0.0f at t = %0.2f\n',i,double(t));
    cout = cin+1;
elseif  n == 2 && double(t) == 4163
    dV_vec = -dcm*[dv(1,2);dv(2,2);dv(3,2)];
%     fprintf('Completed Burn No.%0.0f at t = %0.2f\n',i,double(t));
    cout = cin+1;
elseif double(t) == 1385 && n == 1
    dV_vec = -dcm*[dv(1,3);dv(2,3);dv(3,3)];
%     fprintf('Completed Burn No.%0.0f at t = %0.2f\n',i,double(t));
    cout = cin+1;
elseif  n == 2 && double(t) == 4619
    dV_vec = -dcm*[dv(1,3);dv(2,3);dv(3,3)];
%     fprintf('Completed Burn No.%0.0f at t = %0.2f\n',i,double(t));
    cout = cin+1;
else
    dV_vec = [0;0;0];
    cout = cin;
end

end