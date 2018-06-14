function [cout,dV_vec] = impulse_man(n,cin,u,uM,dv,dcm,t)

i = cin+1;

dv = dv*0;

if i <= 3
    if(abs(uM(i) - u) <= 10^-2)
        if i == 1
            dV_vec = dcm*[dv(1,i);dv(2,i);dv(3,i)]/3;
            fprintf('3 - Completed Burn No.%0.0f at t = %0.2f\n',i,double(t));
        else
            dV_vec = dcm*[dv(1,i);dv(2,i);dv(3,i)]/5;
            fprintf('5 - Completed Burn No.%0.0f at t = %0.2f\n',i,double(t));
        end
        cout = cin+1;
        return;
    else
        dV_vec = [0;0;0];
        cout = cin;
    end
else
    dV_vec = [0;0;0];
    cout = cin;
end

end