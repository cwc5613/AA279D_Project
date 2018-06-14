function [X] = SCHAUBprop(fvec,del_oe,oe,mu)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
a = oe(1);
e = oe(2);
i = oe(3);
Om = oe(4);
w = oe(5);
M0 = oe(6);

del_a = del_oe(1);
del_e = del_oe(2);
del_i = del_oe(3);
del_Om = del_oe(4);
del_w = del_oe(5);
del_M0 = del_oe(6);

n_c = sqrt(mu/a^3);
n_d = sqrt(mu/(a+del_a)^3);

ii = 0;
for f = fvec
    ii = ii + 1;
    k = 1+e*cos(f);
    eta = sqrt(1-e^2);
    numorbs = floor((f+.000001)/(2*pi));
    M = E2M(anom2E(mod(f,2*pi),e),e) + numorbs*2*pi;
    del_M = del_M0 - 3/2*del_a/a*(M-M0);
    
    x(ii) = del_a/a + k*e*del_M/eta^3*sin(f) - k*del_e/eta^2*cos(f);
    y(ii) = k^2*del_M/eta^3 + del_w + (2+e*cos(f))*del_e/eta^2*sin(f) + cos(i)*del_Om;
    z(ii) = sin(f+w)*del_i - cos(f+w)*sin(i)*del_Om;
    
    
    fdot = sqrt(mu/(a^3*(1-e^2)^3))*(1+e*cos(f))^2;
    dx(ii) = -e^2*del_M/eta^3*sin(f)^2 + k*e/eta^3*(n_d-n_c)/fdot*sin(f) + k*e*del_M/eta^3*cos(f) + k*del_e/eta^2*sin(f) - e*del_e/eta^2*sin(f)*cos(f);
    dy(ii) = -2*k*e*sin(f)*del_M/eta^3 + k^2/eta^3*(n_d-n_c)/fdot + 2*del_e/eta^2*cos(f) + e*del_e/eta^2*(-sin(f)^2+cos(f)^2);
    dz(ii) = cos(f)*cos(w)*del_i - sin(f)*sin(w)*del_i + sin(f)*cos(w)*sin(i)*del_Om + cos(f)*sin(w)*sin(i)*del_Om;
end

X = [x;dx;y;dy;z;dz];


end

