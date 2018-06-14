function [X] = YApropagation(fvec,C,f0,e)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

iter = 0;
I_fun = @(fff) (1+e*cos(fff)).^(-2);
X = zeros(6,length(fvec));
for f = fvec
    iter = iter+1;
k = 1+e*cos(f);
c = k.*cos(f);
cprime = -k.*sin(f) - e*sin(f).*cos(f);
s = k.*sin(f);
sprime = k.*cos(f) - e*sin(f).^2;
% eta = sqrt(1-e.^2);
I = integral(I_fun,f0,f);


Phi = [s, c, 2-3.*e.*s.*I, 0, 0, 0;
    sprime, cprime, -3*e*(sprime.*I+s./k.^2), 0, 0, 0;
    c.*(1+1./k), -s.*(1+1./k), -3.*k^2.*I, 1, 0, 0;
    -2*s, e-2*c, -3*(1-2*e.*s.*I), 0, 0, 0;
    0, 0, 0, 0, cos(f), sin(f);
    0, 0, 0, 0, -sin(f), cos(f)];

X(:,iter) = Phi*C;
end

end

