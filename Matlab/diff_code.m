clear;
clc;
close all;

syms f;
syms del_a;
syms del_e;
syms del_w;
syms del_RAAN;
syms del_inc;
syms del_M0;

syms a;
syms e;
syms inc;
syms w;
syms M0;


k = 1+e*cos(f);
eta = sqrt(1-e^2);

ESh = acos((e + cos(f))/(1 + e*cos(f)));
M = ESh - e*sin(ESh);
del_M = del_M0 - 3/2*del_a/a*(M-M0);

x = del_a/a + k*e*del_M/eta^3*sin(f)-k*del_e/eta^2*cos(f);
y = k^2*del_M/eta^3 + del_w + (2+e*cos(f))*del_e/eta^2*sin(f)+...
    cos(inc*pi/180)*del_RAAN;
z = sin(f+w*pi/180)*del_inc - cos(f+w*pi/180)*sin(inc*pi/180)*del_RAAN;

dx = diff(x,f)
dy = diff(y,f)
dz = diff(z,f)