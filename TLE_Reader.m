%%
% Two-line element set
clc;
clear all;

mu = 398600; %  Standard gravitational parameter for the earth
% TLE file name 
fname = 'DRAGON_CRS14_TLE.txt';
% Open the TLE file and read TLE elements
fid = fopen(fname, 'rb');
L1c = fscanf(fid,'%21c%',1);
L2c = fscanf(fid,'%71c%',1);
L3c = fscanf(fid,'%71c%',1);
fprintf(L1c);
fprintf(L2c);
fprintf([L3c,'\n']);
fclose(fid);
% Open the TLE file and read TLE elements
fid = fopen(fname, 'rb');
L1 = fscanf(fid,'%24c%*s',1);
L2 = fscanf(fid,'%d%6d%*c%5d%*3c%*2f%f%f%5d%*c%*d%5d%*c%*d%d%5d',[1,9]);
L3 = fscanf(fid,'%d%6d%f%f%f%f%f%f%f',[1,8]);
fclose(fid);

epoch = L2(1,4)*24*3600;        % Epoch Date and Julian Date Fraction
Db    = L2(1,5);                % Ballistic Coefficient
inc   = L3(1,3);                % Inclination [deg]
RAAN  = L3(1,4);                % Right Ascension of the Ascending Node [deg]
e     = L3(1,5)/1e7;            % Eccentricity 
w     = L3(1,6);                % Argument of periapsis [deg]
M     = L3(1,7);                % Mean anomaly [deg]
n     = L3(1,8);                % Mean motion [Revs per day]
   
% Orbital elements

revs_per_day = n;
n = revs_per_day*2*pi/(24*3600);
a = (mu/n^2)^(1/3);     % Semi-major axis [km]    

%M = 0;

E = M2E(M*pi/180,e,10^-6); %[rad]
anom = E2anom(E, e); %[rad]

E = E*180/pi; %[deg]
anom = anom*180/pi; %[deg]

% Six orbital elements 
OE = [a e inc RAAN w M];
fprintf('\n a [km]   e        inc [deg]     RAAN [deg]   w[deg]      M [deg]')
fprintf('\n %4.2f  %4.4f   %4.4f       %4.4f     %4.4f    %4.4f', OE);

fprintf('\n\nEpoch: %f', epoch)

fprintf('\n\n---------- End of TLE Import ----------\n')