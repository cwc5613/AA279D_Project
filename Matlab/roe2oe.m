function oe_d = roe2oe(roe, oe_c, anomFlag)
% roe2oe Converts relative to absolute orbital elements
%
% Inputs:
%        roe - quasi-non singular relative orbital elements
%       oe_c - chief orbital elements (all angles in degrees)
%              [a, e, i, Om, w, anom]
%   anomFlag - (optional) flag indicating whether the provided anomalies
%              are true (1 - default) or mean (0). The output set of deputy
%              OEs will match.
%
% Outputs:
%       oe_d - deputy orbital elements (all angles in degrees)
%              [ad, ed, id, Omd, wd, anomd]

% Check for anomaly flag
if nargin < 3
    anomFlag = 1;
end

% Parse out the chief orbit elements
a = oe_c(1);
e = oe_c(2);
i = oe_c(3);
Om = oe_c(4);
w = oe_c(5);
anom = oe_c(6);

% Get chief mean longitude
if anomFlag
    anom_rad = deg2rad(anom);
    M_rad = E2M(anom2E(anom_rad, e), e);
    M = rad2deg(M_rad);
else
    M = anom;
end

u = w + M;

% Parse out the ROEs
da = roe(1);
dl = roe(2);
dex = roe(3);
dey = roe(4);
dix = roe(5);
diy = roe(6);

% Compute the deputy's orbital elements
ad = a.*da + a;
ed = sqrt((dex + e.*cosd(w)).^2 + (dey + e.*sind(w)).^2);
id = dix + i;
Omd = mod(diy./sind(i) + Om, 360);
wd = atan2d(dey + e.*sind(w), dex + e.*cosd(w));
ud = mod(dl - (Omd - Om).*cosd(i) + u, 360);
Md = mod(ud - wd, 360);

if anomFlag
    Md_rad = deg2rad(Md);
    anomd_rad = E2anom(M2E(Md_rad, ed, 1E-6), ed);
    anomd = rad2deg(anomd_rad);
else
    anomd = Md;
end

oe_d = [ad; ed; id; Omd; wd; anomd];

end