function roe  = oe2roe(oe_c, oe_d)
% oe2roe Converts absolute to relative orbital elements
%
% Inputs:
%       oe_c - chief orbital elements (all angles in degrees)
%              [a, e, i, Om, w, anom]
%       oe_d - deputy orbital elements (all angles in degrees)
%              [ad, ed, id, Omd, wd, anomd]
%   anomFlag - (optional) flag indicating whether the provided anomalies
%              are true (1 - default) or mean (0). The output set of deputy
%              OEs will match.
%
% Outputs:
%        roe - quasi-non singular relative orbital elements

% Parse out the chief and deputy orbit elements
a = oe_c(1);
e = oe_c(2);
i = oe_c(3);
Om = oe_c(4);
w = oe_c(5);
anom = oe_c(6);

ad = oe_d(1);
ed = oe_d(2);
id = oe_d(3);
Omd = oe_d(4);
wd = oe_d(5);
anomd = oe_d(6);

M = rad2deg(E2M(anom2E(deg2rad(anom), e), e));
Md = rad2deg(E2M(anom2E(deg2rad(anomd), ed), ed));

u = w + M;
ud = wd + Md;

% Compute D'Amico's quasi-nonsingular ROEs
da = (ad - a)./a;
dl = (ud - u) + (Omd - Om).*cosd(i);
dex = ed.*cosd(wd) - e.*cosd(w);
dey = ed.*sind(wd) - e.*sind(w);
dix = id - i;
diy = (Omd - Om).*sind(i);

if dl > 355
    dl = dl-360;
end

roe = [da; deg2rad(dl); deg2rad(dex); deg2rad(dey); deg2rad(dix); deg2rad(diy)];

end

