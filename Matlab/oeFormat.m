function oe1 = oeFormat(oe0, dir)
% oeFormat Converts format of orbital elements
%
% Inputs:
%   oe0 - initial orbit elements (all angles in degrees)
%   dir - (optional) flag indicating which direction to perform the
%         conversion. Keplerian to quasi-nonsingular (1 - default) or
%         quasi-nonsingular to Keplerian (0)
%
% Outputs:
%   oe1 - converted orbit elements (all angles in degrees)
%
% Notes:
%   Keplerian parameterization - [a; e; i; Om; w; v]
%   Quasi-nonsingular parameterization - [a; u; ex; ey; i; Om]

% Check input arguments
if nargin < 2
    dir = 1;
end

if dir
    % Conversion from the classical Keplerian parameterization to
    % quasi-nonsingular parameterization
    a = oe0(1);
    e = oe0(2);
    i = oe0(3);
    Om = oe0(4);
    w = oe0(5);
    v = oe0(6);
    
    v_rad = deg2rad(v);
    M_rad = E2M(anom2E(v_rad, e), e);
    M = rad2deg(M_rad);
    
    u = mod(w + M, 360);
    ex = e.*cosd(w);
    ey = e.*sind(w);
    
    oe1 = [a; u; ex; ey; i; Om];
else
    % Conversion from the quasi-nonsingular parameterization to the
    % classical Keplerian parameterization
    a = oe0(1);
    u = oe0(2);
    ex = oe0(3);
    ey = oe0(4);
    i = oe0(5);
    Om = oe0(6);
    
    e = sqrt(ex.^2 + ey.^2);
    w = atan2d(ey, ex);
    
    M = mod(u - w, 360);
    
    M_rad = deg2rad(M);
    v_rad = E2anom(M2E(M_rad, e, 1E-6), e);
    v = rad2deg(v_rad);
    
    oe1 = [a; e; i; Om; w; v];
end
    

end

