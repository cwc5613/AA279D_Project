function [r_curv, v_curv] = rect2curv(r_rect, v_rect, r0)
% rect2curv Convert between rectilinear and curvilinear RTN states
%
%   This function is designed to handle matrix inputs for the relative
%   state (e.g. an array representing the relative state over time). Note
%   that it will output matrices of column vectors of appropriate size.
% 
% Inputs:
%   r_rect - Deputy relative rectilinear position [km; km; km]
%            Composed of elements [x; y; z] 
%   v_rect - Deputy relative rectilinear velocity [km/s; km/s; km/s]
%            Composed of elements [x_dot; y_dot; z_dot]
%       r0 - Chief absolute distance from central body [km]
%
% Outputs:
%   r_curv - Deputy relative curvilinear position [km; rad; rad]
%            Composed of elements [dr; theta; phi]
%   v_curv - Deputy relative curvilinear velocity [km/s; rad/s; rad/s]
%            Composed of elements [dr_dot; theta_dot; phi_dot]

% Set up the input matrices as a (3 x n) series of column vectors
if size(r_rect, 2) == 3 && size(r_rect, 1) ~= 3
    r_rect = r_rect.';
end

if size(v_rect, 2) == 3 && size(v_rect, 1) ~= 3
    v_rect = v_rect.';
end

if size(r0, 2) == 1 && size(r0, 1) ~= 1
    r0 = r0.';
end

% Check for consistent input dimensions
numR = size(r_rect, 2);
numV = size(v_rect, 2);
numR0 = size(r0, 2);

if ~isequal(numR, numV, numR0)
    error('All inputs must be of the same size');
end

% Unpack rectilinear components
x = r_rect(1, :);
y = r_rect(2, :);
z = r_rect(3, :);

xd = v_rect(1, :);
yd = v_rect(2, :);
zd = v_rect(3, :);

% Compute curvilinear position state elements
dr = sqrt((r0 + x).^2 + y.^2 + z.^2) - r0;
th = atan2(y, r0 + x);
phi = asin(z./(dr + r0));

% Compute curvilinear velocity state elements
drd = ((r0 + x).*xd + y.*yd + z.*zd)./(r0 + dr);
thd = ((r0 + x).*yd - y.*xd)./((r0 + x).^2 + y.^2);
phid = ((r0 + dr).*zd - z.*drd)./((r0 + dr).*sqrt((r0 + dr).^2 - z.^2));

% Re-pack curvilinear position and velocity
r_curv = [dr; th; phi];
v_curv = [drd; thd; phid];

end

