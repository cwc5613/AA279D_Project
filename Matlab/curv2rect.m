function [r_rect, v_rect] = curv2rect(r_curv, v_curv, r0)
% curv2rect Convert between curvilinear and rectilinear RTN states
%
%   This function is designed to handle matrix inputs for the relative
%   state (e.g. an array representing the relative state over time). Note
%   that it will output matrices of column vectors of appropriate size.
% 
% Inputs:
%   r_curv - Deputy relative curvilinear position [km; rad; rad]
%            Composed of elements [dr; theta; phi]
%   v_curv - Deputy relative curvilinear velocity [km/s; rad/s; rad/s]
%            Composed of elements [dr_dot; theta_dot; phi_dot]
%       r0 - Chief absolute distance from central body [km]
%
% Outputs:
%   r_rect - Deputy relative rectilinear position [km; km; km]
%            Composed of elements [x; y; z] 
%   v_rect - Deputy relative rectilinear velocity [km/s; km/s; km/s]
%            Composed of elements [x_dot; y_dot; z_dot]

% Set up the input matrices as a (3 x n) series of column vectors
if size(r_curv, 2) == 3 && size(r_curv, 1) ~= 3
    r_curv = r_curv.';
end

if size(v_curv, 2) == 3 && size(v_curv, 1) ~= 3
    v_curv = v_curv.';
end

if size(r0, 2) == 1 && size(r0, 1) ~= 1
    r0 = r0.';
end

% Check for consistent input dimensions
numR = size(r_curv, 2);
numV = size(v_curv, 2);
numR0 = size(r0, 2);

if ~isequal(numR, numV, numR0)
    error('All inputs must be of the same size');
end

% Unpack curvilinear components
dr = r_curv(1, :);
th = r_curv(2, :);
phi = r_curv(3, :);

drd= v_curv(1, :);
thd = v_curv(2, :);
phid = v_curv(3, :);

% Compute rectilinear position state elements
x = (r0 + dr).*cos(phi).*cos(th) - r0;
y = (r0 + dr).*cos(phi).*sin(th);
z = (r0 + dr).*sin(phi);

% Compute rectilinear velocity state elements
xd = drd.*cos(th).*cos(phi) - (r0 + dr).*(phid.*sin(phi).*cos(th) + ...
    thd.*cos(phi).*sin(th));
yd = drd.*sin(th).*cos(phi) + (r0 + dr).*(-phid.*sin(phi).*sin(th) + ...
    thd.*cos(phi).*cos(th));
zd = drd.*sin(phi) + (r0 + dr).*phid.*cos(phi);

% Re-pack rectilinear position and velocity
r_rect = [x; y; z];
v_rect = [xd; yd; zd];

end

