function vel = chriscub(orient,chi,cij,rho)

% Christoffel Equation. chriscub calculates the velocities for any cubic crystal,
% as a function of its crystallographic orientation and its elastic
% constants. orient is a vector containg the initial orientation of the crystal
% and it is defined as orient = [omega_0 chi_0 theta_0]. chi is a vector
% containg the azymuthal angles (degrees). cij is a vector containing the
% elastic constants in the following order cij = [c11 c12 c44].
% rho (scalar) is the density of the crystal.
% This function in based on the following article: "A.G. Every 1980 - General
% closed-form expression for acoustic waves in elastically anisotropic
% solids., Physical Review B, volume 22, number 4".

% ========________Enrico Bruschini - Roma 19-20/04/2012________========
% ========________     Last modified on 09/09/2012     ________========

% Splits up the orient and cij vectors into their component
om0 = deg2rad(orient(1)); chi0 = deg2rad(orient(2)); th0 = deg2rad(orient(3));
c11 = cij(1);       c12 = cij(2);       c44 = cij(3);

chi = deg2rad(chi);         % Converts azimuthal angles into radians

% Defines the direction cosines and related expressions
n1 =  cos(om0) .* cos(chi0+chi) .* cos(th0) - sin(chi0+chi) .* sin(th0);
n2 = -cos(om0) .* cos(chi0+chi) .* sin(th0) - sin(chi0+chi) .* cos(th0);
n3 =  sin(om0) .* cos(chi0+chi);

n_sq = n1.^2 + n2.^2 + n3.^2;
P = n1.^2.*n2.^2 + n3.^2.*n2.^2 + n1.^2.*n3.^2;
Q = n1.^2.*n2.^2.*n3.^2;

% Defines some variables and functions to be used in the Christoffel equation
C1 = c11 + 2*c44;       C2 = c11 - c44;     K = c11 - c12 - 2*c44;
T = C1 .* n_sq;
G = C2.^2 .* n_sq.^2 - 3*K*(2*C2 - K).*P;
H = C2.^3 .* n_sq.^3 - 9/2*C2*K.*(2.*C2 - K).*P.*n_sq + 27/2.*K.^2.*(3.*C2 - 2.*K).*Q;

HG = H./G.^1.5;

for i = 1:length(H)
if HG(i) > 1
    HG(i) = 1;
end
end

psi = 1/3.*acos(HG);         % Equation (17)

% Trigonometric form of the Christoffel Equation solved for the velocities

n = length(chi);

for i = 1:n
    for j = 1:3
        vel(i,j) = sqrt((T(i) + 2.*G(i).^0.5 .* cos(psi(i) + 2/3*pi*j))/(3*rho));
    end
end     