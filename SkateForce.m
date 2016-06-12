%Code adapted from @sampson by Alex Goldstein and Adam Shaw

% Simulates a pair of long, rectangular airskates
function [thrustForce, massFlow] = SkateForce(gapHeight,internalPressure,skateLength)
% Skate parameters
k = 8e-8*(2.54e-2)^2;		% Air permeability [m^2]
D = 0.187*2.54e-2;			% Thickness of porous layer [m]
P0 = 11e3;					% Default pressure of the skate [Pa]
n = 2;						% Number of skates
W = 0.3048;					% Skate width [m]
L = skateLength;			% Skate length [m]
T = 400;					% Nominal temperature of air through skate [K]


if internalPressure > 0
    P0 = internalPressure;
end


%Assign used variable name to gapHeight
H=gapHeight;

% Physical constants
M_air = 28.97e-3;			% Molecular weight of air [kg/mol]
R = 8.3144598;				% Molar gas constant [J/K*mol]
g = 9.81;					% Acceleration of gravity [m/s^2]
c = 343;					% Speed of sound [m/s]

% Intermediate calculations
mu = 0.01827e-3*(291.15+120)/(T+120)... % Viscosity of air [Pa*s]
	*(T/291.15)^(1.5);
alpha =	sqrt(12*k/(H^3*D));	% Dimensionless parameter "alpha"
A = n*L*W;						% Total skate area [m^2]

% Force as a function of gap height [N]
thrustForce = n*L*P0*(W-2/alpha*tanh(alpha*W/2));

% Flow rate as a function of gap height [kg/s]
massFlow = (P0*M_air/(R*T))*...
	W*P0*alpha/(2*mu)*tanh(alpha*W/2)*H^3*(1/2-1/3)*(2*(W+L));
