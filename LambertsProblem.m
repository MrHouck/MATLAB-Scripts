% Code Name: LambertsProblem.m
% Code Description: Solves Lambert's Problem and shows the convergence
% history for the orbital parameter, along with the orbital elements for
% the orbit described by the two position vectors and time of flight
% Author: Nathan Houck
% Email: houckn@my.erau.edu
% Class: AE313 - Space Mechanics
% Date: 11/11/2025



clear; clc;
format long


%% Definitions
r1_vec = [6472.7 -7470.8 -2469.8];              % [m]
r2_vec = [6864.0 5916.0 -5933.0];               % [m]
r1 = norm(r1_vec);                              % [m]
r2 = norm(r2_vec);                              % [m]
timeElapsed = 50*60;                            % [min]
mu = 398600;                                    % [m^3 * s^-2]

dTheta1 = acos(dot(r1_vec, r2_vec)/(r1*r2));    % [rad]
dTheta2 = (2*pi)-dTheta1;                       % [rad] (unused)

%% Derived values
k = r1*r2*(1-cos(dTheta1));                     
l = (r1+r2);
m = r1*r2*(1+cos(dTheta1));
pplus = k/(l+sqrt(2*m));
pminus = k/(l-sqrt(2*m));
clear k l m;

%% Set up iteration process
pprevious = 0.7*pplus + 0.3 * pminus;
pguess = 0.3*pplus + 0.7 * pminus;
[tprevious] = calculateTimeOfFlightLambert(pprevious, r1, r2, dTheta1, mu);
[tnext] = calculateTimeOfFlightLambert(pguess, r1, r2, dTheta1, mu);
dT_previous = tprevious - timeElapsed;
dT_next = tnext - timeElapsed;
numIter = 1;
iterationHistory{numIter} = struct('pprevious', pprevious, 'pguess', pguess, 'dT', dT_next);

%% Iterate to find the lagrange coefficients by minimizing dT
while abs(dT_next) > 1e-9
   pnew = pguess - dT_next * (pguess - pprevious) / (dT_next - dT_previous);
  
   % Update for next iteration
   pprevious = pguess;
   dT_previous = dT_next;
  
   pguess = pnew;
   [tnext, lagrangeCoeff] = calculateTimeOfFlightLambert(pguess, r1, r2, dTheta1, mu);
   dT_next = tnext - timeElapsed;
   numIter = numIter + 1;
   iterationHistory{numIter} = struct('pprevious', pprevious, 'pguess', pguess, 'dT', dT_next);
end

%% Convergence Table Display
Iteration = (1:length(iterationHistory))';
pprev_array = cellfun(@(x) x.pprevious, iterationHistory)';
pguess_array = cellfun(@(x) x.pguess, iterationHistory)';
dT_array = cellfun(@(x) x.dT, iterationHistory)';
convergenceTable = table(Iteration, pprev_array, pguess_array, dT_array, ...
   'VariableNames', {'Iteration', 'p_previous', 'p_guess', 'Residual_dT'});
fprintf('\n=== Convergence History ===\n');
fprintf('Converged in %d iterations\n\n', numIter);
disp(convergenceTable);

%% Find velocity vectors from lagrange coefficients and the position vectors
f = lagrangeCoeff(1);
g = lagrangeCoeff(2);
fdot = lagrangeCoeff(3);
gdot = lagrangeCoeff(4);
v1_vec = (r2_vec - f * r1_vec)/g;
v2_vec = fdot*r1_vec + (gdot/g)*(r2_vec - f*r1_vec);
fprintf("v1: %.3f %.3f %.3f\n", v1_vec);
fprintf("v2: %.3f %.3f %.3f\n", v2_vec);
%% Find orbital elements
coe_1 = findOE(r1_vec, v1_vec, mu);
coe_2 = findOE(r2_vec, v2_vec, mu); % really only care about the true anomaly from here
fprintf("---------------------------------------------------------------\n");
fprintf("h1: %.3f\ne1: %.3f\ni1: %.3f\nΩ1: %.3f\nω1: %.3f\nθ1: %.3f\n", coe_1);
fprintf("---------------------------------------------------------------\n");
fprintf("θ2: %.3f\n", coe_2(6));
% E = eccentric anomaly
% Me = mean anomaly
% ToF = time of flight
%% Find eccentric anomaly and other info to find ToF
e=coe_1(2);
h=coe_1(1);
theta1 = coe_1(6);
theta2 = coe_2(6);
E1 = trueToEccentricAnomaly(theta1, e);
E2 = trueToEccentricAnomaly(theta2, e);
Me1=E1-(e*sin(E1));
Me2=E2-(e*sin(E2));
t1=(h^3/mu^2)*(1/(1-e^2)^(3/2))*Me1;
t2=(h^3/mu^2)*(1/(1-e^2)^(3/2))*Me2;
ToF=abs(t2-t1);
%% Results
fprintf("---------------------------------------------------------------\n");
fprintf("Eccentric Anomaly at Position 1: %.3f\n", E1);
fprintf("Eccentric Anomaly at Position 2: %.3f\n", E2);
fprintf("Mean Eccentric Anomaly at Position 1: %.3f\n", Me1);
fprintf("Mean Eccentric Anomaly at Position 2: %.3f\n", Me2);
fprintf("ToF: %.3f minutes\n", ToF/60);

%%

function [orbitalElements] = findOE(r, v, mu)
    % FINDOE Determine orbital elements from position and velocity vectors
    % 
    % Description:
    %   Calculates the six orbital elements from the position and velocity
    %   vectors.
    % 
    % Inputs:
    %   r   - Position vector           [m]
    %   v   - Velocity vector           [m/s]
    %   mu  - Gravitaional parameter    [m^3 * s^-2]
    %
    % Outputs:
    %   orbitalElements - 1x6 vector
    %                     h     - Specific angular momentum magnitude   [km^2/s]
    %                     e     - Eccentricity                          
    %                     i     - Inclination                           [degrees]
    %                     Omega - Right Ascension of ascending node     [degrees]
    %                     omega - Argument of periapsis                 [degrees]
    %                     theta - True anomaly                          [degrees]
    %

    h_vec = cross(r, v);
    h=norm(h_vec);
    i = acosd(h_vec(3)/h);
    e_vec = ((1/mu)*cross(v, h_vec)) - (r/norm(r));
    e = norm(e_vec);
    N = cross([0 0 1], h_vec);
    norm(N)
    if N(2) >= 0
        Omega = acosd(N(1)/norm(N));
    else
        Omega = 360 - acosd(N(1)/norm(N));
    end
    if e_vec(3) >= 0
        omega = acosd(dot(N, e_vec)/(norm(N) * e));
    else
        omega = 360-acosd(dot(N, e_vec)/(norm(N) * e));
    end
    if dot(r, v) > 0
        theta = acosd(dot(e_vec, r)/(norm(r)*e));
        disp("Case 1")
    else
        disp("Case 2")
        theta = 360 - acosd(dot(e_vec, r)/(norm(r)*e));
    end
    orbitalElements = [h e i Omega omega theta];
end

%%

function [tof, lagrangeCoeff] = calculateTimeOfFlightLambert(p, r1, r2, dTheta, mu)
    % CALCULATETIMEOFFLIGHTLAMBERT Calculate time of flight and lagrange
    %   coefficients
    % 
    % Description: 
    %  Calculates the time of flight between two positions in the orbit and
    %  the associated lagrange coefficients
    %
    % Inputs:
    %   p               - Orbital parameter (h^2/mu)
    %   r1              - Magnitude of initial position vector
    %   r2              - Magnitude of final position vector
    %   dTheta          - Change in true anomaly
    %   mu              - Gravitational parameter
    %
    % Outputs:
    %   tof             - Time of flight between the two positions
    %   lagrangeCoeff   - 1x4 vector [f, g, fdot, gdot]
    %
    
    k = r1*r2*(1-cos(dTheta));                     
    l = (r1+r2);
    m = r1*r2*(1+cos(dTheta));
    a= (m*k*p)/(((2*m-l^2)*p^2)+(2*k*l*p)-k^2);
    f = 1 - (r2/p)*(1-cos(dTheta));
    g = (r1*r2*sin(dTheta))/sqrt(mu*p);
    fdot = sqrt(mu/p)*((1-cos(dTheta))/(sin(dTheta)))*((1-cos(dTheta))/p - 1/r1 - 1/r2);
    gdot = (1 + fdot*g)/f;
    cosdE = 1-((r1/a)*(1-f));
    sindE = -((r1*r2*fdot)/sqrt(mu*a));
    dE = acos(cosdE);
  
    % dE in quad III IV
    if sindE < 0
           dE = 2*pi - dE;
    end
    tof = g + sqrt(a^3/mu)*(dE-sin(dE));
    lagrangeCoeff = [f g fdot gdot];
end

%%

function E = trueToEccentricAnomaly(theta_deg, e)
    % TRUETOECCENTRICANOMALY Convert true anomaly to eccentric anomaly
    % 
    % Description: converts the true anomaly to the eccentric anomaly,
    % resolving any quadrant amibguities. Returns a value in radians from
    % [0, 2pi)
    % 
    % Inputs:  
    %   theta_deg   - True anomaly in degrees
    %   e           - Eccentricity (0 < e < 1)
    %
    % Outputs:
    %   E           - Eccentric anomaly in radians [0, 2pi)
    %

    if e < 0 || e > 1
        error("e out of bounds (0 < e < 1)");
    end
    theta_rad = deg2rad(theta_deg);
    E = 2*atan(sqrt((1-e)/(1+e))*tan(theta_rad/2));
    if theta_deg > 180
        E = E + 2*pi;
    end
    if E < 0
        E = E + 2*pi;
    end
end
