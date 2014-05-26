%% Euler Method
% projective motion

%% Specify the initial conditions
%
maxSteps = 25; % Maximum number of steps
v = zeros(1,maxSteps);
r = zeros(1,maxSteps);
r(1) = 0; % Initial position
v(1) = 10; % Initial velocity
g = -9.8; % Gravity
h = 0.1; % Choose a time step
t = 0:h:h*(maxSteps-1);
%% Use equations below to compute the new r and v
%
for i = 2:maxSteps
    v(i) = h*g + v(i-1); % Calculate new velocity
    r(i) = h*v(i-1) + r(i-1);
end
%% Plot
plot(t,r)