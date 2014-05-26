% heat_equation - Program to solve the diffusion equation
% using the Backward Euler method
clear; help heat_equation; % Clear memory and print header

%* Initialize parameters (time step, grid spacing, etc.)
tau = 1e-3; % Enter time step
N = 100; % Number of grid points
L = 1; % The system extends from (x)=(0) to (x)=(L)
h = L/N;
i = 0:(N-1); % 1:N
x = h/2 + i*h;
y = x;
w = 0.2;
xs = 0.5;
ys = 0.5;
coeff = tau/h^2;
tfinal = .2;
tau = .5*h;
tsteps = ceil(tfinal/tau);
%% Initialize Source function
S = zeros(N); % Set all elements to zero
xExponent = (x'-xs).^2;
yExponent = (y'-ys).^2;
S = exp(-xExponent'/w^2);
% mesh(x,S);
deltaFunction = zeros(N,1);
deltaFunction(round(N/2))=2;
%% * Set up the Laplacian operator matrix
lap = zeros(N);  % Set all elements to zero
for i=2:(N-1)
    lap(i,i-1) = 1;
    lap(i,i) = -2;  % Set interior rows
    lap(i,i+1) = 1;
end
% Boundary conditions
lap(1,1)=-1;
lap(1,2)=1;
lap(N,N)=-1;
lap(N,N-1)=1;
%% * Initialize Q-matrix
Q = S; 
%% * Compute A-matrix (Tn+1)=ATn
dM = eye(N) - coeff*lap; 
%% * Initialize loop and plot variables
max_iter = .25*tsteps;
time = linspace(0,max_iter*tau,max_iter);      % Record time for plots
Qplot(:,1) = Q; % initial value
%% * Loop over desired number of steps 
for iter=2:max_iter
    %* Compute new temperature
    Q = dM\(Q);
    Qplot(:,iter) = Q(:);
end
%%
for iter=max_iter:tsteps
    %* Compute new temperature
    Q = dM\(Q);
    Qplot(:,iter) = Q(:);
end

%% Plot
figure(2);clf;
mesh(time,x,Qplot);
xlabel('t (s)'); ylabel('x (m)');
%%
figure(1);clf;
for i = 1:tsteps
    plot(x,Qplot(:,1));
    hold on;
    plot(x,Qplot(:,i))
    title(sprintf('%g',time(i)));
    hold off;
    pause(0.001)
end
% %% Print Plots
% saveFigurePath = '/Users/kevin/SkyDrive/KTH Work/LaTeX Reports/Heat Equation/Figures/';
% %% Plot 1
% set(figure(2), 'PaperPositionMode', 'auto');
% print('-depsc2', [saveFigurePath ...
%     sprintf('deltaFunctionPlot')]);