%  schro - Program to solve the Schrodinger equation
%  for a free particle using the Crank-Nicolson scheme
clear all;  help schro;   % Clear memory and print header

%% * Initialize parameters (grid spacing, time step, etc.)
N = 800;
L = 1;              % System extends from x = 0 to L
h = L/N;          % Grid size
x = h*(0:N-1)+h/2;  % Coordinates of grid points
tau = 1e-3;

%% * Set up the Hamiltonian operator matrix
lap = zeros(N);  % Set all elements to zero
coeff = 1/h^2;
for i=2:(N-1)
    lap(i,i-1) = coeff;
    lap(i,i) = -2*coeff;  % Set interior rows
    lap(i,i+1) = coeff;
end
lap(1,1)=-coeff;
lap(1,2)=coeff;
lap(N,N)=-coeff;
lap(N,N-1)=coeff;

%% * Compute the Backward Matrix
dM = eye(N) - tau*lap;

%% * Initialize the source
x0 = 0.5;          % Location of the center of the wavepacket
sigma0 = 0.2;   % Standard deviation of the wavefunction
Q = exp(-(x'-x0).^2/(sigma0^2));
S = Q;


%% * Initialize loop and plot variables
max_iter = .5/tau;
time = linspace(0,max_iter*tau,max_iter);      % Record time for plots
Qplot(:,1) = Q; % initial value
%% * Loop over desired number of steps (wave circles system once)
for iter=2:round(.25/tau)
    %* Compute new temperature
    Q = dM\(Q)+S;
    Qplot(:,iter) = Q(:);
end
for iter=round(.25/tau):max_iter
    %* Compute new temperature
    Q = dM\(Q);
    Qplot(:,iter) = Q(:);
end


%% * Plot probability versus position at various times
figure(2);clf;
mesh(time,x,Qplot);
xlabel('t'); ylabel('x');
%% Print Plots
saveFigurePath = '/Users/kevin/SkyDrive/KTH Work/LaTeX Reports/Heat Equation/Figures/';
%% Plot 1
set(figure(1), 'PaperPositionMode', 'auto');
print('-depsc2', [saveFigurePath ...
    sprintf('SourceFunction')]);