%  schro - Program to solve the Schrodinger equation
%  for a free particle using the Crank-Nicolson scheme
clear all;  help schro;   % Clear memory and print header

%% * Initialize parameters (grid spacing, time step, etc.)
N = 800;
L = 1;              % System extends from -L/2 to L/2
h = L/N;          % Grid size
x = h*(0:N-1)+h/2;  % Coordinates  of grid points
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

%% * Compute the Crank-Nicolson matrix
dCN = eye(N) - tau*lap;

%% * Initialize the wavefunction
x0 = 0.5;          % Location of the center of the wavepacket
sigma0 = 0.2;   % Standard deviation of the wavefunction
Q = exp(-(x'-x0).^2/(sigma0^2));
source_function = Q;

%% * Plot the initial wavefunction
figure(1); clf;
plot(x,real(Q),'-',x,imag(Q),'--');
title('Initial wave function');
xlabel('x');  ylabel('\psi(x)'); legend('Real  ','Imag  ');
drawnow;  pause(1);

%% * Initialize loop and plot variables
max_iter = .5/tau;
time = linspace(0,max_iter*tau,max_iter);      % Record time for plots
Qplot(:,1) = Q; % initial value
%% * Loop over desired number of steps (wave circles system once)
for iter=2:round(.25/tau)
    %* Compute new temperature
    Q = dCN\(Q+source_function);
    Qplot(:,iter) = Q(:);
end
for iter=round(.25/tau):max_iter
    %* Compute new temperature
    Q = dCN\(Q);
    Qplot(:,iter) = Q(:);
end


%% * Plot probability versus position at various times
figure(2);clf;
mesh(time,x,Qplot);
xlabel('t'); ylabel('x');
