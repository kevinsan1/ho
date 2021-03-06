% heat_equation - Program to solve the diffusion equation
% using the Backward Euler method
clear; help heat_equation; % Clear memory and print header

%* Initialize parameters (time step, grid spacing, etc.)
tau = 0.1; % Enter time step
N = 15; % Number of grid points
L = 2; % The system extends from (x,y)=(0,0) to (x,y)=(L,L)
h = L/N;
i = 0:(N-1);
j = 0:(N-1);
x = h/2 + i*h;
y = h/2 + j*h;
w = 0.2;
xs = 1;
ys = 1;
%% Initialize Source function
S = zeros(N); % Set all elements to zero
xExponent = (x'-xs).^2;
yExponent = (y-ys).^2;
S = exp(-xExponent/w^2)*exp(-yExponent/w^2);
%% Check lengths and if Smax = 1;
if h*length(x) == L && h*length(y) == L 
    disp('Lengths are correct')
else 
    disp('Warning lengths')
end
if round(max(max(S))) == 1
    fprintf('max(S) is correct = %1d \n',round(max(max(S))));
else
    disp(fprintf('Warning S max = %1d \n', max(max(S))));
end
%% * Set up the Laplacian operator matrix
lap = zeros(N);  % Set all elements to zero
coeff = -1/(h^2);
for i=2:(N-1)
    for j=2:(N-1)
        lap(i,i-1,j) = coeff;
        lap(i,i+1,j) = coeff;
        lap(i,i,j-1) = coeff;
        lap(i,i,j+1) = coeff;
        lap(i,i,j) = 4*coeff;
    end
end
%% * Set boundary conditions
Q0 = 0;
Q = Q0*ones(N);
fprintf('Source at (%1d,%1d) equal to %1d \n',xs,ys,max(max(S)));
fprintf('Flux is zero on all boundaries\n');


%% * Compute the Crank-Nicolson matrix

			 
%% * Initialize the wavefunction 


%% * Plot the initial wavefunction



%% * Initialize loop and plot variables 


%% * Loop over desired number of steps (wave circles system once)


