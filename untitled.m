%%
% heat_equation - Program to solve the diffusion equation
% using the Backward Euler method
clear; help heat_equation; % Clear memory and print header
close all force;clc;
%* Initialize parameters (time step, grid spacing, etc.)
tau = 1e-3; % Enter time step
N = 100; % Number of grid points
L = 2; % The system extends from (x)=(0) to (x)=(L)
h = L/(N-1);
w = 0.2;
coeff = tau/h^2;
tfinal = 2;
tau = .5*h;
tsteps = ceil(tfinal/tau);
%% Initialize Source function
[y,x]=meshgrid(-1:h:1,-1:h:1);
xp=x;
yp=y;
S = exp((-xp.^2-yp.^2)/.2);
surf(x,y,S)
%% * Set up the Laplacian operator matrix
lapx = zeros(N);  % Set all elements to zero
for i=2:(N-1)
    lapx(i,i-1) = 1;
    lapx(i,i) = -2;  % Set interior rows
    lapx(i,i+1) = 1;
end
% Boundary conditions
lapx(1,1)=-1;
lapx(1,2)=1;
lapx(N,N)=-1;
lapx(N,N-1)=1;

%% * Initialize Q-matrix
Q = S; 
%% * Compute A-matrix (Tn+1)=ATn
dM = (eye(N) - coeff*(lapx));
Qplot(:,:,1) = Q(:,:);
for i = 2:ceil(tfinal/tau)
    Qn = dM\Q(:,:);
    Q(:,:) = Qn;
    Qplot(:,:,i)=Q(:,:);
    figure(1);
    surf(x,y,Qplot(:,:,i));
    title(sprintf('time=%g',tau*i));   
end
%%
disp 'done'