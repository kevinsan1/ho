%% Heat_equation 
%- Program to solve the diffusion equation
% using the Backward Euler method
%% Parameters
clear all; help heat_equation; % Clear memory and print header
close all;clc;
saveFigurePath = ['/Users/kevin/SkyDrive/KTH Work/' ...
'Period 3 2014/DN2255/Homework/1/Heat Equation/Figures/'];
%* Initialize parameters (time step, grid spacing, etc.)
N = 50; % Number of grid points
L = 1; % The system extends from (x)=(0) to (x)=(L)
h = L/N;
i = 1:(N); % 1:N
[x,y] = meshgrid(h/2:h:L,h/2:h:L);
w = 0.2;
xs = 0.5;
ys = 0.5;
tfinal = .5;
tau = .01*h;
coeff = tau/h^2;
tsteps = ceil(tfinal/tau);
%% Initialize Source function
xExponent = (x-xs).^2;
yExponent = (y-ys).^2;
S = exp(-(xExponent)/w^2).*exp(-yExponent/w^2);
deltaFunction = zeros(N);
deltaFunction(round(N/2),round(N/2))=2;
% S = deltaFunction;
%% * Initialize Q-matrix
S = reshape(S,[N^2,1]);
Q = S; 
%% * Compute A-matrix (Tn+1)=ATn
% dM = eye(N) - coeff*lap; 
TN = 2*eye(N) - diag(ones(N-1,1),1) - diag(ones(N-1,1),-1);
TNxN = kron(eye(N),TN) + kron(TN,eye(N));
dM = eye(N^2) + coeff*TNxN;
sparsedM = sparse(dM);
%% Initialize loop and plot variables
Qplot(:,1) = Q; % initial value
stepNumber=round(.25/tau);
%% Main loops
for iter=1:stepNumber
    Q = sparsedM\Q + S;
    Qplot(:,iter) = Q(:);
end
% Loop after source is gone
for iter2=(iter+1):tsteps
    Q = sparsedM\Q;
    Qplot(:,iter2) = Q(:);
end
Qresh = reshape(Qplot,[N,N,tsteps]);
time = linspace(0,tfinal,tsteps);
%% look at dx*dy*Qij for Conservation
cons(1:tsteps,1) = h^2*sum(sum(Qresh(:,:,1:end)));
plot(tau*(1:length(cons)),cons);
xL = xlabel('time (sec)');
yL = ylabel('\Deltax \Deltay Q_{i,j} ');
set([xL,yL], 'FontSize',12);
%%
% figure(1);clf;
% for i = 1:tsteps
%     surf(x,y,Qresh(:,:,i))
%     title(sprintf('%g',time(i)));
%     axis([0 1 0 1 0 max(max(Qplot))]);
%     hold off;
%     pause(0.02)
% end
%% Print Plots
fin = length(Qplot(1,:));
maxZ = max(max(max(Qresh)))
figure('Units', 'pixels', ...
    'Position', [100 100 600 600]);clf
%     n1 = [1 ceil(.25/(2*tau)) ceil(.25/tau)...
%        ceil(.25/tau)+5 ceil(.25/tau)+20 ceil(.25/tau)+50];
n1 = [1 floor(fin*6/24) floor(fin*12/24)...
    floor(fin*13/24) floor(fin*14/24) floor(fin*15/24)];
%     n1 = [1 floor(fin*1/2) floor(fin*3/4) fin];
for i = 1:length(n1)
    s1=subplot(3,2,i);
    mesh(x,y,Qresh(:,:,n1(i)))
    hTitle = title(sprintf('t = %0.4f',time(n1(i))));
    axis([0 1 0 1 0 maxZ]);
	hXLabel = xlabel('x');
	hYLabel = ylabel('y');
	hZLabel = zlabel('Q');
	set( gca                       , ...
	    'FontName'   , 'Helvetica' );
	set([hTitle, hXLabel, hYLabel, hZLabel], ...
	    'FontName'   , 'AvantGarde');
	set( gca             , ...
	    'FontSize'   , 8           );
	set([hXLabel, hYLabel, hZLabel]  , ...
	    'FontSize'   , 10          );
	set( hTitle                    , ...
	    'FontSize'   , 12          , ...
	    'FontWeight' , 'bold'      ); 
	set(gca, ...
	    'Box'         , 'off'         , ...
	    'TickDir'     , 'out'         , ...
	    'TickLength'  , [.02 .02]     , ...
	    'XMinorTick'  , 'on'          , ...
	    'YMinorTick'  , 'on'          , ...
	    'XColor'      , [.3 .3 .3]    , ...
	    'YColor'      , [.3 .3 .3]    , ...
	    'LineWidth'   , 1             );
end
% 	hold off;
%% Plot 1
printYesNo = 0;
if printYesNo == 1
set(figure(1), 'PaperPositionMode', 'auto');
print('-depsc2', [saveFigurePath ...
    sprintf('exponentialFunctionPlot')]);
end