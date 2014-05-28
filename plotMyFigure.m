function plotMyFigure(L, N, Q, Qplot, Qresh, S, TN,  ...
    TNxN, coeff, cons, dM,  ...
    deltaFunction, fin, h, i, iter, iter2,  ...
    maxZ, n1, sparsedM, stepNumber, tau,  ...
    tfinal, time, tsteps, w, x, xExponent,  ...
    xL, xs, y, yExponent, yL, ys)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% load 'matlab.mat';

figure('Units', 'pixels', ...
    'Position', [100 100 600 600]);clf
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

end
