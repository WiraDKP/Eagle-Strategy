clear;
clf;

%% Levy Flight Visualization
% a = [ 0 0 ];
% hold on;
% for i = 1 : 1000
% a = [ a ; a(end,:) + levy(1.9,[ 1 2 ])];
% plot(a(:,1),a(:,2),'-k'); pause(0.001);
% end

problem = @(x) problemfunc(x);
nvar = 2;

bound.xmin = -10;
bound.xmax = 10;

param.itermax = 50;
param.npop = 50;
param.gamma = 1;
param.beta0 = 1;
param.alpha = 0.2;
param.damp = 0.9;
param.scale = (bound.xmax-bound.xmin);
param.lambda = 1.5;

ES(problem, nvar, bound, param)

% [X,Y] = meshgrid(bound.xmin:0.01:bound.xmax, bound.xmin:0.01:bound.xmax);
% surf(X,Y,3*(1-X).^2.*exp(-(X.^2) - (Y+1).^2) - 10*(X/5 - X.^3 - Y.^5).*exp(-X.^2-Y.^2) ... 
%                 - 1/3*exp(-(X+1).^2 - Y.^2), 'EdgeColor','none');
% view(40,55); saveas(gcf, 'Peak.png');