clear;
clf;

%% Levy Flight Visualization
% a = [ 0 0 ];
% hold on;
% for i = 1 : 600
%     a = [ a ; a(end,:) + levy(1.5,[ 1 2 ])];
%     plot(a(:,1),a(:,2),'-k'); %pause(0.00001);
%     %axis([-100 100 -100 100 ]);
% 	framel(i) = getframe(gcf);
% end
% name = 'Levy Flight';
% video = VideoWriter(name,'MPEG-4');
% video.FrameRate = 60;
% open(video)
% writeVideo(video,framel);
% close(video)
% a = [ ];
% for i = 1 : 10000
%     a = [ a ; levy(1.5,[ 1 1 ])];
% end
% histogram (abs(a)); hold off;

problem = @(x) problemfunc(x);
nvar = 2;

bound.xmin = -3;
bound.xmax = 3;

param.itermax = 50;
param.npop = 4;
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