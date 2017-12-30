function [bestff, mincost] = ES(problem, nvar, bound, param)

    func = problem;
    varsize = [1 nvar]; 
    
    itermax = param.itermax;
    npop = param.npop;
    gamma = param.gamma;
    b0 = param.beta0;
    alpha = param.alpha;
    damp = param.damp;
    scale = param.scale;
    lambda = param.lambda;
    
    if (lambda <= 1 | lambda > 3)
        error("=========================================" + "\n" + ...
              "Parameter for Levy Flight. 1 < alpha <= 3" + "\n" + ...
              "=========================================", 0);
    end
    
    xmin = bound.xmin;
    xmax = bound.xmax;
    
    globest.cost = inf;
    
    init.loc = [];
    init.cost = [];

    ff = repmat(init, npop, 1);

    % Initiate fireflies
    for i = 1:npop
        ff(i).loc = unifrnd(xmin, xmax, varsize);
        ff(i).cost = func(ff(i).loc);

        if ff(i).cost < globest.cost
            globest = ff(i);
        end
    end

    bestcost = zeros(itermax,1);
    
    tempX = zeros(itermax,npop);
    tempY = zeros(itermax,npop);

    [X,Y] = meshgrid(xmin:0.1:xmax, xmin:0.1:xmax);
    Z = -20*exp(-0.2*sqrt(0.5*(X.^2+Y.^2)))-exp(0.5*(cos(2*pi*X)+cos(2*pi*Y)))+exp(1)+20;
    for iter = 1:itermax
        contour(X,Y,Z, 25); hold on;
        scatter(0,0,35,'ok','filled');
        title('Ackley Function');
        
        % Levy Flight
        for i = 1:npop
            newff(i).loc = min(max(ff(i).loc + levy(lambda, varsize),xmin),xmax);
            newff(i).cost = func(newff(i).loc);
            if (newff(i).cost < ff(i).cost)
                ff(i).loc = newff(i).loc;
                ff(i).cost = newff(i).cost;
            end
            tempX(iter,i) = ff(i).loc(1);
            tempY(iter,i) = ff(i).loc(2);            
        end
        scatter(tempX(iter,:),tempY(iter,:),'or','filled');
        if (iter > 1)
            for i = 1:npop
                line([tempX(iter,i);tempX(iter-1,i)],[tempY(iter,i);tempY(iter-1,i)],'Color','b', 'LineStyle', '-');
            end        
        end
        hold off;
        frame(2*iter-1) = getframe(gcf); pause(0.0001);
        
        contour(X,Y,Z, 25); hold on;
        scatter(0,0,35,'ok','filled');
        title('Ackley Function');            
        for i = 1:npop
            newff(i).cost = inf;
            for j = 1:npop
                if ff(j).cost < ff(i).cost
                    distance = norm (ff(i).loc - ff(j).loc);                 
                    new.loc = min(max(ff(i).loc + b0*exp(-gamma*distance.^2).*(ff(j).loc-ff(i).loc) + alpha*scale*unifrnd(-1,1,varsize),xmin),xmax);
                    new.cost = func(new.loc);
                    if new.cost < newff(i).cost
                        newff(i).loc = new.loc;
                        newff(i).cost = new.cost;                        
                        if newff(i).cost < globest.cost
                            globest = newff(i);
                        end
                    end
                end
            end          
        end
        
        % Choose first npop of best fireflies
        ff=[ff ; newff'];
        [~, order]=sort([ff.cost]);
        ff=ff(order);
        ff=ff(1:npop);
        
        for i = 1:npop
            tempX(iter,i) = ff(i).loc(1);
            tempY(iter,i) = ff(i).loc(2);  
        end
        scatter(tempX(iter,:),tempY(iter,:),'or','filled');
        frame(2*iter) = getframe(gcf); pause(0.0001);
        hold off;
        
        bestcost(iter) = globest.cost;


        disp(['Iteration ' num2str(iter) ' | Minimum cost = ' num2str(bestcost(iter))] );
        alpha = alpha.*damp;
    end
    
    name = ['Ackley pop=' num2str(npop)];
    video = VideoWriter(name,'MPEG-4');
    video.FrameRate = 24;
    open(video)
    writeVideo(video,frame);
    close(video)
    
    clf
    for iter = 1:itermax
        semilogy(1:iter,bestcost(1:iter));
        axis([1 itermax 1e-7 1 ]);
        title(['Minimum Value Plot | Population = ' num2str(npop) ' | MinVal = ' num2str(bestcost(iter))] );
        framec(iter) = getframe(gcf); pause(0.0001);        
    end
    
    name = ['Convergence Ackley pop=' num2str(npop)];
    video = VideoWriter(name,'MPEG-4');
    video.FrameRate = 24;
    open(video)
    writeVideo(video,framec);
    close(video)
    
    bestff = globest;
    mincost = bestcost;
end