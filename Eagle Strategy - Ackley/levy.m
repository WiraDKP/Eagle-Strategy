function v = levy(alpha,varsize)
    % Generate Levy Distribution using Mantegna's Algorithm (1994)
    % Note: For Levy Flight -> 1 < alpha <= 3
    xvar = ( gamma(1+alpha) .* sin(pi.*alpha/2) / gamma((1+alpha)/2.*alpha.*2.^((alpha-1)/2)) ).^(2/alpha);
    yvar = 1;
    x = random('Normal',0,xvar,varsize);
    y = random('Normal',0,yvar,varsize);
    v = x./ (abs(y).^(1/alpha));
end