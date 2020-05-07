function [alpha, gnew] = lineSearch( f, xk, dk, gk )
% In :  f  ... objectve function (handle)
%       xk ... current point
%       dk ... chosen direction of descent
%       gf ... gradient of f in xk
%
% Out:  alpha ... a parameter satisfying (W1) and (W2)
%       gnew  ... gradient of f in  xk+alpha*dk

% We define global constants
c1 = 10e-4;
c2 = 0.99;
phiPrime0 = dk'*gk;
phi0 = f(xk);
maxIter = 1000;

% We initialize variables necessary throught the itertation process
iter = 1;
aMax = 4;
a0 = 0;

% We limit our iterations to a fixed number
while(iter < maxIter)
    % We select an alpha in the proposed interval
    alpha = 0.5*(a0+aMax);
    % We evaluate the function at the proposed point
    phiA = f(xk + alpha*dk);
    % We approximate the gradient at the new proposed point and calculate
    % the derivative with respect to alpha at said point
    gnew = apGrad(f, xk + alpha*dk);
    phiPrimeA = gnew'*dk;
    
    % We check if any conditions have been violated
    if phiA > phi0 + c1*a0*phiPrime0 || (phiA >= phi0 && iter >1)
        [alpha, gnew] = zoom( a0, alpha, f, xk,dk, phi0, phiPrime0);
        break
    elseif abs(phiPrimeA) <= -c2*phiPrime0
        break
    elseif phiPrimeA >= 0
        [alpha, gnew] = zoom(alpha, a0, f, xk, dk, phi0, phiPrime0);
        break;
    else
        % We refine the interval
        a0 = alpha;
        phi0 = phiA;
        iter = iter + 1;
    end
end

end

% Definition of auxiliary function zoom
function [alpha, gnew] = zoom(alo, ahi, f, xk, dk, phi0, phiPrime0)
    c1 = 10e-4;
    c2 = 0.99;
    
    iter = 1;
    maxIter = 1000;
    phiLo = f(xk + alo*dk);
    
    while(iter < maxIter)
       % We select an alpha within the know interval and calculate the
       % value of the function at the new point, as well as the gradient
       % and derivative
       alpha = 0.5*(alo+ahi);
       phiA = f(xk+ alpha*dk);
       gnew = apGrad(f, xk + alpha*dk);
       phiPrime1 = gnew'*dk;
       
       % We check if any condition is not met
       if((phiA > phi0 + alpha*c1*phiPrime0) || (phiA >=phiLo))
           % If the step was too long, we refine the interval 
           ahi = alpha;
           % If W2 is met, we found a correct alpha
       elseif abs(phiPrime1) <= -c2*phiPrime0
           break
       else
           % We refine the interval depending on the conditions met
           if phiPrime1 * (ahi-alo)>=0
               ahi = alo;
           end
           alo = alpha;
           iter = iter+1;
       end
    end
end

