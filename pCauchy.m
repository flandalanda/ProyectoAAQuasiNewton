function [pC] = pCauchy(B,g, delta)
% In: B     ... (symmetric matrix) approximate of the hessian of f at xk
%     g     ... (vector) gradient of f in xk
%     delta ... trust region radius
%
% Out: pC    ... The Cauchy point
% 
% Closed form of the Cauchy point taken from Numerical Optimization by
% Jorge Nocedal and Stephen J. Wright.


% We calculate the constants we will need
normg = norm(g);
termCuad = g'*B*g;

% We determine the largest step we can take in the direction of steepest
% descent and remain within the trust region
p = -delta/normg*g;

% We preemptively set our constant to the largest step possible
tau = 1;

% If g'*B*g is positive we need to verify our step does not leave the trust
% region
if termCuad > 0
    tau = min(normg^3/(delta*termCuad),1);
end 


pC = 0.99*tau * p;

end

