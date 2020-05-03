function [x,iter] = rcSR1(f, x0, itmax)
% In :  f     ... (handle) function to be optimized
%       x0    ... (vector) initial point
%       itmax ... (natural number) upper bound for number of iterations
%
% Out:  x    ... (vector) last approximation of a stationary point
%       iter ... (natural number) number of iterations

% We set global parameters
tol = 10e-5;
r = 10e-6;
eta = 0.1;
deltaMax = 1.25;

% Initializing initial aproximations and iteration count
iter = 0;
x = x0;
delta = deltaMax;
grad = apGrad(f,x);
B = apHess(f,x);
H = inv(B);


% The loop stops when the maximum number of iterations is reached or when
% we are close enough to a stationary point
while(norm(grad)>tol && iter <itmax)
    
    % We calculate the newton direction
    newt = -H*grad;
    
    % In case the Newton direction is a descent direction (which it might 
    % not be due to B not necessarily being spd) we take it as an
    % approximate solution to the trust region problem
    if dot(newt, grad) < 0
        s = newt;
        tam = norm(s);
        % If the Newton direction fails the boundary condition, we resize
        % the vector to the largest possible size
        if tam > delta
            s = delta * s / tam;
        end
    % If the Newton direction is not a descent direction we take the Cauchy
    % point as the solution to the trust region problem
    else
        s = pCauchy(B,grad,delta);
    end
    
    % Actual reduction is calculated
    x1 = x + s;
    df = f(x)-f(x1);
    % Predicted reduction is calculated as mk(0)-mk(s)
    dm = -grad'*s - 0.5*s'*B*s;
    % We calculate the reduction quotient
    rho = df/dm;
    
    % If the model fits the function well enough, we take the step
    if rho > eta
       x = x1;
       
       % We update auxiliary variables
       grad1 = apGrad(f,x);
       gamma = grad1 - grad;
       grad = grad1;
       aux1 = gamma - B*s;
       aux2 = s - H * gamma;
       
       % If the necessary conditions are met, we calculate the new
       % approximations for B and H
       if abs(dot(aux1, s)) >= r*norm(s)*norm(aux1)
           B = B + aux1*aux1'/dot(aux1, s);
           H = H + aux2*aux2'/dot(aux2, gamma);
       end
       
       % If the model fits the function exceptionally well and the step
       % taken was relatively large, we extend the trust region
       if rho > 0.75 && norm(s) > 0.8 * delta
           delta = min(2 * delta, deltaMax);
       end
       % We increment the number of iterations performed
       iter = iter +1;
       
    % If the model was a poor fit, we reduce the trust region
    else
        delta = 0.5*delta;
    end
    
end

end

