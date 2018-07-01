%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Original contributor: Tin Nguyen. 
% Description: MATLAB function,
% - initialize the polynomial component (Q) of the Lyapunov function, 
% - theoretical basis: Lyapunov theory for linear dynamical system.
% Dependencies: 
% - YALMIP (open-sourced LMI parser that comes with good documentation) 
% - optimizers: MOSEK, SDPT3, etc. 
% Inputs: 
% - n: dimensionality, 
% - x: array of YALMIP variables e.g. [x1, x2], 
% - f: array of YALMIP expressions e.g. [-x1 + x1*x2, -x2], 
% - deg: the proposed degree of the polynomial Q, 
% - scale: multiply each coefficient of the Lyapunov function by this value, 
% - noise: add i.i.d. Gaussian r.v. with this variance to the scaled Lyapunov function.
% - options: YALMIP options, such as printing progress of optimizers etc. 
% Outputs: 
% - success: binary variable indicating success of the search,
% - guess: coefficients array of Q, if the search is a success. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [success, guess] = intQ_Qh(n, x, f, deg, noise, scale, options)    
    
    indices = 1:n;
    singletons = sdisplay(x);
    dict = containers.Map(singletons,indices);

    % linearize f near the origin: build A, the transformation matrix, 
    % of this linearization 
    assign(x,zeros(n,1));
    A = zeros(n);
    linearized_f = linearize(f);
    for i = 1:n
        [c,v] = coefficients(linearized_f(i));
        len = length(v);
        for j = 1:len
            if (isa(v(j),'double'))
                continue
            else
                single = sdisplay(v(j));
                index = dict(single{1});
                A(i,index) = c(j);
            end
        end
    end
    
    % find quadratic Lyapunov function
    V = sdpvar(n); 
    constraints = [V >= 0, A'*V+V*A <= 0, trace(V) == 1];
    optimize(constraints,[],options);
    
    % return coefficients of quadratic Lyapunov function
    try
        logLya = x'*value(V)*x;
    catch
        success = 0
        guess = 0
        return 
    end
    logLya = (logLya)^(deg/2);
    [sparse_coef, sparse_basis] = coefficients(logLya);
    sparse_len = length(sparse_coef);
    
    full_basis = monolist(x,deg,deg);
    full_len = length(full_basis);
    full_coef = zeros(full_len,1); 
     
    for i = 1:sparse_len
        v = sdisplay(sparse_basis(i));
        % index = dict2(v{1});
        for j = 1:full_len
            if (strcmp(v,sdisplay(full_basis(j))))
                    full_coef(j) = sparse_coef(i);
                    break
            end
        end
    end
   
    success = 1;
    guess = scale*full_coef + noise*randn(size(full_coef));
end
