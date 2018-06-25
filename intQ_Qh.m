% n is dimensionlity, x is vector of sdpvar and 
% f is description of dynamical system 

% Initialize Q using linearization of the system, perturbed by noise and
% scale 
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
    % disp(A)
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
