% Code to demonstrate N = 13 is infeasible for 3d kissing problem. 
% Utilize Putinar's p-satz to find lower bound ... 
% Written in SPOTLESS.

function [gama_list, status_list, num_vars, nineq, n0] = function_dsos_putinar_pursuit_ver2(N, on_PC,...
    iter, max_deq, max_dineq, max_d0, sc, method, tol)

    folders = genpath('./spotless-spotless_isos'); % load spotless
    addpath(folders);

    if (on_PC)
        addpath 'C:\Users\mosek\8\toolbox\r2014aom' % load mosek  
    else
        addpath '/home/tdn/kissing/mosek_for_linux/8/toolbox/r2014aom' % load mosek  
    end

    x = msspoly('x',3*N); % Array form
    z = reshape(x,[3 N]); % Matrix form 
    y = msspoly('y',1);

    % Degrees of multipliers
    deq = 0:max_deq; % Degree for multipliers on equality constraints
    dineq = 0:max_dineq; % Degree for multipliers on inequality constraints
    m_dineq = 0:(max_dineq/2); % Degree of sufficient monomial vector 
    d0 = 0:max_d0; % Degree for s0 term
    m_d0 = 0:(max_d0/2); % Degree of sufficient monomial vector 

    dir_name = strcat('function_dsos_putinar_pursuit_ver2_N=',num2str(N),'deq=',num2str(max_deq),...
        'dineq=',num2str(max_dineq),'d0=',num2str(max_d0),'method=',num2str(method),'tol=',num2str(tol));

    % make folder if one does not exist already
    if exist(dir_name) ~= 7
        mkdir(dir_name)
    end

    gama_list = zeros(1,1);
    status_list = cell(1,1);
    min_max_Q0_eig_list = cell(1, 1);
    min_max_Q_eig_list = cell(1, N, N);
    num_vars = zeros(1, 1); 
    num_constrs = zeros(1, 1); 

    % Initialize first iteration 
    nineq = length(monomials([z(:,1);z(:,2);y],m_dineq));
    n0 = length(monomials([x;y],m_d0)); 

    % Initialize if meant to be first iteration 
    if (iter == 1)
        % One initialization method: everything is identity matrix 
        if method == 1
            Uineq = cell(N,N);
            for i = 1:N
                for j = (i+1):N % since we want i < j
                    Uineq{i,j} = eye(nineq);
                end
            end
            U0 = eye(n0);
        elseif method == 2
            rng(method, 'twister') % for reproducibility 
            Uineq = cell(N,N);
            for i = 1:N
                for j = (i+1):N % since we want i < j
                    Uineq{i,j} = eye(nineq) - sc*diag(rand(nineq,1));
                end
            end
            U0 = eye(n0) - sc*diag(rand(n0,1));
        elseif method == 3
            rng(method,'twister') % for reproducibility 
            Uineq = cell(N,N);
            for i = 1:N
                for j = (i+1):N % since we want i < j
                    Uineq{i,j} = eye(nineq) - sc*diag(randn(nineq,1));
                end
            end
            U0 = eye(n0) - sc*diag(randn(n0,1));
        end
    % load U0 and Unieq of previous iteration 
    else
        file_name = strcat(dir_name,'/iter=',num2str(iter-1),'.mat') ;
        load(file_name)
        clear gama
    end

    % Declare program
    prog = spotsosprog;
    prog = prog.withIndeterminate(x);
    prog = prog.withIndeterminate(y);

    % Declare gama (objective)
    [prog,gama]=prog.newFree(1);

    expr = 0; % This will be built incrementally

    % Create multipliers for equality constraints and add terms to expr
    L = cell(N,1); % Cell array for storing multipliers
    for k = 1:N
        % Create multiplier
        v = monomials(z(:,k),deq); % Monomials for multiplier 
        [prog,L{k}] = prog.newFreePoly(v,1); % 1 new free polynomial (multiplier) 

        % Add to expr
        expr = expr + L{k}*(z(1,k)^2 + z(2,k)^2 + z(3,k)^2 - 4);
    end

    s = cell(N,N); % Cell array for storing multipliers (some of these will be empty)
    Q = cell(N,N); % cell array for storing inner most Gram matrix 
    for i = 1:N
        for j = (i+1):N % since we want i < j
            v = monomials([z(:,i);z(:,j);y],m_dineq);
            % Want s{i,j} = v'Uineq'*Q*Uineq*v for some dd Q
            [prog,Q{i,j}] = prog.newDD(nineq);
            new_v = Uineq{i,j}*v;
            s{i,j} = new_v'*Q{i,j}*new_v;

            expr = expr + s{i,j}*( (z(1,i)-z(1,j))^2 + (z(2,i)-z(2,j))^2 + ...
                (z(3,i)-z(3,j))^2  + y - 4 );
        end
    end

    v = monomials([x;y],m_d0); 
    [prog,Q0] = prog.newDD(n0);
    new_v = U0*v;
    s0 = new_v'*Q0*new_v;
    
    expr = expr + s0;
    expr = gama - y + expr;

    prog = prog.withPolyEqs(expr); % y - gama = ... 

    % Run (s)dsos program
    options = spot_sdp_default_options();
    options.solveroptions.MSK_IPAR_BI_CLEAN_OPTIMIZER = 'MSK_OPTIMIZER_INTPNT'; % Use just the interior point algorithm to clean up
    options.solveroptions.MSK_IPAR_INTPNT_BASIS = 'MSK_BI_NEVER'; % Don't use basis identification (it's slow)
    options.verbose = 0;

    disp('Done setting up constraints');

    sol = prog.minimize(-gama, @spot_mosek, options);

    status = sol.status; 
    try
        gama = double(sol.eval(gama));
    catch
        disp('Cannot evaluate current iteration gama because problem is primal infeasible')
        % save progress thus far
        file_name = strcat(dir_name,'/','eigen_info.mat');
        save(file_name, 'min_max_Q0_eig_list','min_max_Q_eig_list')

        disp(gama_list)
        disp(status_list)
        return 
    end

    % Update bases 
    Q0_val = double(sol.eval(Q0));

    [V,D] = eig(Q0_val);

    min_max_Q0_eig = {min(diag(D)), max(diag(D))};

    % change the basis only if eigenvalues are larger than threshold value 
    % not the reason the seg fault is happening 
    if 0
        new_U0 = U0;
        if max(diag(D)) >= tol
            new_U0 = sqrt(D)*V'*U0;
        end
        change_U0 = norm(new_U0-U0,'fro');
        U0 = new_U0;
    else
        new_U0 = sqrt(D)*V'*U0; 
        change_U0 = norm(new_U0-U0,'fro');
        U0 = new_U0;
    end

    Qineq_val = cell(N,N);
    change_Uineq = cell(N,N);
    for i = 1:N
        for j = (i+1):N % since we want i < j
            Qineq_val{i,j} = double(sol.eval(Q{i,j}));
            [V,D] = eig(Qineq_val{i,j});
            max_eig = max(diag(D));
            min_eig = min(diag(D));
            min_max_Q_eig_list{i,j} = {min_eig, max_eig};

            U = Uineq{i,j};
            if 0
                new_Uineq = U;
                if max_eig >= tol
                    new_Uineq = sqrt(D)*V'*U;
                end
                change_Uineq{i,j} = norm(new_Uineq-U,'fro');
                Uineq{i,j} = new_Uineq;
            else
                new_Uineq = sqrt(D)*V'*Uineq{i,j};
                change_Uineq{i,j} = norm(new_Uineq-Uineq{i,j},'fro');
                Uineq{i,j} = new_Uineq;
            end
            % try to explicitly delete variables from the workspace 
        end
    end

    disp('Current iteration is');
    disp(iter);
    disp('Current iteration success in solving LP');
    disp(status);
    disp('Current iteration optimal value of gama');
    disp(gama);
    disp('Current iteration Q0');
    disp(Q0_val);
    disp('Distance between new and old U0 basis matrices');
    disp(change_U0);
    disp('Minimum and maximum eigenvalue of Q0 i.e. gram matrix');
    disp(min_max_Q0_eig);
    disp('Distance between new and old Uineq basis matrices');
    disp(change_Uineq);
    % disp('Minimum and maximum eigenvalue of Q{i,j} i.e. gram matrix');
    % disp(min_max_Q_eig_list{iter,:,:});

    % save U0 and Uineq for next iteration 
    file_name = strcat(dir_name,'/iter=',num2str(iter));
    save(file_name, 'gama','Qineq_val','Q0_val', 'U0', 'Uineq')
end
