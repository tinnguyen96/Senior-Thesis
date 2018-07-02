%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Original contributor: Tin Nguyen. 
% Description: MATLAB script,
% - experiment with different configuration of arguments to searchLogLya.m to find Lyapunov function
% + for the dynamical system in Subsection 2.5.1 of Senior Thesis. 
% Dependencies: 
% - YALMIP (open-sourced LMI parser that comes with good documentation) 
% - optimizers: SDPT3, etc. 
% - files in project: intQ_Qh.m, searchLogLya.m 
% Inputs: None. 
% Outputs: 
% - writes to the working directory a file named Sys0.txt, each line corresponding to 
% + a combination of arguments and the results of the Lyapunov search. 
% Room for improvement: 
% - the current script relies on SDPT3 solver but some simple changes can be made to utilize
% + other solvers such as MOSEK (as long as they play nice with YALMIP). 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

%% Make sure YALMIP and SDP solver are visible to MATLAB 
folders = genpath('../YALMIP and SDPT3');
addpath(folders);

%% setup YALMIP variables and expressions 
sdpvar x1 x2
n = 2;
x = [x1; x2];
f = [-x1 + x1*x2; -x2]; 
options = sdpsettings('verbose',0,'solver','sdpt3','sdpt3.maxit',500);

%% try out many different combinations of arguments to searchLogLya.m 
dp = 2:2:8;
dq = 2:2:8;
P_homogeneous = [0;1];
Q_homogeneous = 1;
max_iter = 10;
tol = [1e-8; 1e-7; 1e-6];
noise = [1e-8; 1e-7; 1e-6];
scale = linspace(0.5,2,3);

try
    results = readtable('Sys0.txt');
catch
    emp = cell(0,12);
    results = cell2table(emp);
    results.Properties.VariableNames = {'dp','dq','P_homogeneous',...
        'Q_homogeneous','max_iter','tol','noise','scale','success',...
        'error','P','Q'};
end
entry = cell(1,12);
for i_dp = 1:length(dp)
    entry(1,1) = {dp(i_dp)};
    entry(1,5) = {max_iter};
    for i_dq = 1:length(dq)
        entry(1,2) = {dq(i_dq)};
        for i_P = 1:length(P_homogeneous)
            entry(1,3) = {P_homogeneous(i_P)};
            entry(1,4) = {Q_homogeneous};
            for i_tol = 1:length(tol)
                entry(1,6) = {tol(i_tol)};
                for i_noise = 1:length(noise)
                    entry(1,7) = {noise(i_noise)};
                    for i_scale = 1:length(scale)
                        entry(1,8) = {scale(i_scale)};
                        % disp(entry)
                        [success, error, P, Q] = searchLogLya(n, x, f, dp(i_dp), dq(i_dq), P_homogeneous(i_P), ...
                            Q_homogeneous, max_iter, tol(i_tol), noise(i_noise), scale(i_scale), options);
                        entry(1,9) = {success};
                        entry(1,10) = {error};
                        entry(1,11) = {sdisplay(P)};
                        entry(1,12) = {sdisplay(Q)};
                        results = [results; entry];
                        writetable(results,'Sys0.txt');
                    end
                end
            end
        end
    end
end
writetable(results,'Sys0.txt');
