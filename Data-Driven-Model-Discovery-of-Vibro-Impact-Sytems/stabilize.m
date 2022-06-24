function [A,B,C,D,data] = stabilize(A,B,C,D,U,X,Y)

    maxit = 1000;
    H = [A,B;C,D];
    prob = makeProblem(U(:,1:end-1),X(:,1:end-1),X(:,2:end),Y(:,1:end-1));
    [halt_log_fn, get_log_fn] = makeHaltLogFunctions(maxit);
    opts = gransoOptions(prob.nvar);
    opts.x0 = H(:);
    opts.maxit = maxit;
    opts.print_level = 0;
    opts.halt_log_fn = halt_log_fn;
    opts.opt_tol = 1e-8;
%     opts.limited_mem_size = 10; % Uncomment to use Limited-Memory-BFGS (not recommended)
    [soln,t] = runGranso(prob,opts);
    soln = rmfield(soln,'H_final'); % This is critical
    data.soln = soln;
    data.log = get_log_fn();
    data.log.t = t;
    [A,B,C,D] = prob.xToSubMatrices(soln.most_feasible.x);
    
end