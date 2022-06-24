function [A,B,C,D] = iodmd(U,X,Y)%,dt_rec,dt_sys)
%% 
    thres = 1e-12;
    m = size(U,1);
    n = size(X,1);

    X_in = X(:,1:end-1);
    X_out = X(:,2:end);

    [u,d,v] = svd([X_in;U],'econ');

    r_tilde = length(find(diag(d)>thres));

    u = u(:,1:r_tilde);
    v = v(:,1:r_tilde);
    di = diag(1.0./diag(d(1:r_tilde,1:r_tilde)));

%% 

%     thres2 = 1e-13;
%     [uhat,dhat,~] = svd(X_out,'econ');
%     r = length(find(diag(dhat)>thres2));
%     uhat = uhat(:,1:r);

%% Matrices

%     D = D.^(floor(dt_rec/dt_sys));
%     A = Phi*D*pinv(Phi);
%     B = uhat*uhat'*X_out*v*di*u(n+1:n+m,:)';
%     C = Y(:,1:end-1) * v * di * u(1:n,:)';
%     D = Y(:,1:end-1) * v * di * u(n+1:n+m,:)';

    A = X_out * v * di * u(1:n,:)';
    
%     A_tilde = u(1:n,:)'* X_out * v * di;
%     [W,Lmb] = eig(A_tilde);
%     Phi = inv(Lmb)*X_out*v*di*W;

    B = X_out * v * di * u(n+1:n+m,:)';
    C = Y(:,1:end-1) * v * di * u(1:n,:)';
    D = Y(:,1:end-1) * v * di * u(n+1:n+m,:)';

end