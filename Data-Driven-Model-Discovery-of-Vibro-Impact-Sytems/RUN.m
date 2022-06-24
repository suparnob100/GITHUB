close all
clear all
clear global
clc;

%% Define variables

global rho A E I EDOF cv cm;
cm = 2;
cv = 1;
E=210e9; % Young Modulus
rho=2700; % Density
A=0.05*0.05; % Cross section area
I=(0.05*0.05^3)/12; % Second moment of inertia
Le=1; % Length
Node_DOF=2;
NOE=31;

GDOF=Node_DOF*(NOE+1);
Nodes_per_Elem = 2;
EDOF=Nodes_per_Elem*Node_DOF;


%% Generate Mesh - 1D

nodeco=linspace(0,Le,NOE+1)';
B = zeros(NOE,EDOF);

for i=1:NOE
    B(i,1:Node_DOF)=i*Node_DOF-[(Node_DOF-1):-1:0];
    B(i,Node_DOF+1:EDOF) =(i+1)*Node_DOF-[(Node_DOF-1):-1:0];
end

% nodes on boundary
bDOF=[1,2]';


%% Determine stiffness and mass matrix

[K,M,C,F,Ida] = assmbl(nodeco,B,bDOF);


%% Forcing function

fs = 5e3;
t = 0:1/fs:100-1/fs;
% f = chirp(t,10,t(end),50);
% frc = @(tp) f(t==tp);
% U = frc(t);
U = (sin(5*t));


%% Full state-space data obtained from simulation

A_sys = full([zeros(size(M)), eye(size(M,1)); -M\K -M\C]);

B_sys = zeros(size(A_sys,1),1);
B_sys(end/2+5) = 1;
B_sys_split = B_sys(size(M,1)+1:end);
B_sys(size(M,1)+1:end) = M\B_sys_split;

C_sys2 = eye(size(A_sys,1));
D_sys2 = zeros(size(C_sys2,1), 1);

full_sys=ss(A_sys,B_sys,C_sys2,D_sys2);
x=lsim(full_sys,U,t);


%% Partial observation data

obsDOF = [5 13 23 35 47 59]-2;
C_sys = zeros(length(obsDOF),size(A_sys,1));
C_sys(sub2ind(size(C_sys), 1:length(obsDOF), obsDOF)) = 1;
D_sys = zeros(size(C_sys,1), 1);

obs_sys=ss(A_sys,B_sys,C_sys,D_sys);
y=lsim(obs_sys,U,t);


%% Reconstructed system using iodmd

[sysA,sysB,sysC,sysD] = iodmd(U(1:end-1),x',y');

z = eigplot(sysA);

if z~=0
    
    simulated_sys = ss(sysA,sysB,sysC,sysD,1/fs);
    y_rec=lsim(simulated_sys,U,t);
    
    %% Compare reconstructed data with originally observed data
    
    subplot(2,1,1)
    plot(t,y(:,end),t,y_rec(:,end))
    error_s = norm(y(:,end)-y_rec(:,end))/norm(y(:,end));
    legend('True','Reconstructed','interpreter','latex');
    ylabel(strcat('$w(1,t)$'),'Interpreter','Latex');
    title(strcat('Relative error $\widehat{e} = $', num2str(error_s)),'Interpreter','Latex');
    
    subplot(2,1,2)
    semilogy(t,abs(y(:,end)-y_rec(:,end))./abs(y(:,end)));
    ylabel('Relative error','Interpreter','Latex');
    xlabel('time ($t$)','Interpreter','Latex');
    set(gca,'FontSize',15);
    
    %% Wavelet-based analysis
    
    Z=[];
    
    for i = 1:size(y,2)
        
        m = modwt(y(:,i)');
        mra = modwtmra(m);
        Z = [Z;mra];
        
    end
    
    [sysA2,sysB2,sysC2,sysD2] = iodmd(U(1:end-1),Z,y');
    figure(3)
    stb = eigplot(sysA2);
    simulated_sys_wt = ss(sysA2,sysB2,sysC2,sysD2,1/fs);
    y_rec_wt=lsim(simulated_sys_wt,U,t);
    
    if(stb==0)
        %% Compare reconstructed data with originally observed data
        
        figure(2)

        subplot(2,1,1)
        plot(t,y(:,end),t,y_rec_wt(:,end))
        error_s = norm(y(:,end)-y_rec_wt(:,end))/norm(y(:,end));
        legend('True','Reconstructed','interpreter','latex');
        ylabel(strcat('$w(1,t)$'),'Interpreter','Latex');
        title(strcat('Relative error $\widehat{e} = $', num2str(error_s)),'Interpreter','Latex');

        subplot(2,1,2)
        semilogy(t,abs(y(:,end)-y_rec_wt(:,end))./abs(y(:,end)));
        ylabel('Relative error','Interpreter','Latex');
        xlabel('time ($t$)','Interpreter','Latex');
        set(gca,'FontSize',15);
        
    end
    
end