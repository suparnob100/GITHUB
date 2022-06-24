function [K,M,C,F,Ida]=assmbl(nodeco,B,bDOF)
%forming ID matrix, global to equation number.

    n = sort(unique(reshape(B,[],1)));

    Id=zeros(size(n));

    for i=1:length(Id)
        if ismember(n(i),bDOF)
            Id(i)=0; % eqn no. is zero corresponding to nodes where diri Boundary cond is defined.
        else
            Id(i)=max(Id)+1; % increase eqn no. corresponding to n.
        end
    end

    Ida=[n,Id];% maps each DOF to equation number.
    [K,M,C,F]=ass(nodeco,B,bDOF,Ida,Id);

end


function [K,M,C,F]=ass(nodeco,B,bDOF,Ida,Id)
    %% Assume rayleigh damping
    

    f=0; % distributed load

    %% element stiffness matrix
    
    global rho A E I EDOF cm cv;

    h=nodeco(2)-nodeco(1);
    element_stiffness =[12 -6*h -12 -6*h;-6*h 4*h^2 6*h 2*h^2;-12 6*h 12 6*h;-6*h 2*h^2 6*h 4*h^2]*E*I/h^3; %element stiffness
    element_mass=[156 -22*h 54 13*h;-22*h 4*h^2 -13*h -3*h^2;54 -13*h 156 22*h;13*h -3*h^2 22*h 4*h^2]*h*rho*A/420;
%     element_damping=cm*element_stiffness + cv*element_mass;
    % elment_force=[f*h -f*h^2/12 f*h/2 f*h^2/12]';

    %% Matrix assemfrc(t)bly
    
    K=sparse(max(Id),max(Id));
    M=K;
%     C=K;
    F=zeros(max(Id),1);

    for e=1:size(B,1)

        lm=Ida(B(e,:)',:); 

        IdIen=[[1:EDOF]',lm(:,2)]; % maps element DOF numbers to equation number.    
        IdIenz=IdIen(IdIen(:,2)~=0,:); % find non zero eqn number

        K(IdIenz(:,2),IdIenz(:,2))=K(IdIenz(:,2),IdIenz(:,2))+element_stiffness(IdIenz(:,1),IdIenz(:,1)); %assemble stiffness matrix based on non zero eqn number.
        M(IdIenz(:,2),IdIenz(:,2))=M(IdIenz(:,2),IdIenz(:,2))+element_mass(IdIenz(:,1),IdIenz(:,1));
%         C(IdIenz(:,2),IdIenz(:,2))=C(IdIenz(:,2),IdIenz(:,2))+element_damping(IdIenz(:,1),IdIenz(:,1));

        Fe=element_stiffness*diri(B(e,:),bDOF);%element force    
        F(IdIenz(:,2))=F(IdIenz(:,2))-Fe(IdIenz(:,1));

    end
    
    C = cm*K + cv*M;
end