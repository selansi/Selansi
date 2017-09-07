function cx=IF_FeedbackMechanism(x,dato,ci)
%%%
% cx=IF_FeedbackMechanism(x,dato,ci)
% Defining a Hill type function representing the feedback mechanism 
%%%

H=dato.H;
K=dato.K;
epsilon=dato.epsilon;
n=dato.n_gene;
IF_Type=dato.IF_Type;
switch IF_Type
    case 'Hill'
        indreg=cell(n,1);
        for i=1:n
            indreg{i}=find(H(i,:));
        end
        RHO=cell(n,n);
        for i=1:n
            for j=1:n
                if H(i,j)>0
                    RHO{j,i}= x{j}.^H(i,j)./(x{j}.^H(i,j)+K(i,j)^H(i,j));
                elseif H(i,j)==0
                    RHO{j,i}=1;
                else
                    RHO{j,i}=K(i,j)^(-H(i,j))./(x{j}.^(-H(i,j))+K(i,j)^(-H(i,j)));
                end
            end
        end
        ON_RHO=cell(n,n);
        for i=1:n
            for j=1:n
                ON_RHO{i,j}=1-RHO{i,j};
            end
        end
        
        S_OFF_ON=cell(n,2*n);
        S_OFF_ON(:,1:2:2*n-1)=RHO;
        S_OFF_ON(:,2:2:2*n)=ON_RHO;
        
        % Calculate the c_i(x) for i=ci 
        ind_J=length(indreg{ci});
        CX=cell(1,2^ind_J);
        CX{1}=S_OFF_ON{indreg{ci}(1),2*ci-1}; 
        CX{2}=S_OFF_ON{indreg{ci}(1),2*ci}; 
        if ind_J>1
            for i=2:ind_J
                AAA=CX(1:2^(i-1));
                BBB=S_OFF_ON(indreg{ci}(i),2*ci-1:2*ci);
                [ja,jb] = meshgrid(1:2^(i-1),[1 2]);
                for jj=1:2^i
                    CX{jj}=AAA{ja(jj)}.*BBB{jb(jj)};
                end
            end
        end
        cx=zeros(size(x{1}));
        for i=1:2^ind_J
            cx=cx+CX{i}*epsilon{ci}(i);
        end
            
    case 'to define'
        
end
