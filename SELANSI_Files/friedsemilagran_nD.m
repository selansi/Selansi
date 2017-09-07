function [x,Xgrid,TT,PX_sol]=friedsemilagran_nD(path_forder_DR,dato,SLdato)
%%%
% [x,Xgrid,TT,PX_sol]=friedsemilagran_nD(path_forder_DR,dato,SLdato)
% Numerical solution of the Friedman equation in general dimension with 
% the semilagrangian method
%%%

PathCurrent_SL = pwd;

% Dimensionless parameters (gamma1 = 0.01 and gamma2 = 4e-4)
b=cell(dato.n_gene,1);
for i=1:dato.n_gene
    b{i}=dato.R_constants(i,2)/dato.R_constants(i,3);
end      

% Spatial discretization
iN=cell(dato.n_gene,1);
x=cell(dato.n_gene,1);
for i=1:dato.n_gene
    iN{i} = SLdato.Prot_mesh(i,3) + 1;
    x{i} = linspace(SLdato.Prot_mesh(i,1),SLdato.Prot_mesh(i,2), iN{i});
    fprintf('\n The discretization step in the dimension %d is: %g \n',i,x{i}(2)-x{i}(1));
end

% Protein "spatial" Mesh
Xgrid=cell(dato.n_gene,1);
[Xgrid{1:dato.n_gene}] = ndgrid(x{1:dato.n_gene});

% Time definition
t0     = SLdato.Time_mesh(1);
tmax   = SLdato.Time_mesh(2);
nt     = SLdato.Time_mesh(3)*SLdato.Time_mesh(4) + 1;
deltat = (tmax-t0)/nt;
tl     = linspace(t0, tmax, nt);
fprintf('\n The time discretization is: %g \n',tl(2)-tl(1));

% Initial conditions (Gaussian density function)
cd(path_forder_DR);
PX0un=IC_Function(x,Xgrid);
cd(PathCurrent_SL);

% Normalized initial condition
auxnor0=PX0un;
for i=1:dato.n_gene
    auxnor0 = trapz(x{i},auxnor0);
end
% initial condition normalized
PX0=PX0un/auxnor0;

% Computation of characteristics curves
xbar=cell(dato.n_gene,1);
xbarlim=cell(dato.n_gene,1);
for i=1:dato.n_gene
    xbar{i}=x{i}*exp(deltat*dato.R_constants(i,4));
    xbarlim{i} = find(xbar{i}>=x{i}(end));
    fprintf('\nThe length of X%g is: %g and there are %g points of Xbar biggest than Xmax \n',i,length(x{i}),length(xbarlim{i}));
end
% Caracteristics grid
Xbargrid=cell(dato.n_gene,1);
[Xbargrid{1:dato.n_gene}] = ndgrid(xbar{1:dato.n_gene});

% Initialization 
PX      = PX0;

% Time Independent functions
e_x=cell(dato.n_gene,1);
e_lx=cell(dato.n_gene,1);
for i=1:dato.n_gene
    e_x{i}=exp(Xgrid{i}/b{i});
    e_lx{i}=exp(-Xgrid{i}/b{i});
    e_x{i}(isinf(e_x{i})==1)=realmax;
    e_lx{i}(isinf(e_lx{i})==1)=realmax;
end

% Input function, c_i(x), construction:
if strcmp(dato.IF_Type,'Hill')==1
    cx=cell(dato.n_gene,1);
    for i=1:dato.n_gene
        cx{i}=IF_FeedbackMechanism(Xgrid,dato,i);
    end
else
    cd(path_forder_DR);
    cx=cell(dato.n_gene,1);
    for i=1:dato.n_gene
        cx{i}=IF_FM_user(Xgrid,i);
    end
    cd(PathCurrent_SL);
end

% Other time independent functions
sumkmcx=0;
for i=1:dato.n_gene
    sumkmcx = sumkmcx + dato.R_constants(i,1)*cx{i};
end
sumprotdeg=0;
for i=1:dato.n_gene
    sumprotdeg = sumprotdeg + dato.R_constants(i,4);
end
expl_den = 1+(sumkmcx - sumprotdeg)*deltat;


% Saving only dato.Time_mesh(4) simulated times 
nt_sol = SLdato.Time_mesh(4) + 1;
PX_sol      = cell(nt_sol,1);
TT          = zeros(nt_sol,1);
PX_sol{1}   = PX;
TT(1)       = tl(1);    
kk_sol = 1;
jjkk=(nt-1)/(nt_sol-1)+1:(nt-1)/(nt_sol-1):nt;


for j = 2:nt
    %tic
    % PX_bar construction using interpolation 
    if dato.n_gene==1
        PX_bar = interp1(Xgrid{1},PX,Xbargrid{1});
        PX_bar(isnan(PX_bar)==1)=0;
    elseif dato.n_gene>1
        PX_bar = interpn(Xgrid{1:dato.n_gene},PX,Xbargrid{1:dato.n_gene});
        PX_bar(isnan(PX_bar)==1)=0;
    end
       
     % Integral term computation by numerical integration
    Lix=0;
    for i=1:dato.n_gene
        Lix = Lix + dato.R_constants(i,1)/b{i}*e_lx{i}.*cumtrapz(x{i},e_x{i}.*cx{i}.*PX,i);
    end
       
    % Explicit method    
    PX = (PX_bar+deltat*Lix)./expl_den; 
    
    % Zero boundary condition
    CFaux=cell(dato.n_gene,1);
    for i=1:dato.n_gene
        CFaux{i}=':';
    end
    for i=1:dato.n_gene
        CF=CFaux;
        CF{i}=iN{i};
        PX(CF{1:dato.n_gene})=zeros(size(PX(CF{1:dato.n_gene})));
    end
    
    % Normalization: int_xmin^xmax(PX)dx=1
    auxnorpx=PX;
    for i=1:dato.n_gene
        auxnorpx = trapz(x{i},auxnorpx);
    end
    PX=PX/auxnorpx;
            
    % Saving the solution fo the current time step
    if j==jjkk(kk_sol)
        fprintf('Time = %f \n',tl(j))
        PX_sol{kk_sol+1} = PX;
        TT(kk_sol+1)       = tl(j);
        kk_sol = kk_sol+1;
    end
    %toc
end

% Plotting the solution at final time
if dato.n_gene>=4
    figure
    hold on
    PX1_aux=PX;
    for i=dato.n_gene:-1:2
        PX1_aux=sum(PX1_aux,i);
    end    
    PX1=PX1_aux/trapz(x{1},PX1_aux);
    plot(x{1},PX1,'k-','LineWidth',1.5)
    xlabel('Protein 1')
    ylabel('Probability')
    hold off
    
elseif dato.n_gene==3
    figure
    hold on
    PX1=sum(sum(PX,3),2)/trapz(x{1},sum(sum(PX,3),2));
    plot(x{1},PX1,'k-','LineWidth',1.5)
    xlabel('Protein 1')
    ylabel('Probability')
    hold off

    figure
    hold on
    PX2=sum(sum(PX,3),1)/trapz(x{2},sum(sum(PX,3),1));
    plot(x{2},PX2,'k-','LineWidth',1.5)
    xlabel('Protein 2')
    ylabel('Probability')
    hold off

    figure
    hold on
    PX123=sum(sum(PX,2),1)/trapz(x{3},sum(sum(PX,2),1));
    PX3=zeros(1,iN{3});
    for i=1:iN{3}
        PX3(i)=PX123(1,1,i);
    end
    plot(x{3},PX3,'k-','LineWidth',1.5)
    xlabel('Protein 3')
    ylabel('Probability')
    hold off
elseif dato.n_gene==2
    figure
    mesh(Xgrid{1},Xgrid{2},PX_sol{end})
    xlabel('Protein_1')
    ylabel('Protein_2')
    zlabel('Probability')

    figure
    hold on
    plot(x{1},sum(PX,2)/trapz(x{1},sum(PX,2)),'r-','LineWidth',1.5)
    xlabel('Protein 1')
    ylabel('Probability')
    hold off

    figure
    hold on
    plot(x{2},sum(PX,1)/trapz(x{2},sum(PX,1)),'r-','LineWidth',1.5)
    xlabel('Protein 2')
    ylabel('Probability')
    hold off
elseif dato.n_gene==1
    if strcmp(dato.IF_Type,'Hill')==1
        % analytical solution
        a      = dato.R_constants(1)/dato.R_constants(4);        % a=k2/gamma1;
        b      = dato.R_constants(2)/dato.R_constants(3);        % b=k1/gamma2;
        H_par  = dato.H;        % H_par  = 0 if there is not feedback
        k_par  = dato.K;
        eps_par= dato.epsilon{1}(1);
        if H_par > 0
            X_a=Xgrid{1};
            PX_un = (k_par^(H_par)+X_a.^(H_par)).^(a*(eps_par-1)/H_par).*X_a.^(a-1).*exp(-X_a/b);
            PX_a_9    = PX_un/trapz(X_a,PX_un);
        elseif H_par == 0
            X_a=Xgrid{1};
            PX_a_9 = gampdf(X_a,a,b)'; % gamma distribution
        else
            if a*eps_par < 1
                X_a=Xgrid{1};
                X_a(1)=0.1;
            else
                X_a=Xgrid{1};
            end
            PX_un = (k_par^(-H_par)+X_a.^(-H_par)).^(a*(eps_par-1)/H_par).*X_a.^(a*eps_par-1).*exp(-X_a/b);
            PX_a_9    = PX_un/trapz(X_a,PX_un);
        end
    
        figure
        hold on
        plot(Xgrid{1},PX)
        plot(X_a,PX_a_9,'r')
        xlabel('Protein')
        ylabel('Probability')
        legend('Numerical at t_f', 'Analytical')
        hold off
    else
         figure
        hold on
        plot(Xgrid{1},PX,'k-','LineWidth',1.5)
        xlabel('Protein')
        ylabel('Probability')
        hold off
    end
end

end



