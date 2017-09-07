function SELANSI_Plot(name)
%%%
% SELANSI_Plot(name)
% Function to plot the solution obtained with the semilagrangian method
%%%

% Obtaining the actual path
PathCurrent = pwd;

% Loading the solution
load(fullfile(PathCurrent,'DATA',name,'Results','Solution'))

% Plotting the solution
n_gene=length(solution.x);

fprintf('\n The solution can be ploted at different times: \n initial time, T0=%g, choose 0 \n final time, Tend=%g, choose %g \n intermediate time, T= T0 + %g*nt, choose an integer nt such that 1 < nt < %g \n',solution.T(1),solution.T(end),length(solution.T)-1,solution.T(2)-solution.T(1),length(solution.T)-1);
nt=input('')+1;

if n_gene==1
    figure
    hold on
    plot(solution.Xgrid{1},solution.PTX{nt},'k-','LineWidth',1.5)
    xlabel('Protein')
    ylabel('Probability')
    title(['Time ',num2str(solution.T(nt))])
    hold off
elseif n_gene==2
    figure
    mesh(solution.Xgrid{1},solution.Xgrid{2},solution.PTX{nt})
    xlabel('Protein_1')
    ylabel('Protein_2')
    zlabel('Probability')
    title(['Time ',num2str(solution.T(nt))])
    shading interp

    figure
    hold on
    PX1un = trapz(solution.x{2},solution.PTX{nt},2);
    plot(solution.x{1},PX1un/trapz(solution.x{1},PX1un),'k-','LineWidth',1.5)
    xlabel('Protein 1')
    ylabel('Probability')
    title(['Time ',num2str(solution.T(nt))])
    hold off

    figure
    hold on
    PX2un = trapz(solution.x{1},solution.PTX{nt},1);
    plot(solution.x{2},PX2un/trapz(solution.x{2},PX2un),'k-','LineWidth',1.5)
    xlabel('Protein 2')
    ylabel('Probability')
    title(['Time ',num2str(solution.T(nt))])
    hold off
elseif n_gene==3
    for i=1:n_gene
        PXiun = solution.PTX{nt};
        for j=1:n_gene
            if isequal(i,j)==0
                PXiun = trapz(solution.x{j},PXiun,j);
            end
        end
        PXi=squeeze(PXiun);
        figure
        hold on
        plot(solution.x{i},PXi/trapz(solution.x{i},PXi),'k-','LineWidth',1.5)
        xlabel(['Protein ',num2str(i)])
        ylabel('Probability')
        title(['Time ',num2str(solution.T(nt))])
        hold off
    end
    PX12=trapz(solution.x{3},solution.PTX{nt},3);
    figure
    mesh(solution.Xgrid{1}(:,:,1),solution.Xgrid{2}(:,:,1),PX12)
    xlabel('Protein_1')
    ylabel('Protein_2')
    zlabel('Probability')
    title(['Time ',num2str(solution.T(nt))])
    shading interp
    PX13=squeeze(trapz(solution.x{2},solution.PTX{nt},2));
    figure
    mesh(squeeze(solution.Xgrid{1}(:,1,:)),squeeze(solution.Xgrid{3}(:,1,:)),PX13)
    xlabel('Protein_1')
    ylabel('Protein_3')
    zlabel('Probability')
    title(['Time ',num2str(solution.T(nt))])
    shading interp
    PX23=squeeze(trapz(solution.x{1},solution.PTX{nt},1));
    figure
    mesh(squeeze(solution.Xgrid{2}(1,:,:)),squeeze(solution.Xgrid{3}(1,:,:)),PX23)
    xlabel('Protein_2')
    ylabel('Protein_3')
    zlabel('Probability')
    title(['Time ',num2str(solution.T(nt))])
    shading interp
elseif n_gene>3
    for i=1:n_gene
        PXiun = solution.PTX{nt};
        for j=1:n_gene
            if isequal(i,j)==0
                PXiun = trapz(solution.x{j},PXiun,j);
            end
        end
        PXi=squeeze(PXiun);
        figure
        hold on
        plot(solution.x{i},PXi/trapz(solution.x{i},PXi),'k-','LineWidth',1.5)
        xlabel(['Protein ',num2str(i)])
        ylabel('Probability')
        title(['Time ',num2str(solution.T(nt))])
        hold off    
    end
end

end