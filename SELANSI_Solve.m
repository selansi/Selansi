% Reactions (i=1,...,n, with n=number of genes involved)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %       0 - c_i(x)km_i -> mRNA_i        %
        %  mRNA_i -    kx_i    -> mRNA_i + X_i  %
        %  mRNA_i - gamma_m_i  -> 0             %
        %     X_i -gamma_x_i(x)-> 0             %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% c_i(x)       := input function which collects the feedback mechanism
% km_i         := mRNA_i production rate
% kx_i         := X_i production rate
% gamma_m_i    := mRNA_i degradation rate
% gamma_x_i(x) := X_i degradation rate (can be a variable function)
function solution=SELANSI_Solve(name)
tic
% Obtaining the actual path
PathCurrent = pwd;

load(fullfile(PathCurrent,'DATA',name,'Reaction_data','parameters'))
load(fullfile(PathCurrent,'DATA',name,'Mesh_data','SL_parameters'))

path_forder_DR=fullfile(PathCurrent,'DATA',name);

% Using the semilagrangian method funtion to simulate the model
cd(fullfile(PathCurrent,'SELANSI_Files'))
fprintf('\n The semilagrangian method is running ... \n ')
[x,Xgrid,T,PX_sol]=friedsemilagran_nD(path_forder_DR,dato,SLdato);
cd(PathCurrent)
solution.T=T;
solution.x=x;
solution.Xgrid=Xgrid;
solution.PTX=PX_sol;

toc
save(fullfile(path_forder_DR,'Results','Solution.mat'),'solution','-v7.3');
%toc
cd(PathCurrent)
end

