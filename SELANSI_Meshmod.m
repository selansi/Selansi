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

function SELANSI_Meshmod(name)

% Obtaining the actual path
PathCurrent = pwd;

load(fullfile(PathCurrent,'DATA',name,'Reaction_data','parameters'))
load(fullfile(PathCurrent,'DATA',name,'Mesh_data','SL_parameters'))

path_forder_DR=fullfile(PathCurrent,'DATA',name);

% Edit initial condition
fprintf('\n Write your Initial Condition and save it \n Check that it is correct \n Push INTRO at the end to continue \n');
open(fullfile(path_forder_DR,'IC_Function.m'))
pause

open(fullfile(path_forder_DR,'Mesh_data','Prot_mesh.txt'))
open(fullfile(path_forder_DR,'Mesh_data','Time_mesh.txt'))
fprintf('\n Write the constant values (substituting values by default) and SAVE \n Push INTRO at the end to continue \n')
pause

SLdato.Prot_mesh=load(fullfile(path_forder_DR,'Mesh_data','Prot_mesh.txt'));
SLdato.Time_mesh=load(fullfile(path_forder_DR,'Mesh_data','Time_mesh.txt'));

save(fullfile(path_forder_DR,'Mesh_data','SL_parameters.mat'),'SLdato');

fprintf('\n Your parameters are saved in %s \n',path_forder_DR)

end
