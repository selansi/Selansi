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

function SELANSI_Gnetmod(name)

% Obtaining the actual path
PathCurrent = pwd;

load(fullfile(PathCurrent,'DATA',name,'Reaction_data','parameters'));
load(fullfile(PathCurrent,'DATA',name,'Mesh_data','SL_parameters'))

path_forder_DR=fullfile(PathCurrent,'DATA',name);

if strcmp(dato.IF_Type,'Hill')==1 
    % Edit the variables
    open(fullfile(path_forder_DR,'Reaction_data','H.txt'))
    open(fullfile(path_forder_DR,'Reaction_data','K.txt'))
    open(fullfile(path_forder_DR,'Reaction_data','epsilon.txt'))
    open(fullfile(path_forder_DR,'Reaction_data','R_constants.txt'))
    fprintf('\n Write the constant values (substituting values by default) and SAVE.\n Push INTRO at the end to continue \n')
    pause

    % Loading the data to be saved in the variables
    dato.H=load(fullfile(path_forder_DR,'Reaction_data','H.txt'));
    dato.K=load(fullfile(path_forder_DR,'Reaction_data','K.txt'));
    n_gene=dato.n_gene;
    indreg=cell(n_gene,1);
    for i=1:n_gene
        indreg{i}=find(dato.H(i,:));
    end
    eps_size=zeros(1,n_gene);
    for i=1:n_gene
        eps_size(i)=length(indreg{i});
    end
    epsilon=cell(n_gene,1);
    aux_epsilon=importdata(fullfile(path_forder_DR,'Reaction_data','epsilon.txt'));
    epsilon{1}=aux_epsilon(1,:);
    if n_gene>1
        aux_row=2;
        for i=2:n_gene
            if eps_size(i)<=eps_size(1)
                epsilon{i}=aux_epsilon(aux_row,1:2^eps_size(i));
                aux_row=aux_row+1;
            else
                n_row=2^(eps_size(i)-eps_size(1));
                for ss=1:n_row
                    epsilon{i}(1+(n_row-1)*2^eps_size(1):(n_row)*2^eps_size(1))=aux_epsilon(aux_row +ss-1,:);
                end  
                aux_row=aux_row+n_row;
            end
        end
    end
    dato.epsilon=epsilon;
    dato.R_constants=load(fullfile(path_forder_DR,'Reaction_data','R_constants.txt'));
else
    %Edit the input function
    fprintf('\n Write your input function and SAVE \n Check that all parameters are correct \n Push INTRO at the end to continue \n');
    open(fullfile(path_forder_DR,'IF_FM_user.m'))
    pause
    
    % Edit the variables
    open(fullfile(path_forder_DR,'Reaction_data','R_constants.txt'))
    fprintf('\n Write the correct constant values (substituting values by default) and SAVE.\n Push INTRO at the end to continue \n')
    pause

    % Loading the data to be saved in the variables
    dato.R_constants=load(fullfile(path_forder_DR,'Reaction_data','R_constants.txt'));
end

save(fullfile(path_forder_DR,'Reaction_data','parameters.mat'),'dato');

fprintf('\n Your parameters are saved in %s \n',path_forder_DR)

end