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
close all
clear all

% Number of genes considered
n_gene=input('Write the number (integer) of genes considered \n');
fprintf('You have selected %d genes \n',n_gene);

% Degradation functions are considered contants
DF_Type=0;
%%%%%%%%%%%%%%%%% this option can be variable in the future %%%%%%%%%%%%%%%
% % DF_Type=input('\n Write 0 if the degradation functions are constant or \n write 1 if the degradation functions are variable \n');
% % if DF_Type==0
% %     fprintf('You have selected constant degradation functions \n');
% % elseif DF_Type==1
% %     fprintf('You have selected variable degradation functions \n');
% % else
% %     fprintf('ERROR!!! Your degradation option is not correct \n')
% % end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Obtain the actual path
PathCurrent = pwd;

% Create a forder to save data and results
Folder_name=input('Write the folder name for saving data and results: \n','s');
fprintf('The name for your forder is : %s \n',Folder_name);
path_forder_DR=fullfile(PathCurrent,'DATA',Folder_name);
mkdir(path_forder_DR);
mkdir(fullfile(path_forder_DR,'Reaction_data'));
mkdir(fullfile(path_forder_DR,'Mesh_data'));
mkdir(fullfile(path_forder_DR,'Results'));

% Feedback mechanism constants cx=IF_FeedbackMechanism(IF_Type,n,x,H,K,epsilon)
IFT_define=input('Choose your feedback mechanism input function \n Insert 1 if you want to use a predefined Hill function \n Insert 0 if you define your input function \n');
if IFT_define==1
    IF_Type = 'Hill'; % this option can be variable in the future
else
    IF_Type = 'user';
end


if strcmp(IF_Type,'Hill')==1 
    % Hill coefficients 
    H=ones(n_gene,n_gene);
    FID=fopen(fullfile(path_forder_DR,'Reaction_data','H.txt'),'w+');
    for i=1:n_gene
        for j=1: n_gene
            fprintf(FID,'%g ',H(i,j));
        end
        fprintf(FID,'\n');
    end
    fclose(FID);
    % Edit the variables
    open(fullfile(path_forder_DR,'Reaction_data','H.txt'))
    fprintf('\n Write the correct Hill coefficients (substitute values by default) and SAVE\n Push INTRO at the end to continue \n')
    pause

     % Loading the data to be saved in the variables
    dato.H=load(fullfile(path_forder_DR,'Reaction_data','H.txt'));
    indreg=cell(n_gene,1);
    for i=1:n_gene
        indreg{i}=find(dato.H(i,:));
    end
    % Equilibrium constants
    K=zeros(n_gene,n_gene);
    for i=1:n_gene
        K(i,indreg{i})=1;
    end
    FID=fopen(fullfile(path_forder_DR,'Reaction_data','K.txt'),'w+');
    for i=1:n_gene
        for j=1: n_gene
            fprintf(FID,'%f ',K(i,j));
        end
        fprintf(FID,'\n');
    end
    fclose(FID);
       
    % Leakage rates
    eps_size=zeros(1,n_gene);
    for i=1:n_gene
        eps_size(i)=length(indreg{i});
    end
    epsilon=cell(n_gene,1);
    for i=1:n_gene
        if eps_size==1
            epsilon{i}=[0.1 1];
        else
            epsilon{i}=[0.1, 0.1*ones(1,2^eps_size(i)-2),1];
        end
    end
    FID=fopen(fullfile(path_forder_DR,'Reaction_data','epsilon.txt'),'w+');
    for i=1:n_gene
        for j=1: 2^eps_size(i)
            fprintf(FID,'%f ',epsilon{i}(j));
        end
        fprintf(FID,'\n');
    end
    fclose(FID);
    
elseif strcmp(IF_Type,'user')==1
    copyfile(fullfile(pwd,'SELANSI_Files','InputFunction','IF_FM_user.m'),path_forder_DR)
end

% Matrix (nx4) of input variables (reaction contants)
R_constants=ones(n_gene,4);
R_constants(:,3)=25*ones(n_gene,1);
R_constants(:,2)=75*ones(n_gene,1);
R_constants(:,1)=2*ones(n_gene,1);
FID=fopen(fullfile(path_forder_DR,'Reaction_data','R_constants.txt'),'w+');
for i=1:n_gene
    fprintf(FID,'%f %f %f %f \n',R_constants(i,:));
end
fclose(FID);

% Mesh for the semilagrangian method
% Number of proteins mesh, [Xi_min, Xi_max, deltaxi]
Prot_mesh=zeros(n_gene,3);
Prot_mesh(:,3)=300*ones(n_gene,1);
Prot_mesh(:,2)=30*ones(n_gene,1);
FID=fopen(fullfile(path_forder_DR,'Mesh_data','Prot_mesh.txt'),'w+');
for i=1:n_gene
    fprintf(FID,'%f %f %f \n',Prot_mesh(i,:));
end
fclose(FID);

% Mesh of the semilagrangian method
Time_mesh=[0 5 10 100];
FID=fopen(fullfile(path_forder_DR,'Mesh_data','Time_mesh.txt'),'w+');
fprintf(FID,'%f %f %f %f \n',Time_mesh);
fclose(FID);

% Initial condition for the semiLagrangian method
copyfile(fullfile(pwd,'SELANSI_Files','InitialCondition','IC_Function.m'),path_forder_DR)
fprintf('\n Write your Initial Condition and SAVE \n Check that it is correct \n Push INTRO at the end to continue \n');
open(fullfile(path_forder_DR,'IC_Function.m'))
pause

if strcmp(IF_Type,'Hill')==1 
    % Edit the variables
    open(fullfile(path_forder_DR,'Reaction_data','K.txt'))
    open(fullfile(path_forder_DR,'Reaction_data','epsilon.txt'))
    open(fullfile(path_forder_DR,'Reaction_data','R_constants.txt'))
    open(fullfile(path_forder_DR,'Mesh_data','Prot_mesh.txt'))
    open(fullfile(path_forder_DR,'Mesh_data','Time_mesh.txt'))
    fprintf('\n Write the constant values (substituting values by default) and SAVE\n Push INTRO at the end to continue \n')
    pause

     % Loading the data to be saved in the variables
    dato.K=load(fullfile(path_forder_DR,'Reaction_data','K.txt'));
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
    dato.n_gene=n_gene;
    dato.DF_Type=DF_Type;
    dato.IF_Type=IF_Type;
    SLdato.Prot_mesh=load(fullfile(path_forder_DR,'Mesh_data','Prot_mesh.txt'));
    SLdato.Time_mesh=load(fullfile(path_forder_DR,'Mesh_data','Time_mesh.txt'));
    save(fullfile(path_forder_DR,'Reaction_data','parameters.mat'),'dato');
    save(fullfile(path_forder_DR,'Mesh_data','SL_parameters.mat'),'SLdato');
else 
    fprintf('\n Write your input function and SAVE \n Check that all parameters are correct \n Push INTRO at the end to continue \n');
    open(fullfile(path_forder_DR,'IF_FM_user.m'))
    pause
     % Edit the variables
    open(fullfile(path_forder_DR,'Reaction_data','R_constants.txt'))
    open(fullfile(path_forder_DR,'Mesh_data','Prot_mesh.txt'))
    open(fullfile(path_forder_DR,'Mesh_data','Time_mesh.txt'))
    fprintf('\n Write the  constant values (substituting values by default) and SAVE.\n Push INTRO at the end to continue \n')
    pause

     % Loading the data to be saved in the variables
    dato.R_constants=load(fullfile(path_forder_DR,'Reaction_data','R_constants.txt'));
    dato.n_gene=n_gene;
    dato.DF_Type=DF_Type;
    dato.IF_Type=IF_Type;
    SLdato.Prot_mesh=load(fullfile(path_forder_DR,'Mesh_data','Prot_mesh.txt'));
    SLdato.Time_mesh=load(fullfile(path_forder_DR,'Mesh_data','Time_mesh.txt'));
    save(fullfile(path_forder_DR,'Reaction_data','parameters.mat'),'dato');
    save(fullfile(path_forder_DR,'Mesh_data','SL_parameters.mat'),'SLdato');
end

fprintf('\n Your parameters are saved in %s \n',path_forder_DR)

