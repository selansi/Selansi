
function IC=IC_Function(x,Xgrid)
%%% 
% IC=IC_Function(x,Xgrid)
% IMPORTANT: The initial Condition is approached by a Gaussian density and must be a 
% multidimensional array of dimension equal to Xgrid{i} with i=1,...,n_gene.
%%%

% Number of genes
n_gene = length(Xgrid);

% Parameters of the Gaussian density funcion
x0g=10*ones(1,n_gene);  % Mean of the Gaussian density
sigmag=5*ones(1,n_gene); % Standard deviation of the Gaussian density 

gausker=0;
for i=1:n_gene
    gausker = gausker+((Xgrid{i}-x0g(i))/sigmag(i)).^2;
end
PX0un=exp(-gausker/2);

% Norm of the initial condition
auxnor0=PX0un;
for i=1:n_gene
    auxnor0 = trapz(x{i},auxnor0);
end

% Normalization of the initial condition to be a density function
IC=PX0un/auxnor0;
end