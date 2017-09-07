function cx=IF_FM_user(Xgrid,ci)
%%%
% cx=IF_FM_user(Xgrid,ci)
% user defined input function 
%%%

% Number of genes
n_gene = length(Xgrid); % DO NOT MODIFY

% Feedback mechanism parameters
H=2*eye(n_gene);
K=10*eye(n_gene);
epsilon=0.1*ones(n_gene,1);

% Saving all c_i(x) in a cell 
cx_cell=cell(n_gene,1); % DO NOT MODIFY

% We define the functions c_i, note that the variable x is a cell of
% dimension equal to the number of genes and each element x{i} is a 
% multidimensional (n_gene) array which the same dimensions of your mesh. 
% Use elementwise multiplication (.*) division (./) power (.^) ... 
% in your input function expression.  
for i=1:n_gene
    if H(i,i)>0
        cx_cell{i}=(epsilon(i)*Xgrid{i}.^H(i,i)+K(i,i)^H(i,i))./(Xgrid{i}.^H(i,i)+K(i,i)^H(i,i));
    elseif H(i,i)<0
        cx_cell{i}=(Xgrid{i}.^(-H(i,i))+epsilon(i)*K(i,i)^(-H(i,i)))./(Xgrid{i}.^(-H(i,i))+K(i,i)^(-H(i,i)));
    end
end

% cx returns the input function (c_i(x)), i.e. the input function for regulated gene i=ci 
cx=cx_cell{ci}; % DO NOT MODIFY
end
