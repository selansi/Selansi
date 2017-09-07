function cx=IF_FM_user(Xgrid,ci)
%%%
% cx=IF_FM_user(Xgrid,ci)
% user defined input function 
%%%

% Number of genes
n_gene = length(Xgrid); % DO NOT MODIFY

% Feedback mechanism parameters
A=5;
B=10;

% Saving all c_i(x) in a cell 
cx_cell=cell(n_gene,1); % DO NOT MODIFY

% We define the functions c_i, note that the variable x is a cell of
% dimension equal to the number of genes and each element x{i} is a 
% multidimensional (n_gene) array which the same dimensions of your mesh. 
% Use elementwise multiplication (.*) division (./) power (.^) ... 
% in your input function expression.  
cx_cell{1}=1./(1+exp(A-B*Xgrid{1}));


% cx returns the input function (c_i(x)), i.e. the input function for regulated gene i=ci 
cx=cx_cell{ci}; % DO NOT MODIFY
end

