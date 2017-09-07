function cx=IF_FM_user(Xgrid,ci)
%%%
% cx=IF_FM_user(Xgrid,ci)
% user defined input function
%%%

% Number of genes
n_gene = length(Xgrid); % DO NOT MODIFY

% Feedback mechanism parameters
kAB = 1;
kBB = 0.01;
kCB = 1;
kAC = 1;
kBC = 1;
kCC = 0.01;

KAB = 100;
KBB = 1;
KCB = 50;
KAC = 0.01;
KBC = 0.01;
KCC = 0.5;

nAB = 2;
nBB = 2;
nCB = 4;
nAC = 2;
nBC = 4;
nCC = 2;

kiA=0.1;
induc=1.1;
KiA=1;


% Saving all c_i(x) in a cell
cx_cell=cell(n_gene,1); % DO NOT MODIFY

% We define the functions c_i, note that the variable x is a cell of
% dimension equal to the number of genes and each element x{i} is a
% multidimensional (n_gene) array which the same dimensions of your mesh.
% Use elementwise multiplication (.*) division (./) power (.^) ...
% in your input function expression.

cx_cell{1}=(kiA*(induc/KiA))/(1 + kiA*(induc/KiA));
cx_cell{2}=(kAB*(Xgrid{1}/KAB).^nAB + kCB*(Xgrid{3}/KCB).^nCB + kBB)./...
    (1 + kAB*(Xgrid{1}/KAB).^nAB + kBB*(Xgrid{2}/KBB).^nBB + kCB*(Xgrid{3}/KCB).^nCB);
cx_cell{3}=(kAC*(Xgrid{1}/KAC).^nAC + kCC*(Xgrid{3}/KCC).^nCC + kCB)./...
    (1 + kAC*(Xgrid{1}/KAC).^nAC + kBC*(Xgrid{2}/KBC).^nBC + kCC*(Xgrid{3}/KCC).^nCC);

% cx returns the input function (c_i(x)), i.e. the input function for regulated gene i=ci
cx=cx_cell{ci}; % DO NOT MODIFY
end
