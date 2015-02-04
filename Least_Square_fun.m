function LeastSqr = Least_Square_fun( params, RealIBD, cutoff )
% This funcion is used to calculate least square of the IBD sharing 
% sharing and the theoretical IBD sharing
% Firstly edit by Xumin, Last edit by Wei
% Input
%    - params           - population size and migration rate 
%    - RealIBD          - real IBD sharing in real data
%    - cutoff           - cutoff
% Output
%    - LeastSqr         - least square of IBD sharing
%--------------------------------------------------------------------------
n = params(1);  % population size
m = params(2);  % migration rate

syms t;

% obtain the theoretical expression of IBD sharing
DiffFun = MIBD_diff(n, m, cutoff);
SameFun = MIBD_same(n, m, cutoff);
DiffSh = quadgk(DiffFun, 0, Inf);
SameSh =quadgk(SameFun,0,Inf);

x = RealIBD;
Error_Vector = (SameSh - x(1)).^2 + (DiffSh - x(2)).^2;
LeastSqr = Error_Vector;