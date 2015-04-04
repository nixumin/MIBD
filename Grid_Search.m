function [popsize, mig] = Grid_Search(N_lo, N_hi, N_step, m_lo, m_hi, m_step, cutoff, IBD_shar)
% This function is used to estimate the parameters using grid search method
% Firstly edit by Xumin, Last edit by Wei
% Input
%    - N_lo             - lower value of N
%    - N_hi             - higher value of N
%    - N_step           - the length of one step of N
%    - m_lo             - lower value of m
%    - m_hi             - higher value of m
%    - m_step           - the length of one step of m   
%    - cutoff           - the shortest length of IBD which we consider
%    - IBD_shar         - two-dimensional vector, the fisrt dimension is the IBD sharing when two
%                         individuals taken from the same population; the sedond dimension is the
%                         IBD sharing when two individuals are taken from the different populations
% Output
%    - popzise          - the effective population size 
%    - mig              - the migration rate
%
% test: 
% [N,m]=Grid_Search(5000,15000,10,0.0001,0.01,0.0001,0.02,[0.0046405 0.0003355])
% Output: N =10030; m=1.0000e-003.
%--------------------------------------------------------------------------
m = m_lo : m_step : m_hi;
n = N_lo : N_step : N_hi;

for i = 1 : length(n)
    
    for j = 1 : length(m)
        
        ls = Least_Square_fun([n(i),m(j)], IBD_shar, cutoff);
        LS(i,j) = ls;
        
    end
    
end

[lmin, k1] = min(LS);
[smin, k2] = min(lmin);

mig = (k2-1)*m_step + m_lo;
popsize = N_lo + (k1(k2)-1)*N_step;

