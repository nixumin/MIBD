function SameFun = MIBD_same(n,m,u)
% The distribution of IBD sharing when the two sampled individuals are taken from the same population
% Firstly edit by Xumin, Last edit by Wei
% Input
%    - n           - population size (haploid)
%    - m           - migration rate 
%    - u           - cutoff (the length is measured in Morgan)
% Output
%    - SameFun     - function for IBD sharing in same population
%--------------------------------------------------------------------------
SameFun = @(t)  ((4*n.*t.*u + 1).*(1./exp(t.*(4*m*n + (16*m^2*n^2 + 1).^(1/2) + 1)) - 1./exp(t.*(4*m*n - (16*m^2*n^2 + 1).^(1/2) + 1)) + ...
         (16*m^2*n^2 + 1)^(1/2)./exp(t.*(4*m*n - (16*m^2*n^2 + 1).^(1/2) + 1)) +...
         (16*m^2*n^2 + 1)^(1/2)./exp(t.*(4*m*n + (16*m^2*n^2 + 1)^(1/2) + 1))))./(exp(4*n.*t.*u).*(16*m^2*n^2 + 1)^(1/2));