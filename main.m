% Main
%% calculate the distribution of TMRCA
[T_Same1, T_Diff, T_Same2]=dsolve('Dx=-2*(1+2*n*m)*x+2*2*n*m*y+2','Dy=2*n*m*x-2*2*n*m*y+2*n*m*z','Dz=2*2*n*m*y-2*(1+2*n*m)*z+2','x(0)=0,y(0)=0,z(0)=0','t');

%% calculate the distribution of the IBD sharing in the same or across the population
syms t; syms u; syms n;

MIBD_Same = (exp(-4*u*n*t)+4*u*n*t*exp(-4*u*n*t))*diff(T_Same1,'t');
MIBD_Diff=(exp(-4*u*n*t)+4*u*n*t*exp(-4*u*n*t))*diff(T_Diff,'t');

simple(MIBD_Same);
simple(MIBD_Diff);

%% grid seaich to estimate parameters
[popsize, mig] = Grid_Search(N_lo, N_hi, N_step, m_lo, m_hi, m_step, cutoff, IBD_shar)