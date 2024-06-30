function nu = mselect(Delta ,Nmax,m0)
%%%% This is the function to select and calculate the perturbate "m" by using fzero. Delta is the incoming energy and
%%%% can be set from 1e-6 to 100.The equation will be integrated into mfun as part of obtaining the coefficient array.
format long

switch m0
    case 0
        nu=fzero(@(m) gammaFunc(Delta,m,Nmax,m0),[1e-10,0.7]);
    case {1,-1}
        nu=fzero(@(m) gammaFunc(Delta,m,Nmax,m0),[0.6,1.4]);
    case {2,-2}
        nu=fzero(@(m) gammaFunc(Delta,m,Nmax,m0),[1.5,2.3]);
    otherwise
        disp("The angular momentum number is out of range!");
        return
end

