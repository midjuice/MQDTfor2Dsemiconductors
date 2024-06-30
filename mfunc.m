function [Marray,m,M] = mfunc(E,m0)
%%%% The calculation for coefficients of M_{m}, Nmax is the max order of
%%%% the perturbation. E is the entrance energy and can be select from
%%%% 0 to 90. m0 is the angular momentum and can be select as 0,+-1,+-2.

Nmax = 10;
Delta = E+1e-6;
m = mselect(Delta,Nmax,m0);





%%%% 0 order 
if Nmax == 0
        Marray = 1;
        garray = 0;
%%%% 1-Nmax order
else
    Marray = zeros(Nmax+1,2*Nmax);
    Marray(1,1) = 1; 
    garray = zeros(1,Nmax+1);
    for i = 2:Nmax+1
        index = nonzeros(Marray(i-1,:));
        [Mmid,garray(i)] = contifunc(index,m,Nmax);
%%%% Calculation of sum_{i=1}^{n}\gamma_{n}M_{n-i}
        sum2 = zeros(1,2*Nmax);
        for ii = 1:(i-1)            
            sum2 = sum2 + garray(ii)*Marray(i-ii+1,:);
        end
        sum1 = dmfunc(Nmax,m).*sum2;
        Marray(i,:) = Delta^(i-1)*(Mmid+sum1);
    end


end
Marray(1,1)=0;
newColumn = zeros(Nmax+1, 1);
newColumn(1,1)=1;
Marray = [Marray(:, 1:Nmax), newColumn, Marray(:, Nmax+1:end)];
M = sum(Marray);


end




%%%% The function of the dm
function dm = dmfunc(size,m)
dm = zeros(1,2*size);
for i = 1:size
    dm(i) = 1/(2*m*(-2*(size-1+i)) + (-2*(size-1+i))^2  );
    dm(i+size) = 1/(2*m*(2*i)+(2*i)^2);
end
end











































