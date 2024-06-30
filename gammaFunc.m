function sumg = gammaFunc(Delta,m,Nmax,m0)
%%%% This is a function to calculate the sum of gamma, where the specific value of m is obtained by solving gammaFunc(m)=0.


%%%% 0 order 
if Nmax == 0
        sumg = 0;
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
%%%% Plus the Delta.
    arrDel = Delta .^ (0:Nmax);
    garray(1,1) = m^2;
    garray1 = arrDel .* garray;
    sumg = sum(garray1)- m0^2;


end






end




%%%% The function of the dm
function dm = dmfunc(size,m)
dm = zeros(1,2*size);
for i = 1:size
    dm(i) = 1/(2*m*(-2*(size-1+i)) + (-2*(size-1+i))^2  );
    dm(i+size) = 1/(2*m*(2*i)+(2*i)^2);
end
end
