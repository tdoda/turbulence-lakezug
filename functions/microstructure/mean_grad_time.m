function grX=mean_grad_time(pres0,x,pres,pint,timef)
%MEAN_GRAD_TIME Computes the average temporal gradient of x for each bin.
%
%   INPUTS:
%   pres0 (1*n double array): pressure profile [dbar]
%   x (1*n double array): variable for which the averaged gradient is
%   computed
%   pres (1*nbin double array): pressure of the bins [dbar]
%   pint (double): bin size [dbar]
%   timef (1*n double array): time values [s]
%
%   OUTPUTS:
%   grX (1*nbin double array): temporal gradient of x in each bin
% 
% T. Doda based on S. Piccolroaz, 19.12.2024
%%  
    grX=nan(size(pres));
    sgrX=nan(size(pres));

    [pres0,is] = sort(pres0);
    x = x(is);

    nn=find(pres>=pres0(end),1,'first')-1;
    if isempty(nn)
        nn = length(pres);
    end

    for i=1:nn
         jj=find(pres0>=pres(i)-0.5*pint & pres0<=pres(i)+0.5*pint);
         if length(jj)>10
             x0=timef(jj);
             y0=x(jj);
             
             [grX(i),sgrX(i)] = slope(-x0,y0);
             
 
         end
         
    end
end