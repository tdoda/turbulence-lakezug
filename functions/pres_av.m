function [X,sX]=pres_av(pres0,x,pres,pint,fact,show)
%PRES_AV Compute mean and std of the variable x for each bin.
%
%   INPUTS:
%   pres0 (1*n double array): pressure profile [dbar]
%   x (1*n double array): variable to average
%   pres (1*nbin double array): middle pressure of the bins [dbar]
%   pint (double): bin size[dbar]
%   fact (double) [optional]: multiplication factor of std(x) above which x data is
%   not considered (i.e., outlier removal). Default value: 0.
%   show (boolean) [optional]: =True to display output. Default value: False.
%
%   OUTPUTS:
%   X (1*nbin double array): mean values of x for each bin.
%   sX (1*nbin double array): std values of x for each bin.

%
% T. Doda based on S. Piccolroaz, 18.12.2024
%%
if nargin<6
    show=false;
elseif nargin<5
    fact=0;
end
[pres0,is] = sort(pres0);
x = x(is);

X=nan(size(pres));
sX=nan(size(pres));

iend = find(isfinite(pres0),1,'last');
nn=find(pres>=pres0(iend),1,'first')-1;
if isempty(nn)
    nn = length(pres);
end
for i=1:nn
    x0=x(pres0>=pres(i)-0.5*pint & pres0<=pres(i)+0.5*pint);
    X(i)=nanmean(x0);
    sX(i)=nanstd(x0);

    if fact>0
        iiout=find(x0>X(i)+fact*sX(i) | x0<X(i)-fact*sX(i));
        x0(iiout)=NaN;
        X(i)=nanmean(x0);
        sX(i)=nanstd(x0);
        if show
            fprintf('\n Depth: %1.2f m, ndel/ntot: %d/%d',pres(i),numel(iiout),numel(x0));
        end
    end
end

end