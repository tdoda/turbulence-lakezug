function SK=Tspec(spc,a,b,c,d,q)

if nargin>=4
    Xi=a;
    KB=b;
    K=c;
    k=d;
elseif nargin==3
    Xi=a(1);
    KB=a(2);
    K=b;
end

% See Sanchez et al. (2011) DOI: 10.1175/JPO-D-11-025.1
if spc=='K'  % Kraichnan
    phi=KB/sqrt(2*q);
    y=K/phi;
    f=y.*exp(-sqrt(3)*y); %Kraichnan
    SK=Xi/(2*k*KB)*sqrt(2*q)*f;
elseif spc=='B' % Batchelor
    phi=KB/sqrt(2*q);
    y=K/phi;
    f=y.*(exp(-y.^2/2)-y.*sqrt(pi()/2).*(1-erf(y/sqrt(2)))); %Batchelor
    SK=Xi/(2*k*KB)*sqrt(2*q)*f;
end
end
