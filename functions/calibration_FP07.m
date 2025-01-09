function uT = calibration_FP07(timef0,Tf0,times,Ts,fsf,fss)

timef = timef0(1:fsf/fss:end);
Tf = Tf0(1:fsf/fss:end);
x = interp1(timef, Tf, times(:));
y = Ts(:);
ii = find(isfinite(x) & isfinite(y));
x = x(ii);
y = y(ii);


pp = polyfit(x,y,1);
uT = pp(1)*Tf0 + pp(2);
