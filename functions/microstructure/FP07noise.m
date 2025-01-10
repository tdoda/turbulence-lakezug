function  [Sn, Sn1, Sn2] = FP07noise(params,fr)
    b1 = params(1);
    m1 = params(2);
    fAA = params(3);
    Sn1=(10.^b1)*fr.^m1;
    
    %b2 = params(3);
    %m2 = params(4);
    %Sn2=(10.^b2)*fr.^m2;
    Sn2 = (1 + (fr/fAA).^8).^-2;
    
    Sn = Sn1.*Sn2;%min([Sn1,Sn2]')';
end