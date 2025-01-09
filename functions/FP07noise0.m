function  Sn = FP07noise0(params,fr)
    b1 = params(1);
    m1 = params(2);
    Sn=(10.^b1)*fr.^m1;
    

end