function Figures (T, X, AG, MA, Fd, Fr, N, n)
    
    nX      = size(X,1);
    
    Swi6Swi4        = X(9,:);
    Swi6Swi4Whi5    = X(10,:);
    Swi6Mbp1        = X(11,:);
    Cln2            = X(12,:);
    Nrm1            = X(13,:);
    Cln1            = X(14,:);
    Clb6            = X(15,:);
    Clb5            = X(16,:);
    Sic1_n          = X(17,:);
    Clb6Sic1        = X(18,:);
    Clb5Sic1        = X(19,:);
    Clb6Sic1p       = X(20,:);
    Clb5Sic1p       = X(21,:);
    Sic1p           = X(22,:);
    Sic1_c          = MA(6) - (Sic1_n + Clb5Sic1 + Clb5Sic1p + Clb6Sic1 + Clb6Sic1p + Sic1p);

    Pd              = X(22+1:22+n(1),:);
    Pr              = X(22+n(1)+1:22+n(1)+n(2),:);

    Whi5_phospho	= X(22+n(1)+n(2)+1:nX,:);
    Whi5_n          = Whi5_phospho(1,:);

    aA              = Fr*Pr;
    aB              = Fd*Pd;
    aC              = aA + aB - aA.*aB;

    Whi5_c          = MA(4) - (Swi6Swi4Whi5 + sum(Whi5_phospho,1) + (N(1) + N(3))*(1-(Pr(1,:))-aA));    
    Mbp1            = MA(5) - (Swi6Mbp1 + (1-(Pd(1,:)))*(N(2) + N(3)));
    Swi6            = MA(2) - (Swi6Swi4 + Swi6Swi4Whi5 + (1-(Pr(1,:)))*(N(1) + N(3)) + (1-(Pd(1,:)))*(N(2) + N(3)));
    Swi4            = MA(3) - (Swi6Swi4 + Swi6Swi4Whi5 + (1-(Pr(1,:)))*(N(1) + N(3)));

    figure
    plot(T,Swi6,'r',T,Swi4,'k',T,Swi6Swi4,'b',T,Swi6Swi4Whi5,'m',T,Whi5_n,'g'),title('Swi6 in RED, Swi4 in BLACK, free unbound Whi5 in GREEN, Swi6Swi4 in BLUE, Swi6Swi4Whi5 in PURPLE')
    figure
    plot(T,Swi6,'r',T,Mbp1,'k',T,Swi6Mbp1,'b'),title('Swi6 in RED, Mbp1 in BLACK, Swi6Mbp1 in BLUE')
    figure
    plot(T,Whi5_n,'r',T,Whi5_c,'k'),title('free Whi5_n in RED, Whi5_{cyt} in BLACK')
    figure
    plot(T,Clb5+Clb6,'b',T,Sic1_n,'k',T,Clb5Sic1+Clb6Sic1,'r',T,Clb5Sic1p+Clb6Sic1p,'g',T,Sic1p,'m',T,Sic1_c,'c'),title('Clb5,6 in BLUE, Sic1_{nuc} in BLACK, Clb5,6Sic1 in RED, Clb5,6Sic1p in GREEN, Sic1p in PURPLE, Sic1_{cyt} in CYAN')
    figure
    plot(T,Cln1,'r',T,Cln2,'k',T,Nrm1,'m',T,Clb5+Clb5Sic1+Clb5Sic1p,'g',T,Clb6+Clb6Sic1+Clb6Sic1p,'b'),title('Sum of all Clb5 in GREEN, Sum of all Clb6 in BLUE, Nrm1 in PURPLE, Cln1 in RED, Cln2 in BLACK')
    figure
    plot(T,(N(1)*aA + N(2)*aB + N(3)*aC)/(N(1)+N(2)+N(3)),'r',T,(MA(4)-Whi5_c)./MA(4),'b'),title('Percentage of activated genes in RED, percentage of nuclear Whi5 in BLUE')
    figure
    plot(T,Cln2,'r',T,Sic1_n,'y',T,Clb5,'g',T,100*AG,'b'),title('Cln2 in RED, free Sic1 in YELLOW, Clb5 in GREEN, % activated genes in BLUE')

end