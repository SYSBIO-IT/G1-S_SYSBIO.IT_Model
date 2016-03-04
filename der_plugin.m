function dX = der_plugin(t,X,n,N,Dd, Dk, Rb, Rt, Rk, Rw, Fd, Fr, Ud, Ur, Hr, G, CYF, G1S, CLNB, S1, Eta, MA, nX)
    
    nd  = n(1);
    nr  = n(2);
    nw1 = n(3);
    
    NA = N(1);
    NB = N(2);
    NC = N(3);
        
    P               = X(1);
    R               = X(2);
    V_cell          = P/G(7);
    V_nuc           = G(6)*V_cell;
    Cln3_tot        = G(8)*P;
    Cln3_ER         = X(3);
    Cln3Ydj1_ER     = X(4);
    Cln3Ydj1_n      = X(5);
    Cln3_n          = X(6);
    Cln3Far1        = X(7);
    Far1            = X(8);
    Cln3_c      	= Cln3_tot - (Cln3_ER + Cln3Ydj1_ER + Cln3Ydj1_n + Cln3_n + Cln3Far1);
    Swi6Swi4        = X(9);
    Swi6Swi4Whi5    = X(10);
    Swi6Mbp1        = X(11);
    Pd              = X(22+1:22+nd);
    Pr              = X(22+nd+1:22+nd+nr);
    Whi5_phospho	= X(22+1+nd+nr:nX);
    Whi5_n          = Whi5_phospho(1);
    
    Swi6            = MA(2) - (Swi6Swi4 + Swi6Swi4Whi5 + (1-Pr(1))*(NA + NC) + (1-Pd(1))*(NB + NC));
    Swi4            = MA(3) - (Swi6Swi4 + Swi6Swi4Whi5 + (1-Pr(1))*(NA + NC));
    Mbp1            = MA(5) - (Swi6Mbp1 + (1-Pd(1))*(NB + NC));

    Cln2            = X(12);
    Nrm1            = X(13);
    Cln1            = X(14);
    Clb6            = X(15);
    Clb5            = X(16);
    Sic1_n          = X(17);
    Clb6Sic1        = X(18);
    Clb5Sic1        = X(19);
    Clb6Sic1p       = X(20);
    Clb5Sic1p       = X(21);
    Sic1p           = X(22);
   
    aA      = Fr*Pr;
    aB      = Fd*Pd;
    aC      = aA + aB - aA*aB;
    f       = (NB+NC)*(Ud*Pd) + (NA+NC)*(Ur*Pr);

    dR	= G(1)*max(0,G(5)*P-R) - R/G(3);
    dP  = G(2)*R               - P/G(4);

    k4_tot	= CYF(5) + (CYF(6)-CYF(5))*(Cln3Ydj1_ER/CYF(7))^CYF(8)/(1+(Cln3Ydj1_ER/CYF(7))^CYF(8));
    if isnan(k4_tot)==1
        if Cln3Ydj1_ER < CYF(7)
            k4_tot = CYF(5);
        else
            k4_tot = CYF(6);
        end    
    end

    k10_tot = CYF(14)*(Cln3_n/Cln3Far1)^CYF(15)/(1+(Cln3_n/Cln3Far1)^CYF(15));
    if isnan(k10_tot)==1
        if Cln3_n < Cln3Far1
            k10_tot = 0;
        else
            k10_tot = CYF(14);
        end    
    end
    
    if Cln3_ER > CYF(4)
        k3_b = CYF(3);
    else
        k3_b = (CYF(3)/CYF(4))*Cln3_ER;
    end
    
    dCln3_ER        = - CYF(1)*Cln3_ER      + CYF(2)*Cln3_c                  - k3_b;
    dCln3Ydj1_ER    =   CYF(1)*Cln3_ER      - k4_tot*Cln3Ydj1_ER;
    dCln3Ydj1_n     = - CYF(9)*Cln3Ydj1_n   + k4_tot*Cln3Ydj1_ER;
    dCln3_n         =   CYF(9)*Cln3Ydj1_n   - (CYF(11)/V_nuc)*Cln3_n*Far1    + CYF(10)*Cln3Far1   - CYF(13)*Cln3_n     + CYF(12)*Cln3_c;
    dFar1           =   CYF(10)*Cln3Far1     - (CYF(11)/V_nuc)*Cln3_n*Far1    - k10_tot*Far1;
    dCln3Far1       = - CYF(10)*Cln3Far1     + (CYF(11)/V_nuc)*Cln3_n*Far1;
    dSwi6Swi4       = - G1S(1)*Swi6Swi4     + (G1S(2)/V_nuc)*Swi6*Swi4       - (G1S(8)/V_nuc)*Swi6Swi4*((NA+ NC)*Pr(1))         - (G1S(4)/V_nuc)*Swi6Swi4*Whi5_n	+ G1S(3)*Swi6Swi4Whi5;
    dSwi6Swi4Whi5   = - G1S(3)*Swi6Swi4Whi5 + (G1S(4)/V_nuc)*Swi6Swi4*Whi5_n - (G1S(9)/V_nuc)*Swi6Swi4Whi5*((NA	+ NC)*Pr(1));
    dSwi6Mbp1       = - G1S(6)*Swi6Mbp1     + (G1S(5)/V_nuc)*Swi6*Mbp1       - (G1S(7)/V_nuc)*Swi6Mbp1*((NB + NC)*Pd(1));

    f_Cln123    = 0.5*Cln1 + 0.5*Cln2 + Cln3_n;
    k20_tot     = G1S(10)*(f_Cln123/G1S(11))^G1S(12)/(1+(f_Cln123/G1S(11))^G1S(12));
    if isnan(k20_tot)==1        
        if f_Cln123 < G1S(11)
            k20_tot = 0;
        else
            k20_tot = G1S(10);
        end    
    end
    k20_tot = k20_tot/(100 + f);

	dPd = ((G1S(7)/V_nuc)*Swi6Mbp1*Dd  + (k20_tot/V_nuc)*Dk)*Pd;
    dPr = ((G1S(8)/V_nuc)*Swi6Swi4*Rb  + (G1S(9)/V_nuc)*Swi6Swi4Whi5*Rt + (k20_tot/V_nuc)*Rk + G1S(14)*Rw)*Pr;

    dWhi5_phospho	= (G1S(3)*Swi6Swi4Whi5 - (G1S(4)/V_nuc)*Swi6Swi4*Whi5_n)*[1;zeros(nw1,1)] + G1S(14)*(NA + NC)*Hr*Pr - diag(G1S(15:end))*Whi5_phospho;
        
    if Eta(2) <= NA*aA
        k23_tot = CLNB(3);     % Proper Cln2 expression by means of gene activation
    elseif t > Eta(7)
        k23_tot = CLNB(2);      % Random Cln2 expression
    else
        k23_tot = 0;
    end

    if Eta(3) <= NA*aA
        k25_tot = CLNB(6);     % Proper Nrm1 expression by means of gene activation
    elseif t > Eta(8)
        k25_tot = CLNB(5);      % Random Nrm1 expression
    else
        k25_tot = 0;
    end

    if Nrm1 < CLNB(10),       % If Nrm1 exceeds the threshold, it inhibits MBF and the activation of class C genes depends of only SBF activation
        Th_Cln1 = NC*aC;        
        Th_Clb5 = NB*aB;        
        Th_Clb6 = NC*aC;        
    else
        Th_Cln1 = NC*aA;
        Th_Clb5 = 0;
        Th_Clb6 = NC*aA;
    end

    if Eta(1) <= Th_Cln1
        k27_tot = CLNB(9);     % Proper Cln1 expression by means of gene activation
    elseif t > Eta(6)
        k27_tot = CLNB(8);      % Random Cln1 expression
    else
        k27_tot = 0;
    end

    dCln2	= k23_tot*V_nuc - CLNB(1)*Cln2;
    dNrm1   = k25_tot*V_nuc - CLNB(4)*Nrm1;
    dCln1   = k27_tot*V_nuc	- CLNB(7)*Cln1;

    if Eta(5) <= Th_Clb6
        k29_tot = CLNB(13);     % Proper Clb6 expression by means of gene activation
    elseif t > Eta(10)
        k29_tot = CLNB(12);      % Random Clb6 expression
    else
        k29_tot = 0;
    end

    if Eta(4) <= Th_Clb5
        k35_tot = CLNB(16);     % Proper Clb5 expression by means of gene activation
    elseif t > Eta(9)
        k35_tot = CLNB(15);      % Random Clb5 expression
    else
        k35_tot = 0;
    end

    k40_tot  = S1(9)*((0.1*Cln3_n + Cln1 + Cln2)/S1(10))^S1(11)/(1 + ((0.1*Cln3_n + Cln1 + Cln2)/S1(10))^S1(11));
    if isnan(k40_tot)==1
        if 0.1*Cln3_n + Cln1 + Cln2 < S1(10)
            k40_tot = 0;
        else
            k40_tot = S1(9);
        end    
    end

    k41_tot  = S1(12)*((0.1*Cln3_n + Cln1 + Cln2)/S1(13))^S1(14)/(1 + ((0.1*Cln3_n + Cln1 + Cln2)/S1(13))^S1(14));
    if isnan(k41_tot)==1
        if 0.1*Cln3_n + Cln1 + Cln2 < S1(13),
            k41_tot = 0;
        else
            k41_tot = S1(12);
        end    
    end

    k42_tot  = S1(15)*((0.1*Cln3_n + Cln1 + Cln2)/S1(16))^S1(17)/(1 + ((0.1*Cln3_n + Cln1 + Cln2)/S1(16))^S1(17));
    if isnan(k42_tot)==1
        if 0.1*Cln3_n + Cln1 + Cln2 < S1(16)
            k42_tot = 0;
        else
            k42_tot = S1(15);
        end    
    end

    dClb6       =   k29_tot*V_nuc   - CLNB(11)*Clb6             - (S1(2)/V_nuc)*Sic1_n*Clb6            + S1(1)*Clb6Sic1  - (S1(4)/V_nuc)*Clb6*Sic1p           + S1(3)*Clb6Sic1p;
    dClb5       =   k35_tot*V_nuc   - CLNB(14)*Clb5             - (S1(6)/V_nuc)*Sic1_n*Clb5            + S1(5)*Clb5Sic1  - (S1(8)/V_nuc)*Clb5*Sic1p           + S1(7)*Clb5Sic1p;
    dSic1_n     =   S1(5)*Clb5Sic1  - (S1(6)/V_nuc)*Sic1_n*Clb5 - (S1(2)/V_nuc)*Clb6*Sic1_n            + S1(1)*Clb6Sic1	 - (k40_tot/V_nuc)*Sic1_n*(Clb5+Clb6);
    dClb6Sic1   = - S1(1)*Clb6Sic1  + (S1(2)/V_nuc)*Sic1_n*Clb6	- (k41_tot/V_nuc)*Clb6Sic1*(Clb5+Clb6);
    dClb5Sic1   = - S1(5)*Clb5Sic1  + (S1(6)/V_nuc)*Sic1_n*Clb5 - (k42_tot/V_nuc)*Clb5Sic1*(Clb5+Clb6);
    dClb6Sic1p  = - S1(3)*Clb6Sic1p + (S1(4)/V_nuc)*Sic1p*Clb6  + (k41_tot/V_nuc)*Clb6Sic1*(Clb5+Clb6);
    dClb5Sic1p	= - S1(7)*Clb5Sic1p + (S1(8)/V_nuc)*Sic1p*Clb5  + (k42_tot/V_nuc)*Clb5Sic1*(Clb5+Clb6);
    dSic1p      = 	S1(7)*Clb5Sic1p	- (S1(8)/V_nuc)*Sic1p*Clb5	- (S1(4)/V_nuc)*Sic1p*Clb6             + S1(3)*Clb6Sic1p + (k40_tot/V_nuc)*Sic1_n*(Clb5+Clb6) - S1(18)*Sic1p;
    
    dX  = [dP;dR;dCln3_ER;dCln3Ydj1_ER;dCln3Ydj1_n;dCln3_n;dCln3Far1;dFar1;dSwi6Swi4;dSwi6Swi4Whi5;dSwi6Mbp1;dCln2;dNrm1;dCln1;dClb6;dClb5;dSic1_n;dClb6Sic1;dClb5Sic1;dClb6Sic1p;dClb5Sic1p;dSic1p;dPd;dPr;dWhi5_phospho];
end