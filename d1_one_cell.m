function [T,X,T1,T2,Ps,AG,MA,Vs,alpha] = d1_one_cell(n, N, Dd, Dk, Rb, Rt, Rk, Rw, Fd, Fr, Ud, Ur, Hr)
    
    nd  = n(1);
    nr  = n(2);
    nw1 = n(3);
    
    %% Growth parameters of the G1/S transition model
    K1          = 1;            % [min^-1]      Ribosome production rate
    K2          = 558;          % [aa/rib/min]  Protein production rate per ribosomes
    D1          = 4000;         % [min]         Ribosome degradation time constant
    D2          = 3000;         % [min]         Protein degradation time constant
    rho         = 1.34e-5;   % [rib/aa]      Ribosome-over-protein ratio
    h           = 0.07;         %               Nuclear-over-cell volume ratio
    H           = 7.1e23;       % [aa/L]        Protein-over-cell volume ratio
    k0          = 3e-8;          % [molecules/aa]        Total Cln3 over protein ratio
    G = [K1, K2, D1, D2, rho, h, H, k0];
    
	%% Cln3 production and nuclear important parameters
    k1         = 100;           % [min^-1]              ON coeff. of Cln3 + Ydj1 --> Cln3Ydj1 in the ER, including the constant contribute at saturation of Ydj1 
    k2         = 1;             % [min^-1]              Cln3 diffusion coeff. from cytoplasm into ER
    k3         = 200;           % [min^-1]              Cln3 diffusion coeff. from ER into cytoplasm (mean value)
    Theta_3    = 10;            % [molecules]           Threshold for Cln3 diffusion rate from ER into cytoplasm
    k4_low     = 3;             % [min^-1]              Cln3Ydj1 diffusion from the ER into the nucleus for low values of ER-Cln3Ydj1
    k4_high    = 200;           % [min^-1]              Cln3Ydj1 diffusion from the ER into the nucleus for high values of ER-Cln3Ydj1
    Theta_4    = 180;           % [molecules]           Threshold for low/high values of Cln3Ydj1 diffusion from the ER into nucleus
    n4         = 30;            %                       Hill coeff. for Cln3Ydj1 diffusion from the ER into nucleus
    k5         = 100;           % [min^-1]              OFF coeff. of Cln3 + Ydj1 <-- Cln3Ydj1 in the nucleus
    k6         = 25;            % [min^-1]              OFF coeff. of Cln3 + Far1 <-- Cln3Far1 in the nucleus
    k7         = 3.2e-15;       % [(L/molecules)/min]	ON coeff. of Cln3 + Far1 --> Cln3Far1 in the nucleus
    k8         = 0.0015;        % [min^-1]              Ydj1-independent Cln3 diffusion coeff. from cytoplasm into nucleus
    k9         = 0.4;           % [min^-1]              Cln3 diffusion coeff. from nucleus into cytoplasm
    k10        = 1;             % [min^-1]              Maximum degradation rate for Far1 after the onset of Timer T1
    n10        = 10;            %                       Hill coeff. for Far1 degradation rate
    CYF = [k1, k2, k3,Theta_3, k4_low, k4_high, Theta_4, n4, k5, k6, k7, k8, k9, k10, n10];
    
    %% Input parameters for G1/S regulon activation
    k11        = 25;            % [min^-1]              OFF coeff. of Swi6 + Swi4 <-- Swi6Swi4
    k12        = 1.6e-15;       % [(L/molecules)/min]   ON coeff. of Swi6 + Swi4 --> Swi6Swi4
    k13        = 2.5;           % [min^-1]              OFF coeff. of Swi6Swi4 + Whi5 <-- Swi6Swi4Whi5
    k14        = 1.6e-14;       % [(L/molecules)/min]   ON coeff. of Swi6Swi4 + Whi5 --> Swi6Swi4Whi5
    k15        = 1.6e-15;       % [(L/molecules)/min]   ON coeff. of Swi6 + Mbp1 --> Swi6Mbp1
    k16        = 25;            % [min^-1]              OFF coeff. of Swi6 + Mbp1 <-- Swi6Mbp1
    k17        = 1e-17;         % [(L/molecules)/min]	 Binding coeff. to MBF-affine genes of Swi6Mbp1
    k18        = 1e-17;         % [(L/molecules)/min]	 Binding coeff. to SBF-affine genes of Swi6Swi4
    k19        = 1e-17;         % [(L/molecules)/min]	 Binding coeff. to SBF-affine genes of Swi6Swi4Whi5
    k20        = 2.4e-12;       % [L*molecules/min]     Maximum phosphorylation rate for Swi6 and Whi5
    Theta_20   = 1e3;           % [molecules]           Threshold for f(Cln1,2,3) in Swi6 and Whi5 phosphorylation
    n20        = 2;             %                       Hill coeff. for f(Cln1,2,3) in Swi6 and Whi5 phosphorylation
    N20        = 100;           % [molecules]           Minimum leakage factor for Swi6 and Whi5 phosphorylation
    k21_0      = 0.07;          % [min^-1]              Diffusion of Whi5 from nucleus into the cytoplasm
    k21p       = 0.7;           % [min^-1]              Diffusion of phosphorylated Whi5 from nucleus into the cytoplasm
    k21_tot    = [k21_0 k21p*ones(1,nw1)];                      
	al_W       = 5;             % [min^-1]              Dissociation rate constant for Whi5
    G1S =[k11, k12, k13, k14, k15, k16, k17, k18, k19, k20, Theta_20, n20, N20, al_W, k21_tot];

    %% Input parameters for Cln1,2, Clb5,6 and Nrm1 dynamics
    k22        = 0.01;          % [min^-1]             Degradation rate of Cln2
    k23_low    = 7.2e14;        % [(molecules/L)/min]  Low Cln2 production rate (without SBF activation)
    k23_high   = 4.3e15;        % [(molecules/L)/min]	Proper Cln2 production rate (due to SBF activation)
    k24        = 0.05;          % [min^-1]             Degradation rate of Nrm1
    k25_low    = 1.8e13;        % [(molecules/L)/min]	Weak Nrm1 production rate (without SBF activation)
    k25_high   = 1.8e15;        % [(molecules/L)/min]  Proper Nrm1 production rate (due to SBF activation)
    k26        = 0.01;          % [min^-1]             Degradation rate of Cln1
    k27_low    = 1.8e15;        % [(molecules/L)/min]	Weak Cln1 production rate (without SBF/MBF activation)
    k27_high   = 1.08e16;       % [(molecules/L)/min]  Proper Cln1 production rate (due to SBF/MBF activation)
    Theta_Nrm1 = 50;            % [molecules]          Number of molecules for Nrm1 to inhibit MBF
    k28        = 0.05;          % [min^-1]             Degradation rate of Clb6
    k29_low    = 1.26e13;       % [(molecules/L)/min]	Weak Clb6 production rate (without MBF activation)
    k29_high   = 4.5e15;        % [(molecules/L)/min]  Proper Clb6 production rate (due to MBF activation)
    k34        = 0.05;          % [min^-1]             Degradation rate of Clb5
    k35_low    = 1.26e13;       % [(molecules/L)/min]	Weak Clb5 production rate (without SBF/MBF activation)
    k35_high   = 6.3e15;        % [(molecules/L)/min]  Proper Clb5 production rate (due to SBF/MBF activation)
    CLNB = ([k22 k23_low k23_high k24 k25_low k25_high k26 k27_low k27_high Theta_Nrm1 k28 k29_low k29_high k34 k35_low k35_high]);
    
    %% Input parameters for Sic1 function
    k30        = 0.35;          % [min^-1]             OFF coeff. of Clb6 + Sic1 <-- Clb6Sic1
    k31        = 1e-15;         % [(L/molecules)/min]	ON coeff. of Clb6 + Sic1 --> Clb6Sic1
    k32        = 50;            % [min^-1]             OFF coeff. of Clb6 + Sic1P <-- Clb6Sic1P
    k33        = 1e-16;         % [(L/molecules)/min]	ON coeff. of Clb6 + Sic1P --> Clb6Sic1P
    k36        = 0.35;          % [min^-1]             OFF coeff. of Clb5 + Sic1 <-- Clb5Sic1
    k37        = 1e-15;         % [(L/molecules)/min]	ON coeff. of Clb5 + Sic1 --> Clb5Sic1
    k38        = 50;            % [min^-1]             OFF coeff. of Clb5 + Sic1P <-- Clb5Sic1P
    k39        = 1e-16;         % [(L/molecules)/min]	ON coeff. of Clb5 + Sic1P --> Clb5Sic1P
    k40        = 1.3e-18;       % [(L/molecules)/min]	Sic1 --> Sic1P phosphorylation rate, by means of g(Cln1,2,3) AND (Clb5 + Clb6)
    Theta_40   = 150;           % [molecules]          Threshold for g(Cln1,2,3) in Sic1 --> Sic1P phosphorylation
    n40        = 5;             %                      Hill coeff. for g(Cln1,2,3) in Sic1 --> Sic1P phosphorylation
    k41        = 1.7e-16;       % [(L/molecules)/min]	Clb6Sic1 --> Clb6Sic1P phosphorylation rate, by means of g(Cln1,2,3) AND (Clb5 + Clb6)
    Theta_41   = 150;           % [molecules]          Threshold for g(Cln1,2,3) in Clb6Sic1 --> Clb5Sic1P phosphorylation
    n41        = 5;             %                      Hill coeff. for g(Cln1,2,3) in Clb6Sic1 --> Clb5Sic1P phosphorylation
    k42        = 1.7e-16;       % [(L/molecules)/min]	Clb5Sic1 --> Clb5Sic1p phosphorylation rate, by means of (Cln1 + Cln2 + 0.1*Cln3_n) AND (Clb5+Clb6)
    Theta_42   = 150;           % [molecules]          Threshold for g(Cln1,2,3) in Clb5Sic1 --> Clb5Sic1P phosphorylation
    n42        = 5;             %                      Hill coeff. for g(Cln1,2,3) in Clb5Sic1 --> Clb5Sic1P phosphorylation
    k43        = 0.7;           % [min^-1]             Diffusion of Sic1P out of the nucleus
    S1 = ([k30 k31 k32 k33 k36 k37 k38 k39 k40 Theta_40 n40 k41 Theta_41 n41 k42 Theta_42 n42 k43]);

    
    %% Proper activation orders and week activation times
	NCln1 = 5;      % Average value for the activation order of Cln1 within the class C genes
    NCln2 = 30;     % Average value for the activation order of Cln2 within the class A genes
    NNrm1 = 110;    % Average value for the activation order of Nrm1 within the class A genes
    NClb5 = 45;     % Average value for the activation order of Clb5 within the class B genes
    NClb6 = 20;     % Average value for the activation order of Clb6 within the class C genes
    tCln1 = 30;     % Average value for the week activation time of Cln1
    tCln2 = 30;     % Average value for the week activation time of Cln1
    tNrm1 = 60;     % Average value for the week activation time of Nrm1
    tClb5 = 60;     % Average value for the week activation time of Clb5
    tClb6 = 60;     % Average value for the week activation time of Clb6
    Eta = ([NCln1 NCln2 NNrm1 NClb5 NClb6 tCln1 tCln2 tNrm1 tClb5 tClb6]); 
    
    %% Molecules amounts
    Far1_tot        = 240;      % [molecules]   Total amount of Far1
    Swi6_tot        = 300;      % [molecules]   Total amount of Swi6
    Swi4_tot        = 200;      % [molecules]   Total amount of Swi4
    Mbp1_tot        = 200;      % [molecules]   Total amount of Mbp1
    Sic1_tot        = 300;      % [molecules]   Total amount of Sic1
    Whi5_tot        = 200;      % [molecules]   Total amount of Whi5
    MA = ([Far1_tot Swi6_tot Swi4_tot Whi5_tot Mbp1_tot Sic1_tot]);    

    %% Initial Conditions of the G1/S transition model
    P0             = 2e10;                      % Initial protein content
    R0             = G(5)*P0;                   % Initial number of ribosomes
    Cln3_ER_0      = k0*P0;                     % Newborn cells have all Cln3 in the ER
    Cln3Ydj1_c_0   = 0;
    Cln3Ydj1_n_0   = 0;
    Cln3_n_0       = 0;
    Cln3Far1_0     = 0;
    Far1_0         = MA(1);                     % Newborn cells have the total amount of Far1 free of Cln3
    Swi6Swi4_0     = 0;
    Swi6Swi4Whi5_0 = 0;
    Swi6Mbp1_0     = 0;
    Whi5_phospho_0 = [MA(4);zeros(nw1,1)];      % Newborn cells have the total amount of Whi5 not-phosphorylated and free of any binding
    Cln2_0         = 0;
    Nrm1_0         = 0;
    Cln1_0         = 0;
    Clb6_0         = 0;
    Clb5_0         = 0;
    Clb6Sic1_0     = 0;
    Clb5Sic1_0     = 0;
    Clb6Sic1p_0    = 0;
    Clb5Sic1p_0    = 0;
    Sic1_n_0       = MA(6);                     % Newborn cells have the total amount of Sic1 not-phosphorylated and free of any binding
    Sic1p_0        = 0;
    
    Vi = P0/G(7);
    
    Pd_0	= [1;zeros(nd-1,1)];
    Pr_0    = [1;zeros(nr-1,1)];
    X(:,1) = [P0;R0;Cln3_ER_0;Cln3Ydj1_c_0;Cln3Ydj1_n_0;Cln3_n_0;Cln3Far1_0;Far1_0;Swi6Swi4_0;Swi6Swi4Whi5_0;Swi6Mbp1_0;Cln2_0;Nrm1_0;Cln1_0;Clb6_0;Clb5_0;Sic1_n_0;Clb6Sic1_0;Clb5Sic1_0;Clb6Sic1p_0;Clb5Sic1p_0;Sic1p_0;Pd_0;Pr_0;Whi5_phospho_0];
    nX      = length(X);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SIMULATION FROM BIRTH TO BUD                                          %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    i  = 1;
    T(1,1)  = 0;
    dt = 1e-3;      % Integration step. Important: do not set the integration time greater than 1e-3
    State = 0;      % Cell phase: if state = 0 cell is in T1; if state = 1 cell is in T2; if state = 2 end simulation. 
    
    while State ~= 2
        dX = der_plugin(T(1,i), X(:,i), n, N, Dd, Dk, Rb, Rt, Rk, Rw, Fd, Fr, Ud, Ur, Hr, G, CYF, G1S, CLNB, S1, Eta, MA, nX);
        
        i = i + 1;
        T(1,i) = T(1,i-1) + dt;
        X(:,i) = X(:,i-1) + dX*dt;
        
        Pr              = X(22+nd+1:22+nd+nr,i);
        aA              = Fr*Pr;
        Whi5_c          = MA(4) - (X(10,i) + sum(X(22+nd+nr+1:nX,i)) + (N(1) + N(3))*(1-Pr(1)-aA));
        Sic1_c = MA(6) - sum(X(17:22,i));
      
        if Whi5_c > 0.5*MA(4) && State == 0      % Computation of t1 and timer T1 length
            T1 = T(i);     % T1 is over (and T2 starts) when cytoplasmic Whi5 definitely exceeds 50% of Whi5_tot
            State = 1;
            Vs = X(1,i)/G(7);
        elseif Sic1_c > 0.5*MA(6) && State == 1  % Computation of t2 and timer T2 length
            T2 = T(i)-T1;     % T1 is over (and T2 starts) when cytoplasmic Whi5 definitely exceeds 50% of Whi5_tot
            Ps = X(1,i);
            State = 2;
        end
        
        
    end
    
    alpha = (Vs -Vi) /T1;
    
    Pd              = X(22+1:22+nd,:);
    Pr              = X(22+nd+1:22+nd+nr,:);
    aA              = Fr*Pr;
    aB              = Fd*Pd;
    aC              = aA + aB - aA.*aB;    
    
    %% Fraction of G1/S regulon activated genes
    AG	= (N(1)*aA + N(2)*aB + N(3)*aC)/(N(1)+N(2)+N(3));
    
end
