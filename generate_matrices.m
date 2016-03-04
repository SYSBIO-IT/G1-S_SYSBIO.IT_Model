function [D_alfa_d,D_k,Fd,Ud,Rb_alfa_t,Rb_alfa_b,Rb_k,Rb_alfa_W,Fbr,Ubr,Hb,nrr]...
          =generate_matrices(ns,nw,no,cell_type)

% Function generate_matrices(ns,nw,nr,cell_type)

% The inputs are the phosphorylations of Swi6 and Whi5
% ns: of Swi6,  nw: useful of Whi5, no: useless of Whi5
% cell_type: 'wt' (wild type), 'sd' (semi-direct) mutato 4e

nspd=ns+2; nwpd=nw+2; nopd=no+2;

% Construction of the matrices of dimer dynamics
% D_alfa_d,D_k,Fd,Ud,

% Dimension of vector d
nd=ns+2;

D_alfa_d=zeros(nd,nd); D_k=zeros(nd,nd);
Fd=zeros(1,nd); Fd(nd)=1;
Ud=zeros(1,nd);

D_alfa_d(1,1)=-1; D_alfa_d(2,1)=1; 

for i=-1:ns
   hd=2+i;
   % Vector Ud (scanned by columns)
   Ud(hd)=(ns-i)*(i>-1);
   % Matrix D_k scanned by rows
   D_k(hd,hd)=-(ns-i)*(i>-1);
   if (i>0), D_k(hd,hd-1)=ns-i+1; end 
end


% Construction of the matrices of trimer dynamics
% R_alfa_t,R_alfa_b,R_k,   R_alfa_W,  Fr,Ur,H

% Dimension of the (not reduced) vector r
nr=nspd*nwpd*nopd;

Znr=zeros(nr,nr);  
R_alfa_t=Znr; R_alfa_b=Znr; R_k=Znr; R_alfa_W=Znr;
Fr=zeros(1,nr); Ur=zeros(1,nr);
Fr(2:nspd)=ones(1,ns+1);  % equazione 10
H=zeros(nw+1,nr);

% Dimension of the reduced vector rr
nrr=(ns+1)*(nw+1)*(no+1)+ns+2;

% Selection matrix
Slz=zeros(nrr,nr);

jrr=0;

for j1=-1:ns
 for j2=-1:nw  % in k5 era nw-1
   for j3=-1:no
    if scheck(j1,j2,j3) % unused stats for scheck(j1,j2,j3)=0
       hr=2+j1+nspd*(j2+1)+nspd*nwpd*(j3+1);
       jrr=jrr+1; % numero di righe della matrice di selezione
       Slz(jrr,hr)=1; % elemento della matrice di selezione
       % Vector Ur (scanned by columns)
       Ur(hr)=((ns-j1)*(j1>-1)+(nw-j2)*(j2>-1)+(no-j3)*(j3>-1));
       % Matrices R (scanned by rows)
       R_k(hr,hr)=-Ur(hr);                                               %1
       if (j1==-1)&&(j2==-1)&&(j3==-1), 
                R_alfa_t(hr,hr)=-1;  R_alfa_b(hr,hr)=-1; end             %2
       if strcmp(cell_type,'wt')
           if (j1==0)&&(j2==0)&&(j3==0),   R_alfa_t(hr,1)=1; end        %3a
       end
       if strcmp(cell_type,'sd')
           if (j1==0)&&(j2==nw)&&(j3==0),   R_alfa_t(hr,1)=1; end       %3b
       end
       if (j1==0)&&(j2==-1)&&(j3==-1), R_alfa_b(hr,1)=1; end             %4
       if (j1>0)&&scheck(j1-1,j2,j3), R_k(hr,hr-1)=ns-j1+1; end          %5
       if (j2>0)&&scheck(j1,j2-1,j3), R_k(hr,hr-nspd)=nw-j2+1; end       %6
       if (j3>0)&&scheck(j1,j2,j3-1), R_k(hr,hr-nspd*nwpd)=no-j3+1; end  %7
       if (j1>-1)&&(j1<ns)&&(j2==-1)&&(j3==-1),
         for j=0:no, 
             jnext=2+j1+nspd*(nw+1)+nspd*nwpd*(j+1);
             R_alfa_W(hr,jnext)=1;  R_alfa_W(jnext,jnext)=-1;            %8
         end
       end
       if (j1==ns)&&(j2==-1)&&(j3==-1),
           for ja=0:nw, 
              for jb=0:no,
                jnext=2+ns+nspd*(ja+1)+nspd*nwpd*(jb+1);
                R_alfa_W(hr,jnext)=1;  R_alfa_W(jnext,jnext)=-1;         %9
              end
           end
       end
     end
   end
 end
end

% Matrix H
for j2=0:nw-1
    for j3=0:no
        H(j2+1,2+ns+nspd*(j2+1)+nspd*nwpd*(j3+1))=1;
    end  
end

% j2=nw;
for j1=0:ns
    for j3=0:no
        H(nw+1,2+j1+nspd*(nw+1)+nspd*nwpd*(j3+1))=1; 
    end
end

% Check dimensioni matrice di Selezione
if not(nrr==jrr), 
    disp('ERROR!!!')
end

% Reduced matrices
Rb_alfa_t=Slz*R_alfa_t*Slz';
Rb_alfa_b=Slz*R_alfa_b*Slz';
Rb_k=Slz*R_k*Slz';
Rb_alfa_W=Slz*R_alfa_W*Slz';   % new entry nel k6
Fbr=Fr*Slz';
Ubr=Ur*Slz';
Hb=H*Slz';


function sout=scheck(j1,j2,j3)
  imp1=(j1==-1)&&((j2>-1)||(j3>-1)); % impossible 1
  imp2=((j2==-1)&&(j3>-1))||((j2>-1)&&(j3==-1)); % impossible 2
  sout=not(imp1||imp2);
end

end
