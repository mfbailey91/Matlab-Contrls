%% Root Locus Lead Lag Compensator Design Coordinator
format compact
clearvars
%% Open loop transfer function
num=[1 ];
deng=[1 6 8];
G=tf(num,deng);
%% Desired Performance Specifications
%transient
pd=[-5+5*1i; -5-5*1i];
%steady state
Kv=50;
%System Type
SysType=0;
%% Lagpole
p2=.01;
%% Lead Design
sd=pd(1);
z1=-real(pd(1))
gs= evalfr(G,pd(1));
PhG= atan2(imag(gs),real(gs));
PhZ= atan2(imag(sd+z1),real(sd+z1) );
Phase1=pi+PhG+PhZ ;
p1=(imag(sd)/tan(Phase1))-real(sd);
p1=round(p1)  %Calculated Pole
%% Solve Gain K
K=abs(-1/([(sd+z1)/(sd+p1)]*gs))
K=round(K)
m=size(num);
%% Compute DC Gain DC gain u lim s->0 sD(s)G(s)
numg=num
if SysType==0
numg=num;
    else if SysType==1
        numg(m+1)=0;  %add s term
    else SysType==2
        numg(m+1)=0; numg(m+2)=0;
        end
end

syms s
numg = poly2sym(numg,s);
deng= poly2sym(deng,s);
G=simplify(numg/deng);
Gv=subs(G,s,0);   %s*G(s) evaluated at s=0
Kvlead=K*(z1)/(p1)*Gv  %Kv gain for lead compensator only

%% Lag Design
p2
z2=Kv/Kvlead*p2
%% Lead Lag Compensator
num = poly2sym(num,s);
G=simplify(num/deng)
syms s
Comp=K*[(s+z1)/(s+p1)]*[(s+z2)/(s+p2)];
fprintf('Closed Loop Transfer Function \n')
OLTF=G*Comp;
CLTF=simplify(OLTF/(1+OLTF))
DCgain=subs(CLTF,s,0)
Kv=subs(Comp,s,0)*Gv
[n,d]=numden(CLTF);
n=sym2poly(n);
d=sym2poly(d);
den=deng*(s+p1)*(s+p2); den=sym2poly(den);
num=1+K*num*(s+z1)*(s+z2); num=sym2poly(num);
H=tf(num, den);
figure
subplot(2,1,1)
rlocus(H)
title('open loop transfer function')

H=tf(n,d)
Zeros=zero(H)
Poles=pole(H)
subplot(2,1,2)
rlocus(H)
title('closed loop transfer function')
        

