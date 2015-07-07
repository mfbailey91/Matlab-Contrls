%% Bode Method of designing Lead Lag Compensator
clearvars
%   1. Specify the desired cross over frequency in rad/sec
       w_cr = 50
%  2. Specify desired phase margin in degrees: 
      PM_d = 1.92
%  3. Specify DC loop gain
      Kv=50;
      SysType=0;
%  4.  Specify Tau2
    t2=100;
% Open loop transfer function G: 

 num=[10 0]
 den=[1   1  0] 
 Gtf=tf(num,den)
 syms w
 numg = poly2sym(num,j*w); %put as syms so we can subsitute w_cr
 deng= poly2sym(den,j*w);  % ' '
 Go=numg/deng;
 %% Calculate necessary phase angle 
 numg=subs(numg,w, w_cr);
 deng=subs(deng,w, w_cr);
 G=numg/deng;
 Gjwcr=angle(G); vpa(Gjwcr,4);
 Djwcr=PM_d-pi-Gjwcr; vpa(Djwcr,4)
 %% Calculate alpha1
 fprintf('Lead Compensator Parameters \n')
 a1=(1+sin(Djwcr))/(1-sin(Djwcr)); 
 a1=round(a1)
 %% From Given w_cr, calculated a1 calculate t1
 t1=1/(sqrt(a1)*w_cr); tau1=vpa(t1)
 z1=1/(a1*t1); Zero1=vpa(z1,3)  %vpa for display purposes
 p1=1/t1; Pole1=vpa(p1,3)  %vpa for display purposes
 %% Calculate loop gain
 K=1/[sqrt(a1)*abs(G)]; vpa(K,3)
 K=ceil(K)  %round up to next interger
 
%% Lag Compensator
fprintf('Lead Compensator Parameters \n')
if SysType==0;
a2=(Kv/(subs(K*abs(Go),w,0)))^-1
    else if SysType==1;
a2=(Kv/(subs(K*abs(j*w*Go),w,0)))^-1
    else SysType==2;
a2=(Kv/(subs(K*abs((j*w)^2*Go),w,0)))^-1
        end
end
tau2=t2
z2=1/(a2*t2); Zero2=vpa(z2,3)
p2=1/t2;      Pole2=vpa(p2,3)
fprintf('====================================== \n Controller Design \n')
syms s
 controller=K*(a1*t1*s+1)/(t1*s+1)*1/a2*(a2*t2*s+1)/(t2*s+1);
 controllerfactor=factor(controller,s);
 K=controllerfactor(1)
numlead=sym2poly(s+z1);denlead=sym2poly(s+p1); 
lead=tf(numlead,denlead)
numlag=sym2poly(s+z2); denlag=sym2poly(s+p2);
lag=tf(numlag,denlag)
 denf= poly2sym(den,s);  % ' '
 numf= poly2sym(num,s);  % ' '
 num=sym2poly(K*numf*(s+z1)*(s+z2)) ; 
 den=sym2poly(denf*(s+p1)*(s+p2))  ;
H=tf(num,den);


figure(1)
bode(Gtf,':')
hold on
bode(lead,'--')
hold on
bode(lag, '--')
hold on
bode(H)
legend('Open Loop TF', 'Lead', 'Lag' ,'Final Design') 
hold off
figure(2)


subplot(2,1,1)
rlocus(H)

 nums=(K*numf*(s+z1)*(s+z2)) ; 
 dens=(denf*(s+p1)*(s+p2))  ; % ' '
 Gs=nums/dens;
OLTF=Gs
CLTF=simplify(OLTF/(1+OLTF));
factor(CLTF)
[n,d]=numden(CLTF);
n=sym2poly(n);
d=sym2poly(d);
H=tf(n,d)
Zeros=zero(H)
Poles=pole(H)
subplot(2,1,2)
rlocus(H)


