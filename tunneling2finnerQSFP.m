%%tunneling2
%T=transfers1;
clear all
clc;
q=-1.60217662*1e-19;% coulombs;%
bins=1;
I_final=[];
cu=16
for Vg=0:1:cu
 Vg
prefix=num2str(Vg);
M = load(['M' prefix '.mat'], '-mat');
eval(['M=M.M' num2str(Vg)]);


M1=M(:,1:18);
M2=M(:,19:36);
M3=M(:,37:54);
%%finding the tunneling paths.
kk=9;
%%%%%
Icuts=[];
for sl=1:1*kk%%%%%%%%%%%%%%%%%%
% CBT=M1(:,2*sl);
% EFT=M2(:,2*sl)*1e2;%V/m
% QFPT=M3(:,2*sl);%
% xt=M1(:,1).*1e-6;
CBT=(finner(M1(:,2*sl),bins))';
 EFT=((finner(M2(:,2*sl),bins)).*1e2)';%V/m
QFPT=(finner(M3(:,2*sl),bins))';%
 xt=((finner(M2(:,1),bins)).*1e-6)';
    


phi_b=max(CBT);%TOB
index_b=find(CBT==phi_b);
index_bTOB= find(CBT<phi_b);
A=index_bTOB(1:index_b-1,1);%primera mitad de la banda de conduccion
B=index_bTOB(index_b:end,1);%2da mitad de la banda de conduccion
%padding

[mb, nb]=size(B);%
[ma, na]=size(A);%
delta=0.005;
newindex=zeros;
delta2=zeros;
iA=[];
jA=[];
%finding the tunneling paths
for i=1:ma
    %i
     for j=1:mb
        % j
    delta3=abs((CBT(A(i,1),:)-CBT(B(j,1),:)));
    delta2=[delta2 delta3];
        if delta3< delta
        iA=[i iA];
        jA=[j jA];
        newindex=[newindex, B(j)];
        else
        end
  
    end

end
figure(33)
 plot(xt,CBT)
 hold on
 plot(xt(iA),CBT(iA),'-s')
 plot(xt(B(jA)),CBT(B(jA)),'-s')
 xlabel('x(nm)')
 ylabel('CB (eV)')


% %%% parameters
 hbar=1.054571800*1e-34; %Js o kg*m^2*s^-1;
 kev=1.6021766208*1e-19; %(kgm^2)/s^2
 me=9.10938356*1e-31; %Kg
 mt=me*0.0516;
 T=300;%k;
 kb=8.6173303*1e-5; %eV.K-1
 A=1.2*1e-6;%Am^-2K^-2
 q=1.60217662*1e-19;
 
 Tpaths=[iA' B(jA)];
 [mp, np]= size(Tpaths);
% get dE for each Tpaths


%for each tunneling path 
dxt=zeros;
diff1=zeros;
dE=zeros;
diff=zeros;    
dx1=zeros;
T0=zeros;
Tg=zeros;
gamma=zeros;
djcc=zeros;
Gp=[];
gamma3=[]
dxp=[];
E_Field=zeros;
v1=zeros;
v2=zeros;

 
 for itp=1:mp
  %all the quantities refered to the tunneling path  
 xp=xt(Tpaths(itp,1):Tpaths(itp,2),:);
 CBp=CBT(Tpaths(itp,1):Tpaths(itp,2),:);
 %QSFP=QFPT(Tpaths(itp,1):Tpaths(itp,2),:);
 QSFP=maxQSF(QFPT);
 EF=EFT(Tpaths(itp,1):Tpaths(itp,2),:);

 %discretization WKB
 gamma1=[];
 lx=length(xp);    
     for ixt=1:lx-1
      diff(ixt,1)=phi_b-CBp(ixt,1);
      dx1(ixt,1)=xp(ixt+1,:)-xp(ixt,:);
      T0(ixt,1)=dx1(ixt,:)'.*(sqrt(2*mt*kev*(diff(ixt,:))));
      gc=2;
      Acc=A*gc;
      Tg=sum(T0);
      tgsize=size(Tg);
      factor=2/hbar;
      gamma=exp(-factor*Tg);
     
     end
       gamma1=[gamma1 gamma];
   gamma2(mp)=gamma;   
  argc2=EF(1,:);
  if argc2>0
       theta_c2=1;
  else
       theta_c2=0;
  end
  
  argc1=-EF(end,:);
  if argc1>0
       theta_c1=1;
  else
      theta_c1=0;
  end
 %epsilon is a constant energy, of the tunneling path
 %can be the average between both extremes of the tunneling
 %path
 epsilon=(CBp(1)+CBp(end))/2;%eV
 v1(1,:)=(epsilon-CBp(end))*(abs(-EF(end)))*theta_c1;
 v2(1,:)=(epsilon-CBp(1))*(abs(EF(1)))*theta_c2;
 gamma3=[gamma3 gamma];
 
 G=+real((Acc*T)*(1/(q*kb))*v1*v2*gamma*(+(1+log((QSFP(end,:)-epsilon)/(kb*T)))));
 R=real((Acc*T)*(1/(q*kb))*v1*v2*gamma*(-(1+log((QSFP(1,1)-epsilon)/(kb*T)))));
 %find dx for the end of the tunneling path
 dx=xp(end,:)-xp(end-1,:);
 %cleaning the discretization
 dx1=zeros;
 diff=zeros;
 %gamma=0;
 v1=zeros;
 v2=zeros;
 QSFP=zeros;
 Gp=[Gp G];
 dxp=[dxp dx];
 end
 
dE=1;%0.008
%no se como  
%dE(itp)=CBp(1,:)-CBp(end,:);
I=-q*sum(Gp.*dxp); % the current density
I_cuts=I.*((4.6/kk).*1e-12); % divide by one micrometer
Icuts=[Icuts I_cuts];
%%cleaning
I=zeros;
Gp=zeros;
dxp=zeros;
CBT=zeros;
EFT=zeros;%V/m
QFPT=zeros;%
xt=zeros;

end
I_x=(sum(Icuts))/kk;

%%cleaning the cut matrices for each VG-point
M=zeros;

I_final=[I_final; I_x]
end

save('perro1bins.mat','I_final')
% my_tunneling=transfers1(:,2)+I_final;
% 
% figure (66)
% hold on
% plot(transfers1(:,1),my_tunneling)
% plot(transfers1(:,1),transfers1(:,4))
% plot(transfers1(:,1),transfers1(:,2))
% hold off
%I_final=[I_final I_x]
%my_tunneling=x+I_tcadnone;
 
 
