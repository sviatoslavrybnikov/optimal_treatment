%====================
% Code for: "Integrating gene drives with established pest controls to mitigate spillover risk"
% Authors: Sviatoslav R. Rybnikov, Adam Lampert, Gili Greenbaum
% Description: The code identifies the optimal multi-treatment strategy (gene drives/pesticides/sterile males)
%====================

clear all;
global s c h r st;

%====================
% PARAMETERS - can be modified
%====================
disp('Biological parameters:');
s=0.7%fitness cost
c=0.9%conversion rate
h=0.5%dominance
r=2%intrinsic growth rate

disp('Initial conditions:');
q0=0.01%initial frequency of gene0drove allele
n0=1%initial population density (carrying capacity)

disp('Per-unit treatment costs:');
alpha=0.6%pesticides (relative ti spillover)
beta=3%sterile males (relative ti spillover)

disp('Computational parameters:');
st=1e-4%accuracy
tmax=40%number of generations

%====================
% FUNCTIONS
%====================

%--------------------
% Rounds population density to the nearest discrete state or to zero 
%--------------------
function res_n=myround(arg_n,arg_st)
 if arg_n<arg_st
  res_n=0;
 else
  res_n=arg_st*round(arg_n/arg_st);
 end 
end%myround

%--------------------
% Calculates population density in the next generation based on the Beverton-Holt growth model
%--------------------
function res_n=BH(arg_n,arg_r)
 res_n=arg_r*arg_n/(1+(arg_r-1)*arg_n);
end%BH

%--------------------
% Calculates dynamics of population with gene drives
% Returns:
% (1) DynQ = dynamics of gene-drive allele frequecny
% (2) DunN = dynamics of population density
% (3) DynW = dynamics of mean fitness
%--------------------
function [Res_DynQ,Res_DynN,Res_DynW]=get_Dyn(arg_q0,arg_n0,arg_tmax)
global s c h r;
 q=arg_q0;
 n=arg_n0;
 SC=[1 (1-h*s)*(1-c) (1-h*s)*c 1-s];
 
 QN=[q n];
 for t=1:arg_tmax
  P=[(1-q)^2 2*q*(1-q) 2*q*(1-q) q^2].*SC;%mating+selection+conversion
  q=P*[0 0.5 1 1]'/sum(P);
  n=BH(n,r)*sum(P);
  QN=[QN;q n];
 end%t
 Res_DynQ=QN(:,1);
 Res_DynN=QN(:,2);
 Res_DynW=sum(SC.*[(1-Res_DynQ).^2 2*Res_DynQ.*(1-Res_DynQ) 2*Res_DynQ.*(1-Res_DynQ) Res_DynQ.^2],2);
end%get_Dyn

%--------------------
% Performs dynamic programming
% Returns:
% (1) Xrem_opt = matrix of the optimal removals, depending on current population densioty and generation 
% (2) Ctot = vector of the minimum total costs depending on initial population size
%--------------------
function [Res_Xrem_opt,Res_Ctot]=DP(arg_Qr,arg_W,arg_alpha,arg_beta)
 global r st;
 
 X=[0:st:1]';%possible states
 T=size(arg_Qr,1);%number of generations+1
 Res_Xrem_opt=zeros(numel(X),T);%optimal actions - removals
 U_opt=zeros(numel(X),T);%optimal costs
 Res_Xrem_opt(:,T)=0;
 U_opt(:,T)=X*arg_Qr(T);
 for t=T-1:-1:1%1%time units
  for i=1:numel(X)
   Xnat=BH(X(i),r)*arg_W(t);%natural BH dynamics
   Xnat=myround(Xnat,st);
   XRem_i=[0:st:Xnat]';%possible removals
   U_i=zeros(size(XRem_i));
   Z_i=zeros(size(XRem_i));
   for j=1:numel(XRem_i)
    Xrem=XRem_i(j);
    Xfin=Xnat-Xrem;%final state
    Usp=X(i).*arg_Qr(t);%cost of spillover
    if Xrem==0
     Utr=0;
    else
     Z1=arg_beta>=2*arg_alpha/(k*log(2));%pesticides can be used
     Z2=arg_beta<=2*arg_alpha/(X(i)*log(2));%sterile males can be used
     Z_i(j)=Z1+2*Z2;%code of treatment (1 = pesticide, 2 = sterile males, 3 = both)
     switch Z_i(j)
      case 1%pesticide alone
       A=log2(Xnat/Xfin);
       Utr=arg_alpha*A;%cost of treatment
      case 2%sterile males alone
       M=0.5*X(i)*(Xnat/Xfin-1);
       Utr=arg_beta*M;%cost of treatment
      otherwise%pesticides+sterile males
       k=X(i)*Xnat/Xfin;
       A=log2(k*log(2)*arg_beta/arg_alpha)-1;
       M=arg_alpha/(arg_beta*log(2))-X(i)/2;
       Utr=arg_alpha*A+arg_beta*M;%cost of treatment
       if isnan(Utr)%handling Ctr=NaN (if alpha=beta=0)
        if Xfin==0
         Utr=inf;
        else
         Utr=0;
        end
       end
     end%switch Z
    end%Utr
    Ufur=U_opt(1+round(Xfin/st),t+1);%cost of further treatment
    U_i(j)=Usp+Utr+Ufur;
   end%j
   [U_opt(i,t),min_pos]=min(U_i);
   Res_Xrem_opt(i,t)=Z_i(min_pos)+XRem_i(min_pos);
  end%i
 end%t
 Res_Ctot=U_opt(:,1);
end%DP

%--------------------
% Pestoes dynamics of population density under optimal treatment strategy
% Returns Nopt = (2T+1)-element vector:
% - initial density no (1st element)
% - T pairs of values showing population density in T generations: before and after treatment
% In odd elements except 1 (i.e., 3, 5 ..., 2T+1)
% - integer part codes the used treatmentL 0 = none, 1 = pesticides, 2 = sterile males, 3 = both
% - decimal part codes population density
% Examples:
% - "1.678" = population density of 0.678 reached by pesticide application (1)
% - "2.345" = population density of 0.345 reached by sterile-male release (2)
%--------------------
function Res_Nopt=get_Nopt(arg_Xrem_opt,arg_W,arg_n0)
 global r st;
 
 N=arg_n0;
 Res_Nopt=N;
 for t=1:numel(arg_W)-1
  Nrem=arg_Xrem_opt(1+round(N/st),t);
  N=BH(N,r)*arg_W(t);
  N=myround(N,st);
  Res_Nopt=[Res_Nopt N];%N before treatment
  N=N-mod(Nrem,1);
  Res_Nopt=[Res_Nopt N+fix(Nrem)];
 end%t 
end%get_Nopt

%====================
% MAIN PART - dynamic programming
%====================
tic;
[DynQ,~,DynW]=get_Dyn(q0,1,tmax);
DynQc=1-(1-DynQ).^2;%dynamics of the gene-drive carriers
DynQc(1)=DynQ(1);%generation 0, gene-drive carriers = q0
disp('Identifying the optimal treatment strategy...');
[Xrem_opt,Ctot]=DP(DynQc,DynW,alpha,beta);
Nopt=get_Nopt(Xrem_opt,DynW,n0);%restores population dynamics followomg the optimal solution

%====================
% RESULTS
% Saves the output row containing:
% - 6 parameters (s, c, h, r, alpha, beta)
% - 2T+1 - dynamics of population density
% - T+1 - dynamics of gene-drive carriers
% - total cost of eradication
%====================
Out=[s c h r alpha beta Nopt DynQc' Ctot(1+round(n0/st))];%output
save('-ascii','out.dat','Out');
disp('The strategy found and saved.');
toc;
