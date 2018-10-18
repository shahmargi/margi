


ns=0; %number of slack bus
ng=0; %number of generator bus
Vm=0; %voltage magnitude
delta=0; 
yload=0; 
deltad=0; %phase angle in degree
basemva=100; %assume
accuracy=0.001; % tolerance
accel=1.8; %alpha acceleration factor
maxiter=100; %maximum iteration number
nbus = totalbus; %  total number of buses (14)
nbr=totalbranch; % total number of lines (20)
nl= branch.Origin  ; %from bus 
nr= branch.Destination  ; % to bus 


bus.RealPowerGeneration= zeros(totalbus,1);
bus.RealPowerGeneration (all_generatordata14(:,1),1) = all_generatordata14(:,2); % bus real power generation (14*1)

bus.ReactivePowerGeneration= zeros(totalbus,1);
bus.ReactivePowerGeneration (all_generatordata14(:,1),1) = all_generatordata14(:,3); % bus reactive power generation (14*1)

bus.ReactivePowermin= zeros(totalbus,1);
bus.ReactivePowermin (all_generatordata14(:,1),1) = all_generatordata14(:,5); % reactive power minimum limit (14*1)

bus.ReactivePowermax= zeros(totalbus,1);
bus.ReactivePowermax (all_generatordata14(:,1),1) = all_generatordata14(:,4); % reactive power maximum limit (14*1)


bus.RealPowerMin= zeros(totalbus,1);
bus.RealPowerMin (all_generatordata14(:,1),1) = all_generatordata14(:,17);


bus.RealPowerMax= zeros(totalbus,1);
bus.RealPowerMax (all_generatordata14(:,1),1) = all_generatordata14(:,16);



%length(busdata(:,1));
for k=1:nbus
n=all_busdata14(k,1);                     % bus number
kb(n)=all_busdata14(k,3);                 % types of buses
Vm(n)=all_busdata14(k,7);                 % voltage magnitude of each bus
delta(n)=all_busdata14(k, 8);             % phase angle of each bus
Pd(n)=bus.RealPowerDemand(k,1);           % real power demand of each bus
Qd(n)=bus.ReactivePowerDemand(k,1);       % reactive power demand of each bus
Pg(n)=bus.RealPowerGeneration(k,1);       % real power generation of each bus
Qg(n) = bus.ReactivePowerGeneration(k,1); % reactive power generation of each bus
Qmin(n)=bus.ReactivePowermin(k,1);        % reactive power maximum limit (14*1)
Qmax(n)=bus.ReactivePowermax(k,1);        % reactive power minimum limit (14*1)
Qsh(n)=branch.ChargingSusceptance(k,1);   % shunt Charging Susceptance of each bus (14*1)
Pmin(n)=bus.RealPowerMin(k,1);
Pmax(n)=bus.RealPowerMax(k,1);

    if Vm(n) <= 0  Vm(n) = 1.0; V(n) = 1 + j*0;
    else delta(n) = pi/180*delta(n);
         V(n) = Vm(n)*(cos(delta(n)) + j*sin(delta(n)));
         P(n)=(Pg(n)-Pd(n))/basemva;         % active power at each bus
         Q(n)=(Qg(n)-Qd(n)+ Qsh(n))/basemva; % reactive power at each bus
         S(n) = P(n) + j*Q(n);               % apparent power at each bus
    end
    
end

for k=1:nbus
if kb(k) == 3, ns = ns+1; 
else 
end
if kb(k) == 2 ng = ng+1; 
else
end
ngs(k) = ng;
nss(k) = ns;
end

Ym=abs(Ybus);
t = angle(Ybus); 
m=2*nbus-ng-2*ns; %size of jacobian
maxerror = 1; 
converge=1;
iter = 0; %iteration number

% Start of iterations
clear A  DC   J  DX
while maxerror >= accuracy & iter <= maxiter % Test for max. power mismatch
for i=1:m
for k=1:m
   A(i,k)=0;      %Initializing Jacobian matrix
end, 
end
size(A);
%pause
iter = iter+1;
for n=1:nbus
nn=n-nss(n);
lm=nbus+n-ngs(n)-nss(n)-ns;
J11=0; J22=0; J33=0; J44=0;
for i=1:nbr
     if nl(i) == n | nr(i) == n
        if nl(i) == n,  l = nr(i); 
        end
        if nr(i) == n,  l = nl(i); 
        end
        J11=J11+ Vm(n)*Vm(l)*Ym(n,l)*sin(t(n,l)- delta(n) + delta(l));
        J33=J33+ Vm(n)*Vm(l)*Ym(n,l)*cos(t(n,l)- delta(n) + delta(l));
     if kb(n)~=3
        J22=J22+ Vm(l)*Ym(n,l)*cos(t(n,l)- delta(n) + delta(l));
        J44=J44+ Vm(l)*Ym(n,l)*sin(t(n,l)- delta(n) + delta(l));
     else
     end
        if kb(n) ~= 3  & kb(l) ~=3
        lk = nbus+l-ngs(l)-nss(l)-ns;
        ll = l -nss(l);
        
      % off diagonalelements of J1
        A(nn, ll) =-Vm(n)*Vm(l)*Ym(n,l)*sin(t(n,l)- delta(n) + delta(l));
        
              if kb(l) == 1  % off diagonal elements of J2
              A(nn, lk) =Vm(n)*Ym(n,l)*cos(t(n,l)- delta(n) + delta(l));
              end
              
              if kb(n) == 1  % off diagonal elements of J3
              A(lm, ll) =-Vm(n)*Vm(l)*Ym(n,l)*cos(t(n,l)- delta(n)+delta(l)); 
              end
              
              if kb(n) == 1 & kb(l) == 1  % off diagonal elements of  J4
              A(lm, lk) =-Vm(n)*Ym(n,l)*sin(t(n,l)- delta(n) + delta(l));
              end
              
        else
        end
     else  
    end
end
 
   Pk = Vm(n)^2*Ym(n,n)*cos(t(n,n))+J33; %active power at each bus
   Qk = -Vm(n)^2*Ym(n,n)*sin(t(n,n))-J11; %reactive power at each bus
   
   if kb(n) == 3 P(n)=Pk; Q(n) = Qk; 
   end   % Swing bus P
     if kb(n) == 2  Q(n)=Qk;
%          if Qmax(n) ~= 0
%            Qgc = Q(n)*basemva + Qd(n) - Qsh(n);
%            if iter <= 7                  % Between the 2th & 6th iterations
%               if iter > 2                % the Mvar of generator buses are
%                 if Qgc  < Qmin(n),       % tested. If not within limits Vm(n)
%                 Vm(n) = Vm(n) + 0.01;    % is changed in steps of 0.01 pu to
%                 elseif Qgc  > Qmax(n),   % bring the generator Mvar within
%                 Vm(n) = Vm(n) - 0.01;    % the specified limits.
%                 end 
%               else 
%               end
%            else
%            end
%          else
%          end
     end
   if kb(n) ~= 3
     A(nn,nn) = J11;  %diagonal elements of J1
     DC(nn) = P(n)-Pk;
   end
   if kb(n) == 1
     A(nn,lm) = 2*Vm(n)*Ym(n,n)*cos(t(n,n))+J22;  %diagonal elements of J2
     
     A(lm,nn)= J33;        %diagonal elements of J3
     
     A(lm,lm) =-2*Vm(n)*Ym(n,n)*sin(t(n,n))-J44;  %diagonal of elements of J4
     
     DC(lm) = Q(n)-Qk;
   end
end

DX=A\DC';

for n=1:nbus
  nn=n-nss(n);
  lm=nbus+n-ngs(n)-nss(n)-ns;
    if kb(n) ~= 3
    delta(n) = delta(n)+DX(nn); 
    end
    if kb(n) == 1
    Vm(n)=Vm(n)+DX(lm); 
    end
 end
  maxerror=max(abs(DC));
     if iter == maxiter & maxerror > accuracy 
   fprintf('\nWARNING: Iterative solution did not converged after ')
   fprintf('%g', iter), fprintf(' iterations.\n\n')
   fprintf('Press Enter to terminate the iterations and print the results \n')
   converge = 0; pause, else, end
   
end

if converge ~= 1
   tech= ('                      ITERATIVE SOLUTION DID NOT CONVERGE'); else, 
   tech=('                   Power Flow Solution by Newton-Raphson Method');
end  

V = Vm.*cos(delta)+j*Vm.*sin(delta);
deltad=180/pi*delta;
i=sqrt(-1);
k=0;

for n = 1:nbus
  if kb(n) == 3
      
     k=k+1;
     S(n)= P(n)+j*Q(n);
     Pg(n) = P(n)*basemva + Pd(n);
     Qg(n) = Q(n)*basemva + Qd(n) - Qsh(n);
     Pgg(k)=Pg(n);
     Qgg(k)=Qg(n);  
     
     elseif  kb(n) ==2
     k=k+1;
     S(n)=P(n)+j*Q(n);
     Qg(n) = Q(n)*basemva + Qd(n) - Qsh(n);
     Pgg(k)=Pg(n);
     Qgg(k)=Qg(n);  
     
  end
yload(n) = (Pd(n)- j*Qd(n)+j*Qsh(n))/(basemva*Vm(n)^2);
end

busdata(:,3)=Vm'; 
busdata(:,4)=deltad';
Pgt = sum(Pg);  %total active power generation of the system
Qgt = sum(Qg);  %total reactive power generation of the system
Pdt = sum(Pd) ; %total active power demand of the system
Qdt = sum(Qd) ; %total reactive power demand of the system
Qsht = sum(Qsh); %total shunt reactive power generation of the system




%%%%%%%%%  This program is for the computation of line flow and line losses.

a=branch.TapRatio; %transformer tap ratio 
SLT = 0;
% fprintf('\n')
% fprintf('                           Line Flow and Losses \n\n')
% fprintf('     --Line--  Power at bus & line flow    --Line loss--  Transformer\n')
% fprintf('     from  to    MW      Mvar     MVA       MW      Mvar      tap\n')

for n = 1:nbus
busprt = 0;
Bc=branch.ChargingSusceptance;
   for L = 1:nbr;
%        if busprt == 0
%        fprintf('   \n'), fprintf('%6g', n), fprintf('      %9.3f', P(n)*basemva)
%        fprintf('%9.3f', Q(n)*basemva), fprintf('%9.3f\n', abs(S(n)*basemva))
% 
%        busprt = 1;
%        else 
%        end
       
  if nl(L)==n      
       k = nr(L);
       In = (V(n) - a(L)*V(k))*Ys(L)/a(L)^2; % nth bus current     %(V(n) - a(L)*V(k))*Ys(L)/a(L)^2 + Bc(L)/a(L)^2*V(n);
       Ik = (V(k) - V(n)/a(L))*Ys(L) ;       % kth bus current     %(V(k) - V(n)/a(L))*Ys(L) + Bc(L)*V(k);
       Snk = V(n)*conj(In)*basemva; % apparent power flow from bus n to k
       PF(L,:)=real(Snk);
       QF(L,:)=imag(Snk);
       Skn = V(k)*conj(Ik)*basemva;   % apparent power flow from bus k to n
       PF1(L,:)=real(Skn);
       QF1(L,:)=imag(Skn);
       SL  = Snk + Skn;               % loss of apparent power 
       SLT = SLT + SL;                % total apparent power loss
%        PF=real(Snk)
%        QF=imag(Snk)
       elseif nr(L)==n  k = nl(L); 
           
       In = (V(n) - V(k)/a(L))*Ys(L);         % nth bus current            %(V(n) - V(k)/a(L))*Ys(L) + Bc(L)*V(n);
       Ik = (V(k) - a(L)*V(n))*Ys(L)/a(L)^2 ; % kth bus current            %(V(k) - a(L)*V(n))*Ys(L)/a(L)^2 + Bc(L)/a(L)^2*V(k);
       Snk = V(n)*conj(In)*basemva; % apparent power flow from bus n to k
       PF(L,:)=real(Snk);
       QF(L,:)=imag(Snk);
       Skn = V(k)*conj(Ik)*basemva;           % apparent power flow from bus k to n
       PF1(L,:)=real(Skn);
       QF1(L,:)=imag(Skn);
       SL  = Snk + Skn;                       % loss of apparent power
       SLT = SLT + SL;                        % total apparent power loss
%        PF=real(Snk)
%        QF=imag(Snk)
    else
  end
    
%      if nl(L)==n | nr(L)==n
%          fprintf('%12g', k),
%          fprintf('%9.3f', PF), fprintf('%9.3f', QF)
%          %fprintf('%9.3f', real(Snk)), fprintf('%9.3f', imag(Snk))
%          fprintf('%9.3f', abs(Snk)),
%          fprintf('%9.3f', real(SL)),
%              if nl(L) ==n & a(L) ~= 1
%              fprintf('%9.3f', imag(SL)), fprintf('%9.3f\n', a(L))
%              else fprintf('%9.3f\n', imag(SL))
%              end
%          else 
%       end
   end
end
  Pij=PF;
  Qij=QF;
  Pji=PF1;
  Qji=QF1;
  PLoss=Pij+Pji;
  QLoss=Qij+Qji;
%   Qloss=QL
% SLT = SLT/2;
% fprintf('   \n'), fprintf('    Total loss                         ');
% fprintf('%9.3f', real(SLT)), fprintf('%9.3f\n', imag(SLT));
% clear Ik In SL SLT Skn Snk
   
