
function [best_fit , meanPar ] = CE_EPSO_opt(bus,Pg,Qg,totalbus,totalbranch,all_busdata14, all_generatordata14,branch,generator)

D=9;

    Pmin=bus.RealPowerMin;
    Pmax=bus.RealPowerMax;
    %Active_PowerGeneration=Pg' % active power output of all generators
    
    Qmin=bus.ReactivePowermin;
    Qmax=bus.ReactivePowermax;
    Reactive_PowerGeneration=Qg'; % active power output of all generators
   % Pgi=Pg';
    

nbus = totalbus; %  total number of buses (14)
bus.VoltageMagnitudeMin = all_busdata14(:,10);%min limit of voltage
bus.VoltageMagnitudeMax = all_busdata14(:,9); %max limit of voltage
VMAX=bus.VoltageMagnitudeMax; %min limit of voltage
VMIN=bus.VoltageMagnitudeMin; %max limit of voltage

Xmin=[Pmin(2,1) Pmin(3,1) Pmin(6,1) Pmin(8,1) 0.9 0.9 0.9 0.9 0.9];
Xmax=[Pmax(2,1) Pmax(3,1) Pmax(6,1) Pmax(8,1) 1.10 1.10 1.10 1.10 1.10];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CE Initializations
n_samples =100;
ro = 0.1;
sigma = 0.85;
beta = 0.9;
q = 5;
indexRo = floor( ( 1 - ro ) * n_samples );
% indexRo = n_samples - 10;
H = zeros( 1, n_samples );
pos = zeros( n_samples, D );
gammaIter = zeros( 1, 1e6 );
x_min = zeros( 1, D );
x_max = zeros( 1, D );
n_trials = ones( 1, D );
meanPar = zeros( 1, D );
stdPar = zeros( 1, D );

for i = 1 : D
     % Continuous variables
        x_min( 1, i ) = Xmin( 1, i );
        x_max( 1, i ) = Xmax( 1, i );
        meanPar( 1, i ) = x_min( 1, i ) + rand() * ( x_max( 1, i ) - x_min( 1, i ) );
        stdPar( 1, i ) = ( x_max( 1, i ) - x_min( 1, i ) ) * 5;
end

iter1 = 0;
iter1_max=30;
memMeanPar( iter1 + 1, : ) = meanPar;
memStdPar( iter1 + 1, : ) = stdPar;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fitness Function Initializations
global ff
ff.coeff = ones( 1, 4 );
ff.weights = zeros( n_samples, 4 );
ff.i_sample = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%
while ( iter1 < iter1_max )
           
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Generate population
        for i = 1 : D
             % Continuous variables + OLTC
                tmp_stat = normrnd( meanPar( i ), stdPar( i ), 1, n_samples );
                tmp_stat( tmp_stat < x_min( 1, i ) ) = x_min( 1, i );
                tmp_stat( tmp_stat > x_max( 1, i ) ) = x_max( 1, i );
                tmp_eval = tmp_stat;
                pos( :, i ) = tmp_eval;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for z=1:n_samples
    
    
     all_generatordata14(2,2)=pos(z,1);
     all_generatordata14(3,2)=pos(z,2);
     all_generatordata14(4,2)=pos(z,3);
     all_generatordata14(5,2)=pos(z,4);
     all_busdata14(1,7)=pos(z,5);
     all_busdata14(2,7)=pos(z,6);
     all_busdata14(3,7)=pos(z,7);
     all_busdata14(6,7)=pos(z,8); 
     all_busdata14(8,7)=pos(z,9);
     
    

 [Ybus, Yf, Yt, Ys] = Ybus_1(totalbus,totalbranch,all_busdata14,branch,bus);
 LFNEWTON;
  Pgggggg(z,:)=Pg;
  Qgggggg(z,:)=Qg;
 
 %penalty value calculation for bus voltage violation
    Voltage_magnitude=Vm; % voltage mag. of each bus
   for count_1=1:nbus
        if Voltage_magnitude(count_1)>VMAX(count_1,:)
            penalty_volt(count_1)=10000*(Voltage_magnitude(count_1)-VMAX(count_1,:))^2;
        elseif Voltage_magnitude(count_1)<VMIN(count_1,:)
            penalty_volt(count_1)=10000*(VMIN(count_1,:)-Voltage_magnitude(count_1))^2;
        else
            penalty_volt(count_1)=0;
        end
   end
   penalty_volt=sum(penalty_volt);%summation of penalty for bus voltage violation
   
   

     Active_PowerGeneration=Pg';
    for count_2=1:nbus
        if  Active_PowerGeneration(count_2)>Pmax(count_2,:)
            penalty_active(count_2)=10000*((Active_PowerGeneration(count_2)-Pmax(count_2,:))^2);
        elseif Active_PowerGeneration(count_2)<Pmin(count_2,:)
            penalty_active(count_2)=10000*((Pmin(count_2,:)-Active_PowerGeneration(count_2))^2);
        else
            penalty_active(count_2)=0;
        end
    end
    penalty_active=sum(penalty_active);  %summation of penalty for Generator Real Power Bounds Violtaion
    
    

  % penalty value calculation for Generator Reactive Power Bound violation 
    Reactive_PowerGeneration=Qg';   % active power output of all generators
    for count_3=1:nbus
        if  Reactive_PowerGeneration(count_3)>Qmax(count_3,:)
            penalty_reactive(count_3)=10000*((Reactive_PowerGeneration(count_3)-Qmax(count_3,:))^2);
        elseif Reactive_PowerGeneration(count_3)<Qmin(count_3,:)
            penalty_reactive(count_3)=10000*((Qmin(count_3,:)-Reactive_PowerGeneration(count_3))^2);
        else
            penalty_reactive(count_3)=0;
        end
    end
    penalty_reactive=sum(penalty_reactive);   %summation of penalty for Generator Real Power Bounds Violtaion


    
     APF=Pij(:,1).^2;
     RPF=Qij(:,1).^2;
     SPF=sqrt(APF+RPF); %actual values of apparent powers "From bus injection" 
     APF1=Pji(:,1).^2;
     RPF1=Qji(:,1).^2;
     SPF1=sqrt(APF1+RPF1); %actual values of apparent powers "To bus injection" 
     jjj=[];
     nbr=totalbranch;

%    %penalty value calculation for line flow violations
  for count_1=1:nbr; % total number of lines (20)
      if SPF(count_1)>branch.PowerMagnitudeMax(count_1,:)
          kkk(count_1)=10000*((SPF(count_1)-branch.PowerMagnitudeMax(count_1,:))^2); % if limit violates,penalty is imposed
      else
          kkk(count_1)=0; %if no violation, penalty is zero
      end
  end
penalty_line=sum(kkk);  % summation of all penalty of line flow violations%



  for count_1=1:nbr; % total number of lines (20)
      if SPF1(count_1)>branch.PowerMagnitudeMax(count_1,:)
          kkk(count_1)=10000*((SPF1(count_1)-branch.PowerMagnitudeMax(count_1,:))^2); % if limit violates,penalty is imposed
      else
          kkk(count_1)=0; %if no violation, penalty is zero
      end
  end
penalty_line1=sum(kkk);


Total_penalty= penalty_volt+ penalty_active + penalty_reactive + penalty_line + penalty_line1;
%Active_PowerGeneration= Active_PowerGeneration';
%% Calculation of Objective Function
 Pgi=Pg';
for p=1:nbus;
    
generator.a1= zeros(totalbus,1);
generator.a1(all_generatordata14(:,1),1) = generator.RealPowerCostCoefficient(:,1);
a1=generator.a1;

generator.b1= zeros(totalbus,1);
generator.b1(all_generatordata14(:,1),1) = generator.RealPowerCostCoefficient(:,2);
b1=generator.b1;

generator.c1= zeros(totalbus,1);
generator.c1(all_generatordata14(:,1),1) = generator.RealPowerCostCoefficient(:,3);
c1=generator.c1;


 F1(p,1)=(Pgi(p,1).*Pgi(p,1))*c1(p,1)+Pgi(p,1)*b1(p,1)+(a1(p,1));

end

% global sum_Fuelcost
% sum_Fuelcost=(sum(F1));
global Total_penalty
fun_Swarm=(sum(F1))+(Total_penalty);
% Objective function= sum of active power losses of the transmission lines
fit(z,:)= sum(F1)+Total_penalty;
end 

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Evaluate population
        %[ fit,~,~,pos,~ ] = feval( fhd,ii,jj,kk,args,pos );
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Performs CE correction
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%our code
        [ Saux, Iaux ] = sort( fit', 2, 'ascend' ); % lowest to highest

        H = zeros( 1, n_samples ); % 1*N
        count=n_samples-indexRo;
        for i=1:count
            H(Iaux(i))=1;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %our code starts
       
    old_mu=meanPar;
    old_std=stdPar;
    pos=pos(Iaux,:); 
    best_fit=Saux(1);
    xx=pos(1:count, :);
    mu=mean(xx); % Consider top best 100 x, columnwise mean of x1 and x2
    
    mu = mu * sigma + old_mu * ( 1 - sigma );
    meanPar=mu;
    sigmaStd1 = beta - beta * ( 1 - 1 / ( iter + 1 ) )^q;
    sigma2=sqrt(var(xx,1)); % same as above column wise var of x1 and x2
    sigma2 = sigma2 * sigmaStd1 + old_std * ( 1 - sigmaStd1 );
    stdPar=sigma2;
    
        
        
        
        % Computes gammaIter
        gamma = Saux( 1, indexRo );
        gammaIter( 1, iter1 + 1 ) = gamma;
        
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Update the weights of the fitness function
        tmpMat = ff.weights( 1:ff.i_sample, : );
        tmp = max( tmpMat );
        ff.coeff = 10 .^ round( log10( max( tmp )./ tmp )  );
        ff.coeff( isinf( ff.coeff ) ) = 1;
        ff.i_sample = 0;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Update CE iteration
        iter1 = iter1 + 1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        memMeanPar( iter1 + 1, : ) = meanPar;
        memStdPar( iter1 + 1, : ) = stdPar;
  
    
    
end
meanPar
best_fit