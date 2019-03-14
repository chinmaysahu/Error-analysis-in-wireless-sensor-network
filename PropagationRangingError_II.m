%two anchors and several targets are on same line. The two anchors are used
%to identify 1st target. Then the first anchor is discarded. The second
%anchor and target is used to estimate the next target and so on...
clc
clear all
close all
EstimatedTarget=zeros(15,1000000); % intitializng targets

for trials=1:100000000% number of trials
    trials
    
Targets=[2,5,10:40:600]; % first two targets are anchors included and rest are targets
% count=1;


anchor1=[Targets(1),0];
anchor2=[Targets(2),0];

H=[anchor2-anchor1];

for n=1:(length(Targets)-2) % target number: 1st target starts from 3rd element in Targets array
 rho=EstimateRho(Targets(n+2),n+1,Targets);
 noise=4*randn;  

if(n>1)
H=[H;[[estimatedTarget{n-1}(1),0]-anchor1]];
end

%Estimate target location
H_pseudo = pinv(H); % Use pseudo inverse to compute localization estimate

if(n==1)
b_e = 1/2.*[norm(anchor2).^2-norm(anchor1).^2-(rho(2)^2-rho(1)^2+noise)]; % b estimate matrix with introduced error
else
b_e=EstimateB_b_e_matrix(rho,noise,Targets,n,estimatedTarget);   
end
estimate=H_pseudo*b_e;% compute localization estimate
estimatedTarget{n}=estimate;
EstimatedTarget(n,trials)=estimate(1); % assignng estimate
% EstimatedTarget=[EstimatedTarget;estimatedTarget];
 error_at_target_estimation=(estimatedTarget{n}(1)-Targets(n+2)); % error
end
end

% Variance Calculation
Target=10:40:600;
for n=1:15
   Variance_of_estimates(n)=var(EstimatedTarget(n,:));
end

stem(10:40:600,Variance_of_estimates);
hold on 
title('Variance graph with Noise N(0,1000000) & targets 40 units apart with first target at 10.');


function [rho]=EstimateRho(target,n,Targets)
for j=1:n
   rho(j)=pdist2([Targets(j),0],[target,0],'euclidean');
end
end

% function to update b_e
function [b_e]=EstimateB_b_e_matrix(rho,noise,Targets,n,estimatedTarget)
b_e=[];
for j=1:n    
    noise=1000*randn;
    if j==1
       e=1/2.*[norm(Targets(j+1)).^2-norm(Targets(1)).^2-(rho(j+1)^2-rho(1)^2+noise)]; 
    else
        anchorn=estimatedTarget{j-1};
        e=1/2.*[norm(anchorn).^2-norm(Targets(1)).^2-(rho(j+1)^2-rho(1)^2+noise)]; 
    end % end if else
       b_e=[b_e;e]; % update b_e
end % for loop
end % function
