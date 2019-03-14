%two anchors and several targets are on same line. The two anchors are used
%to identify 1st target. Then the first anchor is discarded. The second
%anchor and target is used to estimate the next target and so on...
clc
clear all
close all
EstimatedTarget=zeros(10,50); % intitializng targets
for trials=1:100000 % number of times experiments is run
    trials  % showing code progress
x1=2;
y1=0;
anchor1=[x1,y1];  %first time anchor1

x2=5;
y2=0;
anchor2=[x2,y2]; %first time anchor2

x0=10;
y0=0;
target=[x0,y0]; %first time target

dist=40; % define the distance between targets 
Error_at_target_estimation=0;
for n=1:15  % number of targets involved at dist=40
    
[estimatedTarget,Error_at_target_estimation(n)] = EstimateTarget(anchor1,anchor2,target,n,Error_at_target_estimation); % estimation call

EstimatedTarget(n,trials)=estimatedTarget(1); % assignng estimate

% plotting target and estimated targets for each new target
% stem(EstimatedTarget(n,trials),zeros(1,1)); 
% hold on
% stem(target(1),0,'filled',':diamondr');

anchor1=anchor2; %anchor2 is new anchor1
anchor2(1)=estimatedTarget(1); % EstimatedTarget is new anchor2 
target(1)=target(1)+dist; % new target spaced 
end

% TargetCollected{trials}=EstimatedTarget;
ErrorCollection{trials}=Error_at_target_estimation;

end
Target=10:40:600;
for n=1:15
   Variance_of_estimates(n)=var(EstimatedTarget(n,:));
end

stem(10:40:600,Variance_of_estimates);
hold on 
title('Variance graph with Noise N(0,16) & targets 40 units apart with first target at 10.');



 function [estimatedTarget,error_at_target_estimation]=EstimateTarget(anchor1,anchor2,Target,n,Error_at_target_estimation)
 % for n=3 and higher anchor1 needs to be orginal anchor2 location not the estimate 
 if n<3
 rho1=pdist2(anchor1,Target,'euclidean');% euclidean distances rho1=|anchor1-target|
 else
 rho1=pdist2(anchor1-[Error_at_target_estimation(n-2),0],Target,'euclidean');% euclidean distances rho1=|anchor1-target|
 end

%for n=2 and higher anchor2 needs to be orginal anchor2 location not the
%estimate
if n<2
rho2=pdist2(anchor2,Target,'euclidean');% euclidean distances rho1=|anchor1-target|
else
rho2=pdist2(anchor2-[Error_at_target_estimation(n-1),0],Target,'euclidean');% euclidean distances rho1=|anchor1-target|
end

%error constant
k=rho2^2-rho1^2;
delk=4*randn;
%Set up matrices
H=[anchor2-anchor1]';% Hermitian matrix

%Estimate target location
H_pseudo = pinv(H); % Use pseudo inverse to compute localization estimate

b_e = 1/2.*[norm(anchor2).^2-norm(anchor1).^2-(k+delk)]; % b estimate matrix with introduced error

estimatedTarget=H_pseudo*b_e;% compute localization estimate

error_at_target_estimation=(estimatedTarget(1)-Target(1)); %squared error
%  end
 
 %find the mse
% mse= mean(error_at_target_estimation);

end

               



