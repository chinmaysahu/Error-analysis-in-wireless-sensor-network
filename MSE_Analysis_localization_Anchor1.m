%MSE analysis%
clear
clc
%%
%Declare anchor nodes

%source location
x0=10;
y0=0;
r0=[x0, y0];% source/target location

%Anchor node 1
%x1=5;
y1=0;
%r1=[x1; y1]; % moving anchor node 1

%Anchor node 2
x2=0;
y2=6;
r2=[x2, y2]; %anchor node 2

%node 3
x3=0;
y3=-5; %anchor node 3
r3=[x3, y3]; %anchor node 3
%%
%Get distances with and without noise
rho2=pdist2(r2,r0,'euclidean'); % euclidean distances rho2=|r2-r0|

rho3=pdist2(r3,r0,'euclidean'); % euclidean distances rho3=|r3-r0|



rho1=[];% initializing rho3 vector which will store arays of rho1
SE=[];%initializing Squared Error(SE) vector which will store arays of se
r0_E=[];%initializing target estimates
Estimated_roots=[];
MSE=[];
var=1;

% creates sets of rho1 when x1 moves from 0 to 10 at interval of 0.1
loc = linspace(1,500,100);
for x1=loc %  x1 varies from 0 - 0.0005 at the interval of 0.00005
    r0_E=[];% intiate every iteration to flush out old information of estimates for a given loc
    r1=[x1, y1]; % anchor node 1 taking x3 as i at interval 0.1
    
    rho1_temp=pdist2(r1,r0,'euclidean');% euclidean distances rho3=|r3-r0|
    
    rho1=[rho1 ; rho1_temp]; %collecting all rho1
    %error constant
    k=rho2^2-rho1_temp^2;
    
     %error constant
    k1=rho3^2-rho1_temp^2;

    %write loop here for finding mse when we use randn as errors
     parfor iter=1:100
      
    %introduced error in estimate
     delk=4*randn; %normally distributed error
      

    %introduced error in estimate
    delk1=4*randn; %normally distributed error
    
    %Set up matrices
    H=[r2-r1;r3-r1];% Hermitian matrix for given r3
    
   %Estimate target location
    H_pseudo = pinv(H); % Use pseudo inverse to compute localization estimate
               
    b_e = (1/2)*[norm(r2)^2-norm(r1)^2-(k+delk);
               norm(r3)^2-norm(r1)^2-(k1+delk1)]; % b estimate matrix with introduced error
            
    r0_estimate=H_pseudo*b_e;% compute localization estimate
%     r0_estimate = [r0_estimate(2,1); r0_estimate(1,1)];
    r0_E=[r0_E;r0_estimate'];
    
    diff=sum((r0_estimate'-r0).^2);% finding |r0-r0_estimate|
    
    SE=[SE;diff];
    end
     % find avg of the roots
    avg_root_loc=mean(r0_E); % avg root location
    Estimated_roots=[Estimated_roots;avg_root_loc];
    
    % find the mse
    mse= mean(SE);
    MSE=[MSE;mse];% update the MSE
end

%mse_vector=(1/length(SE))*SE; % normalizing weight by divinding squared error(SE) vector with total no. of samples  

% mse_total=mean(SE);% mse by finding avg.

figure
grid on
box on
grid minor
hold on
% xlim([0,5e-4]);
plot(loc,MSE,'LineWidth',2);
title('MSE','FontName','Times New Roman','FontSize',12,'FontWeight','bold');
xlabel('x1','FontName','Times New Roman','FontSize',12,'FontWeight','bold');
ylabel('mse','FontName','Times New Roman','FontSize',12,'FontWeight','bold');
set(gca,'FontName','Times New Roman','FontSize',12,'FontWeight','bold','linewidth',2);







