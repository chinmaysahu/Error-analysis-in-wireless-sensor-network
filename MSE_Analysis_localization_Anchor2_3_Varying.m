%MSE analysis%
clear
clc
%%
%Declare anchor nodes

%source/target location
x0=10;
y0=0;
r0=[x0, y0];% source/target location

% Anchor node 1
x1=5;
y1=0;
r1=[x1, y1]; % moving anchor node 1

%node 2
x2=0;
%y2=5;
%r2=[x2; y2]'; %anchor node 2

%Anchor node 3
x3=0;
%y3=-5; %anchor node 3
%r3=[x3; y3]'; %anchor node 3
%%
%Get distances with and without noise
rho1=pdist2(r1,r0,'euclidean');% euclidean distances rho1=|r1-r0|


rho2=[];% initializing rho3 vector which will store arays of rho2
rho3=[];% initializing rho3 vector which will store arays of rho2
SE=[];%initializing Squared Error(SE) vector which will store arays of SE
r0_E=[];%initializing target estimates
MSE=[];


% creates sets of rho2,rho3 when y2,y3 moves from 0 to 10/-10 at interval of 0.1
 loc = linspace(1,10,1000);
for y2=loc
    r2=[x2, y2]; % anchor node 2 taking y2 as i at interval 0.1
    r3=[x3, -y2]; % anchor node 3 taking y3 as i at interval 0.1
    
    rho2_temp=pdist2(r2,r0,'euclidean'); % euclidean distances rho2=|r2-r0|
    rho3_temp=pdist2(r3,r0,'euclidean'); % euclidean distances rho3=|r3-r0|
    
    %error constant
    k=rho2_temp^2-rho1^2;
    %error constant
    k1=rho3_temp^2-rho1^2;
    %rho2 array
    rho2=[rho2 ; rho2_temp]; %collecting all rho2
    %rho3 array
    rho3=[rho3 ; rho3_temp]; %collecting all rho3
    
    %write loop here for finding mse when we use randn as errors
    for iter=1:1:100
    %introduced error in estimate
    delk=randn; %normally distributed error
    %introduced error in estimate
    delk1=randn; %normally distributed error
    %Set up matrices
    H=[r2-r1;r3-r1];% Hermitian matrix for given r3
    
    %Estimate target location
    H_pseudo = pinv(H); % Use pseudo inverse to compute localization estimate
               
    b_e = 1/2.*[norm(r2).^2-norm(r1).^2-(k+delk);
               norm(r3).^2-norm(r1).^2-(k1+delk1)]; % b estimate matrix with introduced error
            
    r0_estimate=H_pseudo*b_e;% compute localization estimate
%     r0_estimate=flip(r0_estimate);
%     r0_E=[r0_E;r0_estimate'];% location estimates
    
    diff=norm(r0_estimate'-r0).^2;% finding |r0-r0_estimate|
    
    SE=[SE;diff];
    end
    % find the mse
    mse= mean(SE);
    MSE=[MSE;mse];% update the MSE
end

%mse_vector=(1/length(SE))*SE; % normalizing weight by divinding squared error(SE) vector with total no. of samples  

mse_total=mean(SE);% mse by finding avg.

figure
subplot(1,2,2)
grid on
box on
grid minor
hold on
plot(loc,MSE);
title('MSE','FontName','Times New Roman','FontSize',12,'FontWeight','bold');
xlabel('y2','FontName','Times New Roman','FontSize',12,'FontWeight','bold');
ylabel('mse','FontName','Times New Roman','FontSize',12,'FontWeight','bold');
set(gca,'FontName','Times New Roman','FontSize',12,'FontWeight','bold');

hold on

subplot(1,2,1)
grid on
box on
grid minor
hold on
plot(-loc,MSE);
title('MSE','FontName','Times New Roman','FontSize',12,'FontWeight','bold');
xlabel('y3','FontName','Times New Roman','FontSize',12,'FontWeight','bold');
ylabel('mse','FontName','Times New Roman','FontSize',12,'FontWeight','bold');
set(gca,'FontName','Times New Roman','FontSize',12,'FontWeight','bold');






