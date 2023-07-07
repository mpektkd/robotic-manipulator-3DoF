clear all;  %clear all the variables
close all;  %close all the previous plots
%define the varibles for indexing 
syms x y z ; 
x=1; y=2; z=3;
%define the coordinates
h = 30; zA = h; zB = h; xA = -7; xB = 5; yA = 7; yB = -10;
%define the links(cm)
l2 = 20; l4 = 15; l5 = 23;
%define the parameters
D = 0.5;      %2D = time interval of each of 5th order polynomials
half_tf = 3;    %time of each part oh the motion
tf = 2*half_tf; %final time
dt = 0.001; %sampling rate
t = 0:dt:(tf);  %time sampling
P = zeros(3, tf/dt + 1);    %initialize position
V = zeros(3, tf/dt + 1);    %initialize velocity
A = zeros(3, tf/dt + 1);    %initialize acceleration
coorA = [xA yA];    
coorB = [xB yB];
fig_num = 1;    %variable counting figures
%Task B.7a
P_coeff = {{},{},{},{},{},{}};  %coefficient initialization for each physical quantity
V_coeff = {{},{},{},{},{},{}};
A_coeff = {{},{},{},{},{},{}};
for i = [x y]
    P_coeff{i} = parameters(coorA(i), coorB(i), D, half_tf);   
    %P_coeff is a cell matrix that contains in each element(to the corresponding axis) the coefficients of
    %each motion(5th order polinomial-linear-5th order polinomial-
    %5th order polinomial-linear-5th order polinomial). Each of
    %the element is a cell of k matrices(k indicates the number of motions e.g. k = 6)
    %with different size because of the different order of the
    %polynomials.
    
    %filling P(i) data
   
    P(i, 1:2*D/dt) = polyval(P_coeff{i}{1}, t(1:2*D/dt));    %for the time interval[0,2*D]
    P(i, 2*D/dt+1:(half_tf-2*D)/dt) = polyval(P_coeff{i}{2}, t(2*D/dt+1:(half_tf-2*D)/dt));  %for the time interval[2*D,half_tf-2*D]
    P(i, (half_tf-2*D)/dt+1:half_tf/dt) = polyval(P_coeff{i}{3}, t((half_tf-2*D)/dt+1:half_tf/dt));   %for the time interval[half_tf-2*D,half_tf]
    P(i, (half_tf)/dt+1:(half_tf+2*D)/dt) = polyval(P_coeff{i}{4}, t((half_tf)/dt+1:(half_tf+2*D)/dt));   %for the time interval[half_tf,half_tf+2*D]
    P(i, (half_tf+2*D)/dt+1:(tf-2*D)/dt) = polyval(P_coeff{i}{5}, t((half_tf+2*D)/dt+1:(tf-2*D)/dt));   %for the time interval[half_tf+2*D,tf-2*D]
    P(i, (tf-2*D)/dt+1:tf/dt+1) = polyval(P_coeff{i}{6}, t((tf-2*D)/dt+1:tf/dt+1));   %for the time interval[tf-2*D,tf]
    
    %differentiation of polynomial function P for i-axis
    V_coeff{i}{1} = polyder(P_coeff{i}{1});
    V_coeff{i}{2} = polyder(P_coeff{i}{2});
    V_coeff{i}{3} = polyder(P_coeff{i}{3});
    V_coeff{i}{4} = polyder(P_coeff{i}{4});
    V_coeff{i}{5} = polyder(P_coeff{i}{5});
    V_coeff{i}{6} = polyder(P_coeff{i}{6});
    %filling V(i) data
    V(i, 1:2*D/dt) = polyval(V_coeff{i}{1}, t(1:2*D/dt));    %for the time interval[0,2*D]
    V(i, 2*D/dt+1:(half_tf-2*D)/dt) = polyval(V_coeff{i}{2}, t(2*D/dt+1:(half_tf-2*D)/dt));   %for the time interval[2*D,half_tf-2*D]
    V(i, (half_tf-2*D)/dt+1:half_tf/dt) = polyval(V_coeff{i}{3}, t((half_tf-2*D)/dt+1:half_tf/dt));   %for the time interval[half_tf-2*D,half_tf]
    V(i, (half_tf)/dt+1:(half_tf+2*D)/dt) = polyval(V_coeff{i}{4}, t((half_tf)/dt+1:(half_tf+2*D)/dt));   %for the time interval[half_tf,half_tf+2*D]
    V(i, (half_tf+2*D)/dt+1:(tf-2*D)/dt) = polyval(V_coeff{i}{5}, t((half_tf+2*D)/dt+1:(tf-2*D)/dt));   %for the time interval[half_tf+2*D,tf-2*D]
    V(i, (tf-2*D)/dt+1:tf/dt+1) = polyval(V_coeff{i}{6}, t((tf-2*D)/dt+1:tf/dt+1));   %for the time interval[tf-2*D,tf]
    %differentiation of polynomial function V for i-axis
    A_coeff{i}{1} = polyder(V_coeff{i}{1});
    A_coeff{i}{2} = polyder(V_coeff{i}{2});
    A_coeff{i}{3} = polyder(V_coeff{i}{3});
    A_coeff{i}{4} = polyder(V_coeff{i}{4});
    A_coeff{i}{5} = polyder(V_coeff{i}{5});
    A_coeff{i}{6} = polyder(V_coeff{i}{6});
    %filling Vx data
    A(i, 1:2*D/dt) = polyval(A_coeff{i}{1}, t(1:2*D/dt));    %for the time interval[0,2*D]
    A(i, 2*D/dt+1:(half_tf-2*D)/dt) = polyval(A_coeff{i}{2}, t(2*D/dt+1:(half_tf-2*D)/dt));  %for the time interval[2*D,half_tf-2*D]
    A(i, (half_tf-2*D)/dt+1:half_tf/dt) = polyval(A_coeff{i}{3}, t((half_tf-2*D)/dt+1:half_tf/dt));   %for the time interval[half_tf-2*D,half_tf]
    A(i, (half_tf)/dt+1:(half_tf+2*D)/dt) = polyval(A_coeff{i}{4}, t((half_tf)/dt+1:(half_tf+2*D)/dt));  %for the time interval[half_tf,half_tf+2*D]
    A(i, (half_tf+2*D)/dt+1:(tf-2*D)/dt) = polyval(A_coeff{i}{5}, t((half_tf+2*D)/dt+1:(tf-2*D)/dt));  %for the time interval[half_tf+2*D,tf-2*D]
    A(i, (tf-2*D)/dt+1:tf/dt+1) = polyval(A_coeff{i}{6}, t((tf-2*D)/dt+1:tf/dt+1));   %for the time interval[tf-2*D,tf]
end
%As at the z-axis the motion is constant, we fill P,V,A(i) with constant
%values
P_coeff{z} = h;
P(z,:) = polyval(P_coeff{z}, t(:));
V_coeff{z} = polyder(P_coeff{z});
V(z,:) = polyval(V_coeff{z}, t(:));
A_coeff{z} = polyder(V_coeff{z});
A(z,:) = polyval(A_coeff{z}, t(:));

titles = {' X-axis' ; ' Y-axis' ; ' Z-axis'};   %cell array for easier indexing of the titles

for i = 1:3
    
    figure(fig_num);
    subplot(3,1,1); %we create subplot in order to display 3 different plots about each of the axes(x,y,z)
    plot(t(:), P(i,:));
    title(strcat('Trajectory of the Robot at the ' ,titles{i}));
    xlabel('Time(sec)');
    ylabel('Position(m)');
    
    subplot(3,1,2);
    plot(t(:), V(i,:));
    title(strcat('Velocity of the Robot at the ', titles{i}));
    xlabel('Time(sec)');
    ylabel('Velocity(m/sec)');
    
    subplot(3,1,3);
    plot(t(:), A(i,:));
    title(strcat('Acceleration of the Robot at the ', titles{i}));
    xlabel('Time(sec)');
    ylabel('Acceleration(m/sec^2)');
    %saveas(figure(fig_num), strcat(num2str(fig_num),'.jpg'));
    fig_num = fig_num + 1;
end
    
%Task B.7b
%Such we have mentioned to the report, we define the distance Pxz
%and A23 = +- Pxz
Pxz = sqrt(P(x,:).^2 + P(z,:).^2 - l2^2);
A23 = [Pxz ; -Pxz];
%We are going to calculate the values of cos and sin for each angle
%and then the values of each angle, corresponding to the 4 different
%solutions, the robotic configuration can approach
%We have to remind that:
%q3 -> 2 solutions(2 double)
%q2 -> 4 solutions 
%q1 -> 2 solutions 
%Totally we are going to elaborate the 4 different solutions for the 
%given movement
%First Solution: A23(1,:) = Pxz(:) -> 
%   { (c1(1,:),s1(1,:)) , (c2(1,:),s2(1,:)) ,(c3(1,:),s3(1,:)) } ->
%   { q1(1,:) , q2(1,:) , q3(1,:) }    
%Second Solution: A23(1,:) = Pxz(:) -> 
%   { (c1(1,:),s1(1,:)) , (c2(2,:),s2(2,:)) ,(c3(2,:),s3(2,:)) } ->
%   { q1(1,:) , q2(2,:) , q3(2,:) }    
%Third Solution: A23(2,:) = -Pxz(:) -> 
%   { (c1(2,:),s1(2,:)) , (c2(3,:),s2(3,:)) ,(c3(1,:),s3(1,:)) } ->
%   { q1(2,:) , q2(3,:) , q3(1,:) }    
%Fourth Solution: A23(2,:) = -Pxz(:) -> 
%   { (c1(2,:),s1(2,:)) , (c2(4,:),s2(4,:)) ,(c3(2,:),s3(2,:)) } ->
%   { q1(2,:) , q2(4,:) , q3(2,:) }    
c3 = zeros(length(t)); s3 = zeros(2,length(t)); 
c2 = zeros(4,length(t)); s2 = zeros(4,length(t)); 
c1 = zeros(2,length(t)); s1 = zeros(2,length(t)); 
q3 = zeros(2, length(t));
q2 = zeros(4, length(t));
q1 = zeros(2, length(t));
%%%%% We start with A23(1) = Pxz %%%%%
for i = 1:length(t)
    %analysis for q3
    c3(i) = (P(y,i)^2 + A23(1,i)^2 - l5^2 - l4^2)/(2*l4*l5);
    s3(1,i) = sqrt(1 - c3(i)^2);
    s3(2,i) = -sqrt(1 - c3(i)^2);
    q3(1,i) = atan2(s3(1,i),c3(i));    %first solution of q3
    q3(2,i) = atan2(s3(2,i),c3(i));    %second solution of q3
end
k = 0;
%we use variablr k in order to index the corresponding 
%rows of the arrays
for j = 1:2
    for i = 1:length(t)
        %analysis for q1
        s1(j,i) = (P(z,i)*l2 +P(x,i)*A23(j,i))/(A23(j,i)^2 + l2^2);    
        c1(j,i) = (P(x,i) - s1(j,i)*A23(j,i))/l2;
        q1(j,i) = atan2(s1(j,i),c1(j,i));    
        %analysis for q2
        c2(j+k,i) = (P(y,i)*s3(1,i)*l5+A23(j,i)*(c3(i)*l5+l4))/((c3(i)*l5+l4)^2 + (s3(1,i)*l5)^2);
        c2(j+k+1,i) = (P(y,i)*s3(2,i)*l5+A23(j,i)*(c3(i)*l5+l4))/((c3(i)*l5+l4)^2 + (s3(2,i)*l5)^2);
        s2(j+k,i) = (c2(j+k,i)*(c3(i)*l5+l4)-A23(j,i))/(s3(1,i)*l5);
        s2(j+k+1,i) = (c2(j+k+1,i)*(c3(i)*l5+l4)-A23(j,i))/(s3(2,i)*l5);
        q2(j+k,i) = atan2(s2(j+k,i), c2(j+k,i));       
        q2(j+k+1,i) = atan2(s2(j+k+1,i), c2(j+k+1,i));       
    end
    k = 1;
    
end

%Titles and messages to make code and plots more interactive
titles = {' for A23 = Pxz',' for A23 = -Pxz' };
solutions = {' (1st Solution)', ' (2nd Solution)', ' (3rd Solution)', ' (4th Solution)' };
disp('Now we are gonna plot the four different solutions for the motion of angles');

for i = 1:4
    
    if i == 1 || i == 2
        j = 1;
    else
        j = 2;
    end
    
    figure(fig_num);
    subplot(3,1,1);
    plot(t, q1(j,:));
    title(strcat('Motion of Angle q1 ', titles{j}, solutions{i}));
    ylabel('Amplitude(rad)');
    xlabel('Time(sec)');
    
    subplot(3,1,2);
    plot(t, q2(i,:));
    title(strcat('Motion of Angle q2 ', titles{j}, solutions{i}));
    ylabel('Amplitude(rad)');
    xlabel('Time(sec)');
    
    subplot(3,1,3);
    plot(t, q3(mod(i-1,2)+1,:));
    title(strcat('Motion of Angle q3 ', titles{j}, solutions{i}));
    ylabel('Amplitude(rad)');
    xlabel('Time(sec)');
    %saveas(figure(fig_num), strcat(num2str(fig_num),'.jpg'));
    fig_num = fig_num + 1;
end

%%%% Inv-Diff Kinematic for Angles %%%%
%For each solution, we have the coresponding velocity
q1dot = zeros(2,length(t));
q2dot = zeros(4,length(t));
q3dot = zeros(4,length(t));

%%%% Computation of Angles Velocity %%%%
%First we calculate the  elements of the J^(-1) from the inverse-diff kinematics 
%As we have 4 solutions, there are different values of each element
%according to the variables , that subject to.
Jinv_11 = zeros(2,length(t)); Jinv_21 = zeros(4,length(t)); Jinv_31 = zeros(4,length(t)); 
Jinv_12 = zeros(2,length(t)); Jinv_22 = zeros(4,length(t)); Jinv_32 = zeros(4,length(t)); 
Jinv_13 = zeros(2,length(t)); Jinv_23 = zeros(4,length(t)); Jinv_33 = zeros(4,length(t));
k = 0;
for j = 1:2
    for i = 1:length(t)

        Jinv_11(j,i) = c1(j,i)/A23(j,i);
        Jinv_12(j,i) = 0;
        Jinv_13(j,i) = s1(j,i)/A23(j,i);
        Jinv_21(j+k,i) = (c2(j+k,i)*c3(i)-s2(j+k,i)*s3(1,i))*(s1(j,i)*A23(j,i)+c1(j,i)*l2)/(A23(j,i)*s3(1,i)*l4);
        Jinv_21(j+k+1,i) = (c2(j+k+1,i)*c3(i)-s2(j+k+1,i)*s3(2,i))*(s1(j,i)*A23(j,i)+c1(j,i)*l2)/(A23(j,i)*s3(2,i)*l4);
        Jinv_22(j+k,i) = (s2(j+k,i)*c3(i)+c2(j+k,i)*s3(1,i))/(s3(1,i)*l4);
        Jinv_22(j+k+1,i) = (s2(j+k+1,i)*c3(i)+c2(j+k+1,i)*s3(2,i))/(s3(2,i)*l4);
        Jinv_23(j+k,i) = (c2(j+k,i)*c3(i)-s2(j+k,i)*s3(1,i))*(c1(j,i)*A23(j,i)-s1(j,i)*l2)/(A23(j,i)*s3(1,i)*l4);
        Jinv_23(j+k+1,i) = (c2(j+k+1,i)*c3(i)-s2(j+k+1,i)*s3(2,i))*(c1(j,i)*A23(j,i)-s1(j,i)*l2)/(A23(j,i)*s3(2,i)*l4);
        Jinv_31(j+k,i) = -(s1(j,i)*A23(j,i)+c1(j,i)*l2)/(s3(1,i)*l4*l5);
        Jinv_31(j+k+1,i) = -(s1(j,i)*A23(j,i)+c1(j,i)*l2)/(s3(2,i)*l4*l5);
        Jinv_32(j+k,i) = -((s2(j+k,i)*c3(i)+c2(j+k,i)*s3(1,i))*l5+s2(j+k,i)*l4)/(s3(1,i)*l5*l4);
        Jinv_32(j+k+1,i) = -((s2(j+k+1,i)*c3(i)+c2(j+k+1,i)*s3(2,i))*l5+s2(j+k+1,i)*l4)/(s3(2,i)*l5*l4);
        Jinv_33(j+k,i) = (c1(j,i)*A23(j,i)-s1(j,i)*l2)/(s3(1,i)*l4*l5);
        Jinv_33(j+k+1,i) = (c1(j,i)*A23(j,i)-s1(j,i)*l2)/(s3(2,i)*l4*l5);

    end
    k=1;
end
%Then we calculate the coresponding values of the angles
k = 0;
for i = 1:2  
    for j = 1:length(t)
        
        q1dot(i,j) =  Jinv_11(i,j)*V(x,j) + Jinv_12(i,j)*V(x,j) + Jinv_13(i,j)*V(z,j);
        q2dot(i+k,j) =  Jinv_21(i+k,j)*V(x,j) + Jinv_22(i+k,j)*V(y,j) + Jinv_23(i+k,j)*V(z,j);
        q2dot(i+k+1,j) =  Jinv_21(i+k+1,j)*V(x,j) + Jinv_22(i+k+1,j)*V(y,j) + Jinv_23(i+k+1,j)*V(z,j);
        q3dot(i+k,j) = Jinv_31(i+k,j)*V(x,j) + Jinv_32(i+k,j)*V(y,j) + Jinv_33(i+k,j)*V(z,j);
        q3dot(i+k+1,j) = Jinv_31(i+k+1,j)*V(x,j) + Jinv_32(i+k+1,j)*V(y,j) + Jinv_33(i+k+1,j)*V(z,j);

    end
    k = 1;
end

disp('Now we are gonna plot the four corresponding velocities of angles');
for i = 1:4
    
    if i == 1 || i == 2
        j = 1;
    else
        j = 2;
    end
    
    figure(fig_num);
    subplot(3,1,1);
    plot(t, q1dot(j,:));
    title(strcat('Velocity of Angle q1 ', titles{j}, solutions{i}));
    ylabel('Amplitude(rad/sec)');
    xlabel('Time(sec)');
    
    subplot(3,1,2);
    plot(t, q2dot(i,:));
    title(strcat('Velocity of Angle q2 ', titles{j}, solutions{i}));
    ylabel('Amplitude(rad/sec)');
    xlabel('Time(sec)');
    
    subplot(3,1,3);
    plot(t, q3dot(i,:));
    title(strcat('Velocity of Angle q3 ', titles{j}, solutions{i}));
    ylabel('Amplitude(rad/sec)');
    xlabel('Time(sec)');
    %saveas(figure(fig_num), strcat(num2str(fig_num),'.jpg'));
    fig_num = fig_num + 1;
end

%%%% Task B.7c %%%%
%%%% We are going to define the coordinates of the each link %%%%
%%%   Because of the multiple solutions we are computing different %%%
%%%   values of each link   %%%
xd0 = zeros(2,length(t));
yd0 = zeros(2,length(t));
zd0 = zeros(2,length(t));
xd1 = zeros(2,length(t)); 
yd1 = zeros(2,length(t)); 
zd1 = zeros(2,length(t)); 
xd2 = zeros(4,length(t)); 
yd2 = zeros(4,length(t)); 
zd2 = zeros(4,length(t));
xd3 = zeros(4,length(t));
yd3 = zeros(4,length(t));
zd3 = zeros(4,length(t));

%%%% Coordinates of the First Link l2 %%%%
for i = 1:2
    for j = 1:length(t)
    xd1(i,j) = cos(q1(i,j))*l2;
    zd1(i,j) = sin(q1(i,j))*l2;
    end
end
k = 0;
for j = 1:2

    for i = 1:length(t)
    %%%% Coordinates of the Second Link l4 %%%%
        xd2(j+k,i) = sin(q1(j,i))*cos(q2(j+k,i))*l4 + cos(q1(j,i))*l2;
        xd2(j+k+1,i) = sin(q1(j,i))*cos(q2(j+k+1,i))*l4 + cos(q1(j,i))*l2;
        yd2(j+k,i) = sin(q2(j+k,i))*l4;
        yd2(j+k+1,i) = sin(q2(j+k+1,i))*l4;
        zd2(j+k,i) = -cos(q1(j,i))*cos(q2(j+k,i))*l4 + sin(q1(j,i))*l2;
        zd2(j+k+1,i) = -cos(q1(j,i))*cos(q2(j+k+1,i))*l4 + sin(q1(j,i))*l2;
        
    %%%% Coordinates of the Third Link l5 are the same as the end-effector's%%%%
        xd3(j+k,i) = sin(q1(j,i))*((cos(q2(j+k,i))*cos(q3(1,i))-sin(q2(j+k,i))*sin(q3(1,i)))*l5+cos(q2(j+k,i))*l4)+cos(q1(j,i))*l2;
        xd3(j+k+1,i) = sin(q1(j,i))*((cos(q2(j+k+1,i))*cos(q3(2,i))-sin(q2(j+k+1,i))*sin(q3(2,i)))*l5+cos(q2(j+k+1,i))*l4)+cos(q1(j,i))*l2;
        yd3(j+k,i) =  (sin(q2(j+k,i))*cos(q3(1,i))+cos(q2(j+k,i))*sin(q3(1,i)))*l5+sin(q2(j+k,i))*l4;
        yd3(j+k+1,i) = (sin(q2(j+k+1,i))*cos(q3(2,i))+cos(q2(j+k+1,i))*sin(q3(2,i)))*l5+sin(q2(j+k+1,i))*l4;
        zd3(j+k,i) = -cos(q1(j,i))*((cos(q2(j+k,i))*cos(q3(1,i))-sin(q2(j+k,i))*sin(q3(1,i)))*l5+cos(q2(j+k,i))*l4)+sin(q1(j,i))*l2;
        zd3(j+k+1,i) = -cos(q1(j,i))*((cos(q2(j+k+1,i))*cos(q3(2,i))-sin(q2(j+k+1,i))*sin(q3(2,i)))*l5+cos(q2(j+k+1,i))*l4)+sin(q1(j,i))*l2;
        
    end
    k = 1;
end

