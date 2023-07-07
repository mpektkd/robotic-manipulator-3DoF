x_lower = zeros(4);
x_upper = zeros(4);
y_lower = zeros(4);
y_upper = zeros(4);
z_lower = zeros(4);
z_upper = zeros(4);
disp('Now we are gonna plotting the Animation of the Robot');

dtk=tf/dt/150; %% plot robot position every dtk samples, to animate its motion 
titles = {' 1st Solution',' 2nd Solution',' 3rd Solution',' 4th Solution' };
j = 1;
%%%% Now we are going to plot the 4 differents solutions in: %%%%
    %%%%    1)xyz space   %%%%
    %%%%    2)xz plane    %%%%
    %%%%    3)xy plane    %%%%
    %%%%    4)yz plane    %%%%
    
for i = 1:4
    if i == 3 || i == 4
        j = 2;
    end
%%%% Find the boundaries of each plot %%%%
    x_lower(i) = (min([xd0(j,:)';xd1(j,:)';xd2(i,:)';xd3(i,:)']));
    x_upper(i) = (max([xd0(j,:)';xd1(j,:)';xd2(i,:)';xd3(i,:)']));
    y_lower(i) = (min([yd0(j,:)';yd1(j,:)';yd2(i,:)';yd3(i,:)']));
    y_upper(i) = (max([yd0(j,:)';yd1(j,:)';yd2(i,:)';yd3(i,:)']));
    z_lower(i) = (min([zd0(j,:)';zd1(j,:)';zd2(i,:)';zd3(i,:)']));
    z_upper(i) = (max([zd0(j,:)';zd1(j,:)';zd2(i,:)';zd3(i,:)']));

    
    figure(fig_num); 
    fig_num = fig_num + 1;
    grid();
    axis([x_lower(i) x_upper(i) y_lower(i) y_upper(i) z_lower(i) z_upper(i)])
    axis on 
    hold on %hold the image
    xlabel('x (cm)'); 
    ylabel('y (cm)'); 
    zlabel('z (cm)'); 
    plot3(P(x,:),P(y,:),P(z,:),'rs'); 
    title(strcat('Animation of the Robot Motion in x-y-z space for the', titles{i}));
    for tk=1:dtk:tf/dt   
       pause(0.01);	%% pause motion to view successive robot configurations    
       plot3([0],[0],[0],'ro'); 
       
       plot3([0,xd1(j,tk)], [0,yd1(j,tk)], [0,zd1(j,tk)]);					
       plot3([xd1(j,tk)],[yd1(j,tk)], [zd1(j,tk)], 'ko');    
       plot3([xd1(j,tk),xd2(i,tk)],[yd1(j,tk),yd2(i,tk)],[zd1(j,tk),zd2(i,tk)]);	
       plot3([xd2(i,tk)],[yd2(i,tk)], [zd2(i,tk)], 'go'); 
       plot3([xd2(i,tk),xd3(i,tk)],[yd2(i,tk),yd3(i,tk)],[zd2(i,tk),zd3(i,tk)]);	
       plot3([xd3(i,tk)],[yd3(i,tk)], [zd3(i,tk)],'y*');    
    end       
    %saveas(figure(fig_num), strcat(num2str(fig_num),'.jpg'));
    hold off 
 
    figure(fig_num); 
    fig_num = fig_num + 1;
    grid();
    axis([x_lower(i) x_upper(i) z_lower(i) z_upper(i)]) 
    axis on 
    hold on 
    xlabel('x (cm)'); 
    ylabel('z (cm)'); 
    plot(P(x,:),P(z,:),'rs'); 

    title(strcat('Animation of the Robot Motion in x-z plane for the', titles{i}));
    for tk=1:dtk:tf/dt  	
       pause(0.01);	%% pause motion to view successive robot configurations    
       plot([0],[0],'ro'); 
       plot([0,xd1(j,tk)],[0,zd1(j,tk)]);					
       plot([xd1(j,tk)],[zd1(j,tk)], 'ko');    
       plot([xd1(j,tk),xd2(i,tk)],[zd1(j,tk),zd2(i,tk)]);	
       plot([xd2(i,tk)],[zd2(i,tk)], 'go'); 
       plot([xd2(i,tk),xd3(i,tk)],[zd2(i,tk),zd3(i,tk)]);	
       plot([xd3(i,tk)],[zd3(i,tk)],'y*');    
    end       
    saveas(figure(fig_num), strcat(num2str(fig_num),'.jpg'));
    hold off
    
    figure(fig_num); 
    fig_num = fig_num + 1;
    grid();
    axis([x_lower(i) x_upper(i) y_lower(i) y_upper(i)]) 
    axis on 
    hold on 
    xlabel('x (cm)'); 
    ylabel('y (cm)'); 
    plot(P(x,:),P(y,:),'rs'); 


    title(strcat('Animation of the Robot Motion in x-y plane for the', titles{i}));
    for tk=1:dtk:tf/dt   	
       pause(0.01);	%% pause motion to view successive robot configurations    
       plot([0],[0],'ro'); 
       plot([0,xd1(j,tk)], [0,yd1(j,tk)]);					
       plot([xd1(j,tk)],[yd1(j,tk)], 'ko');    
       plot([xd1(j,tk),xd2(i,tk)],[yd1(j,tk),yd2(i,tk)]);	
       plot([xd2(i,tk)],[yd2(i,tk)], 'go'); 
       plot([xd2(i,tk),xd3(i,tk)],[yd2(i,tk),yd3(i,tk)]);	
       plot([xd3(i,tk)],[yd3(i,tk)],'y*');    
    end       
    %saveas(figure(fig_num), strcat(num2str(fig_num),'.jpg'));
    hold off
    
    figure(fig_num); 
    fig_num = fig_num + 1;
    grid();
    axis([y_lower(i) y_upper(i) z_lower(i) z_upper(i)]) 
    axis on 
    hold on 
    xlabel('y (cm)'); 
    ylabel('z (cm)'); 
    plot(P(y,:),P(z,:),'rs'); 
    
    title(strcat('Animation of the Robot Motion in y-z plane for the', titles{i}));
    for tk=1:dtk:tf/dt   
       pause(0.01);	%% pause motion to view successive robot configurations    
       plot([0],[0],'ro'); 
       plot([0,yd1(j,tk)], [0,zd1(j,tk)]);					
       plot([yd1(j,tk)], [zd1(j,tk)], 'ko');    
       plot([yd1(j,tk),yd2(i,tk)],[zd1(j,tk),zd2(i,tk)]);	
       plot([yd2(i,tk)], [zd2(i,tk)], 'go'); 
       plot([yd2(i,tk),yd3(i,tk)],[zd2(i,tk),zd3(i,tk)]);	
       plot([yd3(i,tk)], [zd3(i,tk)],'y*');    
    end
    %saveas(figure(fig_num), strcat(num2str(fig_num),'.jpg'));
    hold off
    
end


  