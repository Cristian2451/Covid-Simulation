clear
clc

N = input("number of individuals:");
%Initially there will be 1 infected and the rest are healthy
%An array "Indiv" of length N will be created which will indicate the health status of individuals
%Where Healthy=0 Sick=1 Infected=0.5 and Recovered=-1 !!!!!
Indiv = cat(2,ones(1),zeros(1,N-1)); 

%Define Parameters
W = 1000;      %Width of space (m)
H = 1000;      %Height of space (m)
v_min = 0.1;   %Minimum velocity of individuals (m/s)
v_max = 0.2;   %Maximum velocity of individuals (m/s)
T=10*24*3600;  %Total time for simulation (s)
delta_t = 50;  %Time step for simulation (s)

%Define Initial Variables
x = W*rand(1,N);                   %Initial x-position of individuals (m)
y = H*rand(1,N);                   %Initial y-position of individuals (m)
V = v_min+(v_max-v_min)*rand(1,N); %Initial speed of individuals (m/s)
theta = 2*pi*rand(1,N);            %Initial direction of individuals (m/s)
u = V.*cos(theta);                 %Horizontal velocity component (m/s)
v = V.*sin(theta);                 %Vertical velocity component (m/s)

%Pre-allocation of variables that change size on every loop iteration to optimise code
Healthy_info = zeros(T/delta_t,1);
Infected_info = zeros(T/delta_t,1);
Sick_info = zeros(T/delta_t,1);
Recovered_info = zeros(T/delta_t,1);
x(2:T/delta_t+1,N) = zeros(T/delta_t,1);
y(2:T/delta_t+1,N) = zeros(T/delta_t,1);
t = zeros(T/delta_t,1);
Inf_t = zeros(1,N); %Time of infection of individuals (s)
Inf_t(1) = -2*24*3600; %Time of infection of initial sick individual
History = struct('Healthy',cell(1,T/delta_t),'Infected',cell(1,T/delta_t), ...
    'Sick',cell(1,T/delta_t),'Recovered',cell(1,T/delta_t),'Time',cell(1,T/delta_t));

%Simulation has to iterate through time in steps of delta_t until total time T is reached
%Therefore for-loop has to iterate from 1 to T/delta_t to achieve all the required iterations
for j = 1:T/delta_t
    t(j,1)=j*delta_t; %Time elapsed (s)
    [~,~,d,h,m,s] = datevec(seconds(t(j,1))); %Time elapsed in days, hours, minutes and seconds
    
    Inf_sick = find(Indiv == 1);     %Find sick individuals
    Inf_asym = find(Indiv == 0.5);   %Find Infected individuals
    Inf = cat(2,Inf_sick,Inf_asym);  %Store positions from Indiv array of sick and infected individuals
    
    %Iterate through each individual to update position, direction, health status and plot individuals
    for i = 1:N
        
        %Iterate through each sick/infected individual
        for l = 1:length(Inf)
            %Find gap between sick/infected individual and individual
            gap = sqrt((x(j,Inf(l))-x(j,i))^2+(y(j,Inf(l))-y(j,i))^2); % d=sqrt((x1-x2)^2+(y1-y2)^2)
            
            if Indiv(i)==0 && gap<=2 %If individual healthy and gap<=2m
                w = rand;
                if Indiv(Inf(l)) == 1 && w > 0.5        %If infected individual is sick and probability of infection (50%) satisfied
                    Indiv(i) = 0.5;                     %Update health status of healthy individual to infected 
                    Inf_t(i) = j*delta_t;               %Record time of infection
                elseif Indiv(Inf(l)) == 0.5 && w > 0.7  %If infected individual is only infected and probability of infection (30%) satisfied
                    Indiv(i) = 0.5;                     %Update health status of healthy individual to infected 
                    Inf_t(i) = j*delta_t;               %Record time of infection
                end  
            end
        end
        
        %Update health status of infected individuals depending on time elapsed
        if Indiv(i) == 0.5 %If individual infected
            if t(j,1) == (Inf_t(i) + 2*24*3600) %If 2 days elapsed since infection, individual becomes sick
                Indiv(i) = 1;
            end   
        elseif Indiv(i) == 1 %If individual sick    
            if t(j,1) == (Inf_t(i) + 5*24*3600) %If 5 days elapsed since infection, individual becomes recovered
                Indiv(i) = -1;
            end
        end
        %Update x and y position of individual
        x(j+1,i) = x(j,i) + u(i)*delta_t; %Xnew=Xold+u*delta_t
        y(j+1,i) = y(j,i) + v(i)*delta_t; %Ynew=Yold+v*delta_t
        %Boundary function updates direction of individual if comes within 2m of wall 
        [u(i),v(i)] = boundary(x(j+1,i),y(j+1,i),W,H,u(i),v(i),V(i));
        
        %Plot positions of all individuals every 6 minutes (every 36 j's) to optimize code
        if mod(j,36) == 0 || j == 1 
%             if i == 1
%                 subplot(2,2,1); %Only call subplot when plotting first individual for code efficiency
%             end
            if Indiv(i) == 1 %if sick
                plot(x(j,i),y(j,i),'o','MarkerEdgeColor','r','MarkerFaceColor','r');
            elseif Indiv(i) == 0 %if healthy
                plot(x(j,i),y(j,i),'o','MarkerEdgeColor','g','MarkerFaceColor','g');
            elseif Indiv(i) == 0.5 %if infected
                plot(x(j,i),y(j,i),'o','MarkerEdgeColor',[1 0.5 0],'MarkerFaceColor',[1 0.5 0]); %RGB colour is orange
            elseif Indiv(i) == -1 %if recovered
                plot(x(j,i),y(j,i),'o','MarkerEdgeColor','b','MarkerFaceColor','b');
            end
            %hold on, to plot all individuals in same figure
            %Call 'hold on' only when first individual plotted to optimize code
            if i == 1
                hold on
            end
        end
    end
    %Once all individuals plotted 
    if mod(j,36) == 0 || j == 1
        set(gca,'Color','k')
        set(gca,'xlim',[0 1000],'ylim',[0 1000]); %Set limits of axis
        ax = gca; % Get handle to current axes.
        ax.XColor = 'w'; % Red
        ax.YColor = 'w'; % Blue
        title(['time elapsed = ', num2str(d),' days and ', num2str(h), ' hours'],'FontSize',14,'Color','w'); %Add title indicating days and hours elapsed
        hold off
        pause(0.01) %Pauses the plot by 0.01s
    end
    
    %Create array for each health status, recording total number of individuals with that health status at each time step
    Healthy_info(j,1) = length(find(Indiv == 0));
    Infected_info(j,1) = length(find(Indiv == 0.5));
    Sick_info(j,1) = length(find(Indiv == 1));
    Recovered_info(j,1) = length(find(Indiv == -1));
    
    %Stucture array 'History' with fields for time elapsed and every health status at every time step
    History(j).Healthy = Healthy_info(j,1);      % n° of healthy individuals
    History(j).Infected = Infected_info(j,1);    % n° of infected individuals
    History(j).Sick = Sick_info(j,1);            % n° of sick individuals
    History(j).Recovered = Recovered_info(j,1);  % n° of recovered individuals
    History(j).Time = [d,h,m,s];                 %Total number of days, hours, minutes and seconds elapsed
    
%     if mod(j,36) == 0 || j == 1
%         %Bar plot of health status
%         if j == 1
%             sub2 = subplot(2,2,2); %Calling subplot only at first time step
%         end
%         barchart = bar(sub2,[History(j).Healthy; History(j).Infected; History(j).Sick; History(j).Recovered],'FaceColor','flat');
%         barchart.CData = [0 1 0 ; 1.0 0.5 0.0 ; 1 0 0 ; 0 0 1]; %Add RGB colours to bars 
%         set(sub2,'xticklabel',{'Healthy';'Infected';'Sick';'Recovered'},'ylim',[0 N],'FontSize',12); %Label bars and set limit to axis        
%         ylabel(sub2,'Number of Individuals','FontSize',12); %Label y axis of barchart
%         title(sub2,['time elapsed = ', num2str(d),' days and ', num2str(h), ' hours'],'FontSize',14); %Add title to bar chart
%         
%         if j == 1
%             %Set up animated lines for line plot of health status
%             %Animatedline used instead of plot as it is more code efficient 
%             sub3 = subplot(2,2,[3 4]);
%             p1 = animatedline('Color','g','LineWidth',2);        %Healthy individuals
%             p2 = animatedline('Color',[1 0.5 0],'LineWidth',2);  %Infected individuals
%             p3 = animatedline('Color','r','LineWidth',2);        %Sick individuals
%             p4 = animatedline('Color','b','LineWidth',2);        %Recovered individuals
%             %Label axes
%             xlabel('time (hours)','FontSize',12);
%             ylabel('Number of Individuals','FontSize',12);
%             
%             set(gca,'xlim',[0 T/3600],'ylim',[0 N]); %Set limits of axes
%             set(gcf,'position',[100 500 1200 800]); %Set dimensions of subplot
%             legend('Healthy','Infected','Sick','Recovered','Location','EastOutside','FontSize',14); %Adds legend to line plot
%         end
%         
%         title(sub3,['time elapsed = ', num2str(d),' days and ', num2str(h), ' hours'],'FontSize',14); %Adds title to line plot
%         %Add points to each animatedline
%         addpoints(p1,t(j,1)/3600,History(j).Healthy);
%         addpoints(p2,t(j,1)/3600,History(j).Infected);
%         addpoints(p3,t(j,1)/3600,History(j).Sick);
%         addpoints(p4,t(j,1)/3600,History(j).Recovered);
%         drawnow limitrate 
%     end
    
    %Save subplot after 0, 2, 4 and 6 days
    if j == 1 || j == 17280 || j == 34560 || j == 51840
        filename = ['Subplot_' num2str(d) '_days.elapsed' '.jpeg']; %Set file name and set file as jpeg
        saveas(figure(1),filename)
    end
end

Day = [1;2;3;4;5;6;7;8;9;10];
%Find n° of individuals with each health status at the end of each day
k = 8640; %One day in j's
Healthy = Healthy_info(k:k:k*10);
Infect = Infected_info(k:k:k*10);
Sick = Sick_info(k:k:k*10);
Recovered = Recovered_info(k:k:k*10);

%Create table that displays the total n° of individuals with each health status at the end of each day
DailySummaries = table(Day,Healthy,Infect,Sick,Recovered);
disp(DailySummaries)
%Save table as a txt
writetable(DailySummaries,'DailySummaries.txt','Delimiter','\t');
