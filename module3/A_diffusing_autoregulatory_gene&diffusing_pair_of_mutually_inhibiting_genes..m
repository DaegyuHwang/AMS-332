%%%% Part A

%%%%%%%%%%%%%%%
%% 1
%%%%%%%%%%%%%%%

u=1; w=1; CXprot=1; CXrna=1; K=0.33; Dprot=1*10^(-4); Drna=1*10^(-4);
%over space and time.
totaltime = 30;
totaldistance = 3;
deltat = 0.01;
deltad = 0.02;

numtime = totaltime/deltat;
numposition = totaldistance/deltad;
timearray = zeros(numtime,1);
distancearray = zeros(numposition,1);
Xrna = zeros(numtime, numposition);
Xprot = zeros(numtime, numposition);
Xrna(1,:) = 0.5;
Xprot(1,:) = 0.5;

for t = 1:numtime-1
    for x = 1:numposition
        if x == 1
            laplacian_xp = (Xprot(t,x+1)-2*Xprot(t,x))/(deltad^2);
            laplacian_xr = (Xrna(t,x+1)-2*Xrna(t,x))/(deltad^2);
        elseif x == numposition
            laplacian_xp = (Xprot(t,x-1)-2*Xprot(t,x))/(deltad^2);
            laplacian_xr = (Xrna(t,x-1)-2*Xrna(t,x))/(deltad^2);
        else
            laplacian_xp = (Xprot(t,x+1)+Xprot(t,x-1)-2*Xprot(t,x))/(deltad^2);
            laplacian_xr = (Xrna(t,x+1)+Xrna(t,x-1)-2*Xrna(t,x))/(deltad^2);
        end
        
        dXrnadt = u*(Xprot(t,x)^2)/(K^2+Xprot(t,x)^2)-CXrna*Xrna(t,x)+Drna*laplacian_xr;
        dXprotdt = w*Xrna(t,x)-CXprot*Xprot(t,x)+Dprot*laplacian_xp;
        Xrna(t+1,x) = Xrna(t,x) + dXrnadt*deltat;
        Xprot(t+1,x) = Xprot(t,x) + dXprotdt*deltat;
        
        if x ~= numposition
            distancearray(x+1) = distancearray(x) + deltad;
        end
    end
    timearray(t+1) = timearray(t) + deltat;    
end

figure(1)
subplot(2,2,1)
hold on
plot(timearray, Xrna(:,1))
plot(timearray, Xrna(:,2))
plot(timearray, Xrna(:,3))
hold off
xlabel('Time (s)'); ylabel('Concentration (mM)')
title('A diffusing autoregulatory gene (con vs time)')
legend('RNA first','RNA second','RNA third')

subplot(2,2,2)
hold on 
plot(timearray, Xprot(:,1))
plot(timearray, Xprot(:,2))
plot(timearray, Xprot(:,3))
hold off
xlabel('Time (s)'); ylabel('Concentration (mM)')
title('A diffusing autoregulatory gene (con vs time)')
legend('Protein first','Protein second','Protein third')

subplot(2,2,3)
hold on
plot(distancearray, Xrna(1,:))
plot(distancearray, Xrna(numtime,:))
hold off
xlabel('Distance (um)'); ylabel('Concentration (mM)')
legend('RNA first','RNA last')
title('A diffusing autoregulatory gene (con vs dis)')

subplot(2,2,4)
hold on
plot(distancearray, Xprot(1,:))
plot(distancearray, Xprot(numtime,:))
hold off
xlabel('Distance (um)'); ylabel('Concentration (mM)')
title('A diffusing autoregulatory gene (con vs dis)')
legend('Protein first','Protein last')
%%%%%%%%%%%%%%%
%% 2
%%%%%%%%%%%%%%%

u=1; w=1; CXprot=1; CXrna=1; K=0.33; Dprot=1*10^(-4); Drna=1*10^(-4);
totaltime = 30;
totaldistance = 3;
deltat = 0.01;
deltad = 0.02;

numtime = totaltime/deltat;
numposition = totaldistance/deltad;
timearray = zeros(numtime,1);
distancearray = zeros(numposition,1);
Xrna = zeros(numtime, numposition);
Xprot = zeros(numtime, numposition);
Xrna(1,1) = 0.5;
Xprot(1,1) = 0.5;

for t = 1:numtime-1
    for x = 1:numposition
        if x == 1
            laplacian_xp = (Xprot(t,x+1)-2*Xprot(t,x))/(deltad^2);
            laplacian_xr = (Xrna(t,x+1)-2*Xrna(t,x))/(deltad^2);
        elseif x == numposition
            laplacian_xp = (Xprot(t,x-1)-2*Xprot(t,x))/(deltad^2);
            laplacian_xr = (Xrna(t,x-1)-2*Xrna(t,x))/(deltad^2);
        else
            laplacian_xp = (Xprot(t,x+1)+Xprot(t,x-1)-2*Xprot(t,x))/(deltad^2);
            laplacian_xr = (Xrna(t,x+1)+Xrna(t,x-1)-2*Xrna(t,x))/(deltad^2);
        end
        
        dXrnadt = u*(Xprot(t,x)^2)/(K^2+Xprot(t,x)^2)-CXrna*Xrna(t,x)+Drna*laplacian_xr;
        dXprotdt = w*Xrna(t,x)-CXprot*Xprot(t,x)+Dprot*laplacian_xp;
        Xrna(t+1,x) = Xrna(t,x) + dXrnadt*deltat;
        Xprot(t+1,x) = Xprot(t,x) + dXprotdt*deltat;
        
        if x ~= numposition
            distancearray(x+1) = distancearray(x) + deltad;
        end
    end
    timearray(t+1) = timearray(t) + deltat;    
end

figure(1)
subplot(2,2,1)
hold on
plot(timearray, Xrna(:,1))
plot(timearray, Xrna(:,2))
plot(timearray, Xrna(:,3))
hold off
xlabel('Time (s)'); ylabel('Concentration (mM)')
title('A diffusing autoregulatory gene (con vs time)')
legend('RNA first','RNA second','RNA third')

subplot(2,2,2)
hold on 
plot(timearray, Xprot(:,1))
plot(timearray, Xprot(:,2))
plot(timearray, Xprot(:,3))
hold off
xlabel('Time (s)'); ylabel('Concentration (mM)')
title('A diffusing autoregulatory gene (con vs time)')
legend('Protein first','Protein second','Protein third')

subplot(2,2,3)
hold on
plot(distancearray, Xrna(1,:))
plot(distancearray, Xrna(numtime,:))
hold off
xlabel('Distance (um)'); ylabel('Concentration (mM)')
legend('RNA first','RNA last')
title('A diffusing autoregulatory gene (con vs dis)')

subplot(2,2,4)
hold on
plot(distancearray, Xprot(1,:))
plot(distancearray, Xprot(numtime,:))
hold off
xlabel('Distance (um)'); ylabel('Concentration (mM)')
title('A diffusing autoregulatory gene (con vs dis)')
legend('Protein first','Protein last')

%%%%%%%%%%%%%%%
%% 3
%%%%%%%%%%%%%%%

u=1; w=1; CXprot=1; CXrna=1; K=0.33; Dprot=1*10^(-4); Drna=1*10^(-4);
totaltime = 30;
totaldistance = 3;
deltat = 0.01;
deltad = 0.02;

numtime = totaltime/deltat;
numposition = totaldistance/deltad;
timearray = zeros(numtime,1);
distancearray = zeros(numposition,1);
Xrna = zeros(numtime, numposition);
Xprot = zeros(numtime, numposition);
Xrna(1,1) = 1;
Xprot(1,1) = 1;

for t = 1:numtime-1
    for x = 1:numposition
        if x == 1
            laplacian_xp = (Xprot(t,x+1)-2*Xprot(t,x))/(deltad^2);
            laplacian_xr = (Xrna(t,x+1)-2*Xrna(t,x))/(deltad^2);
        elseif x == numposition
            laplacian_xp = (Xprot(t,x-1)-2*Xprot(t,x))/(deltad^2);
            laplacian_xr = (Xrna(t,x-1)-2*Xrna(t,x))/(deltad^2);
        else
            laplacian_xp = (Xprot(t,x+1)+Xprot(t,x-1)-2*Xprot(t,x))/(deltad^2);
            laplacian_xr = (Xrna(t,x+1)+Xrna(t,x-1)-2*Xrna(t,x))/(deltad^2);
        end
        
        dXrnadt = u*(Xprot(t,x)^2)/(K^2+Xprot(t,x)^2)-CXrna*Xrna(t,x)+Drna*laplacian_xr;
        dXprotdt = w*Xrna(t,x)-CXprot*Xprot(t,x)+Dprot*laplacian_xp;
        Xrna(t+1,x) = Xrna(t,x) + dXrnadt*deltat;
        Xprot(t+1,x) = Xprot(t,x) + dXprotdt*deltat;
        
        if x ~= numposition
            distancearray(x+1) = distancearray(x) + deltad;
        end
    end
    timearray(t+1) = timearray(t) + deltat;    
end

figure(1)
subplot(2,2,1)
hold on
plot(timearray, Xrna(:,1))
plot(timearray, Xrna(:,2))
plot(timearray, Xrna(:,3))
hold off
xlabel('Time (s)'); ylabel('Concentration (mM)')
title('A diffusing autoregulatory gene (con vs time)')
legend('RNA first','RNA second','RNA third')

subplot(2,2,2)
hold on 
plot(timearray, Xprot(:,1))
plot(timearray, Xprot(:,2))
plot(timearray, Xprot(:,3))
hold off
xlabel('Time (s)'); ylabel('Concentration (mM)')
title('A diffusing autoregulatory gene (con vs time)')
legend('Protein first','Protein second','Protein third')

subplot(2,2,3)
hold on
plot(distancearray, Xrna(1,:))
plot(distancearray, Xrna(numtime,:))
hold off
xlabel('Distance (um)'); ylabel('Concentration (mM)')
legend('RNA first','RNA last')
title('A diffusing autoregulatory gene (con vs dis)')

subplot(2,2,4)
hold on
plot(distancearray, Xprot(1,:))
plot(distancearray, Xprot(numtime,:))
hold off
xlabel('Distance (um)'); ylabel('Concentration (mM)')
title('A diffusing autoregulatory gene (con vs dis)')
legend('Protein first','Protein last')

%%%%%%%%%%%%%%%
%% 4
%%%%%%%%%%%%%%%

u=1; w=1; CXprot=1; CXrna=1; K=0.33; Dprot=1*10^(-4); Drna=1*10^(-4);
totaltime = 50;
totaldistance = 3;
deltat = 0.01;
deltad = 0.02;

numtime = totaltime/deltat;
numposition = totaldistance/deltad;
timearray = zeros(numtime,1);
distancearray = zeros(numposition,1);
Xrna = zeros(numtime, numposition);
Xprot = zeros(numtime, numposition);
Xrna(1,1) = 1;
Xprot(1,1) = 1;

for t = 1:numtime-1
    for x = 1:numposition
        if x == 1
            laplacian_xp = (Xprot(t,x+1)-2*Xprot(t,x))/(deltad^2);
            laplacian_xr = (Xrna(t,x+1)-2*Xrna(t,x))/(deltad^2);
        elseif x == numposition
            laplacian_xp = (Xprot(t,x-1)-2*Xprot(t,x))/(deltad^2);
            laplacian_xr = (Xrna(t,x-1)-2*Xrna(t,x))/(deltad^2);
        else
            laplacian_xp = (Xprot(t,x+1)+Xprot(t,x-1)-2*Xprot(t,x))/(deltad^2);
            laplacian_xr = (Xrna(t,x+1)+Xrna(t,x-1)-2*Xrna(t,x))/(deltad^2);
        end
        
        dXrnadt = u*(Xprot(t,x)^2)/(K^2+Xprot(t,x)^2)-CXrna*Xrna(t,x)+Drna*laplacian_xr;
        dXprotdt = w*Xrna(t,x)-CXprot*Xprot(t,x)+Dprot*laplacian_xp;
        Xrna(t+1,x) = Xrna(t,x) + dXrnadt*deltat;
        Xprot(t+1,x) = Xprot(t,x) + dXprotdt*deltat;
        
        if x ~= numposition
            distancearray(x+1) = distancearray(x) + deltad;
        end
    end
    timearray(t+1) = timearray(t) + deltat;    
end

figure(1)
subplot(2,2,1)
hold on
plot(timearray, Xrna(:,1))
plot(timearray, Xrna(:,2))
plot(timearray, Xrna(:,3))
hold off
xlabel('Time (s)'); ylabel('Concentration (mM)')
title('A diffusing autoregulatory gene (con vs time)')
legend('RNA first','RNA second','RNA third')

subplot(2,2,2)
hold on 
plot(timearray, Xprot(:,1))
plot(timearray, Xprot(:,2))
plot(timearray, Xprot(:,3))
hold off
xlabel('Time (s)'); ylabel('Concentration (mM)')
title('A diffusing autoregulatory gene (con vs time)')
legend('Protein first','Protein second','Protein third')

subplot(2,2,3)
hold on
plot(distancearray, Xrna(1,:))
plot(distancearray, Xrna(numtime,:))
hold off
xlabel('Distance (um)'); ylabel('Concentration (mM)')
legend('RNA first','RNA last')
title('A diffusing autoregulatory gene (con vs dis)')

subplot(2,2,4)
hold on
plot(distancearray, Xprot(1,:))
plot(distancearray, Xprot(numtime,:))
hold off
xlabel('Distance (um)'); ylabel('Concentration (mM)')
title('A diffusing autoregulatory gene (con vs dis)')
legend('Protein first','Protein last')

%%%% 5

%%%% Part B
%%%%%%%%%%%%%%%
%% 1
%%%%%%%%%%%%%%%

u=1; w=1; CXprot=1; CXrna=1; K=0.33; Dprot=1*10^(-4); Drna=1*10^(-4);
totaltime = 5; 
totaldistance = 3;
deltat = 0.01;
deltad = 0.02;

numtime = totaltime/deltat;
numposition = totaldistance/deltad;
timearray = zeros(numtime,1);
distancearray = zeros(numposition,1);

Xrna = zeros(numtime, numposition);
Xprot = zeros(numtime, numposition);
Yrna = zeros(numtime, numposition);
Yprot = zeros(numtime, numposition);

Xrna(1,1) = 1; Xrna(1,numposition) = 0;
Xprot(1,1) = 1; Xprot(1,numposition) = 0;
Yrna(1,1) = 0; Yrna(1,numposition) = 1;
Yprot(1,1) = 0; Yprot(1,numposition) = 1;


for t = 1:numtime-1
    for x = 1:numposition
        if x == 1
            laplacian_xp = (Xprot(t,x+1)-2*Xprot(t,x))/(deltad^2);
            laplacian_xr = (Xrna(t,x+1)-2*Xrna(t,x))/(deltad^2);
            laplacian_yp = (Yprot(t,x+1)-2*Yprot(t,x))/(deltad^2);
            laplacian_yr = (Yrna(t,x+1)-2*Yrna(t,x))/(deltad^2);              
        elseif x == numposition
            laplacian_xp = (Xprot(t,x-1)-2*Xprot(t,x))/(deltad^2);
            laplacian_xr = (Xrna(t,x-1)-2*Xrna(t,x))/(deltad^2);
            laplacian_yp = (Yprot(t,x-1)-2*Yprot(t,x))/(deltad^2);
            laplacian_yr = (Yrna(t,x-1)-2*Yrna(t,x))/(deltad^2); 
        else
            laplacian_xp = (Xprot(t,x+1)+Xprot(t,x-1)-2*Xprot(t,x))/(deltad^2);
            laplacian_xr = (Xrna(t,x+1)+Xrna(t,x-1)-2*Xrna(t,x))/(deltad^2);
            laplacian_yp = (Yprot(t,x+1)+Yprot(t,x-1)-2*Yprot(t,x))/(deltad^2);
            laplacian_yr = (Yrna(t,x+1)+Yrna(t,x-1)-2*Yrna(t,x))/(deltad^2);              
        end
        
        dXrnadt = u*(1-(Yprot(t,x)^2)/(K^2+Yprot(t,x)^2))-CXrna*Xrna(t,x)+Drna*laplacian_xr;
        dXprotdt = w*Xrna(t,x)-CXprot*Xprot(t,x)+Dprot*laplacian_xp;
        dYrnadt = u*(1-(Xprot(t,x)^2)/(K^2+Xprot(t,x)^2))-CXrna*Yrna(t,x)+Drna*laplacian_yr;
        dYprotdt = w*Yrna(t,x)-CXprot*Yprot(t,x)+Dprot*laplacian_yp;
        
        Xrna(t+1,x) = Xrna(t,x) + dXrnadt*deltat;
        Xprot(t+1,x) = Xprot(t,x) + dXprotdt*deltat;
        Yrna(t+1,x) = Yrna(t,x) + dYrnadt*deltat;
        Yprot(t+1,x) = Yprot(t,x) + dYprotdt*deltat;
  
        if x ~= numposition
            distancearray(x+1) = distancearray(x) + deltad;
        end
    end
    timearray(t+1) = timearray(t) + deltat;    
end

%%%%%%%%%%%%%%%
%% 2
%%%%%%%%%%%%%%%


u=1; w=1; CXprot=1; CXrna=1; K=0.33; Dprot=1*10^(-4); Drna=1*10^(-4);
totaltime = 5; 
totaldistance = 3;
deltat = 0.01;
deltad = 0.02;

numtime = totaltime/deltat;
numposition = totaldistance/deltad;
timearray = zeros(numtime,1);
distancearray = zeros(numposition,1);

Xrna = zeros(numtime, numposition);
Xprot = zeros(numtime, numposition);
Yrna = zeros(numtime, numposition);
Yprot = zeros(numtime, numposition);

Xrna(1,1) = 1; Xrna(1,numposition) = 0;
Xprot(1,1) = 1; Xprot(1,numposition) = 0;
Yrna(1,1) = 0; Yrna(1,numposition) = 1;
Yprot(1,1) = 0; Yprot(1,numposition) = 1;


for t = 1:numtime-1
    for x = 1:numposition
        if x == 1
            laplacian_xp = (Xprot(t,x+1)-2*Xprot(t,x))/(deltad^2);
            laplacian_xr = (Xrna(t,x+1)-2*Xrna(t,x))/(deltad^2);
            laplacian_yp = (Yprot(t,x+1)-2*Yprot(t,x))/(deltad^2);
            laplacian_yr = (Yrna(t,x+1)-2*Yrna(t,x))/(deltad^2);              
        elseif x == numposition
            laplacian_xp = (Xprot(t,x-1)-2*Xprot(t,x))/(deltad^2);
            laplacian_xr = (Xrna(t,x-1)-2*Xrna(t,x))/(deltad^2);
            laplacian_yp = (Yprot(t,x-1)-2*Yprot(t,x))/(deltad^2);
            laplacian_yr = (Yrna(t,x-1)-2*Yrna(t,x))/(deltad^2); 
        else
            laplacian_xp = (Xprot(t,x+1)+Xprot(t,x-1)-2*Xprot(t,x))/(deltad^2);
            laplacian_xr = (Xrna(t,x+1)+Xrna(t,x-1)-2*Xrna(t,x))/(deltad^2);
            laplacian_yp = (Yprot(t,x+1)+Yprot(t,x-1)-2*Yprot(t,x))/(deltad^2);
            laplacian_yr = (Yrna(t,x+1)+Yrna(t,x-1)-2*Yrna(t,x))/(deltad^2);              
        end
        
        dXrnadt = u*(1-((Yprot(t,x)^2)/(K^2+Yprot(t,x)^2)))-CXrna*Xrna(t,x)+Drna*laplacian_xr;
        dXprotdt = w*Xrna(t,x)-CXprot*Xprot(t,x)+Dprot*laplacian_xp;
        dYrnadt = u*(1-((Xprot(t,x)^2)/(K^2+Xprot(t,x)^2)))-CXrna*Yrna(t,x)+Drna*laplacian_yr;
        dYprotdt = w*Yrna(t,x)-CXprot*Yprot(t,x)+Dprot*laplacian_yp;
        
        Xrna(t+1,x) = Xrna(t,x) + dXrnadt*deltat;
        Xprot(t+1,x) = Xprot(t,x) + dXprotdt*deltat;
        Yrna(t+1,x) = Yrna(t,x) + dYrnadt*deltat;
        Yprot(t+1,x) = Yprot(t,x) + dYprotdt*deltat;
  
        if x ~= numposition
            distancearray(x+1) = distancearray(x) + deltad;
        end
    end
    timearray(t+1) = timearray(t) + deltat;    
end

figure(1)
subplot(4,2,1)
hold on
plot(timearray, Xrna(:,1))
plot(timearray, Xrna(:,numposition/2))
plot(timearray, Xrna(:,numposition))
hold off
xlabel('Time (s)'); ylabel('Concentration (mM)')
title('A diffusing pair of mutually-inhibiting genes (con vs time)')
legend('X RNA first','X RNA middle','X RNA last')

subplot(4,2,2)
hold on 
plot(timearray, Xprot(:,1))
plot(timearray, Xprot(:,numposition/2))
plot(timearray, Xprot(:,numposition))
hold off
xlabel('Time (s)'); ylabel('Concentration (mM)')
title('A diffusing pair of mutually-inhibiting genes (con vs time)')
legend('X Protein first','X Protein middle','X Protein last')

subplot(4,2,3)
hold on
plot(timearray, Yrna(:,1))
plot(timearray, Yrna(:,numposition/2))
plot(timearray, Yrna(:,numposition))
hold off
xlabel('Time (s)'); ylabel('Concentration (mM)')
title('A diffusing pair of mutually-inhibiting genes (con vs time)')
legend('Y RNA first','Y RNA middle','Y RNA last')

subplot(4,2,4)
hold on 
plot(timearray, Yprot(:,1))
plot(timearray, Yprot(:,numposition/2))
plot(timearray, Yprot(:,numposition))
hold off
xlabel('Time (s)'); ylabel('Concentration (mM)')
title('A diffusing pair of mutually-inhibiting genes (con vs time)')
legend('Y Protein first','Y Protein middle','Y Protein last')

subplot(4,2,5)
hold on
plot(distancearray, Xrna(1,:))
plot(distancearray, Xrna(numtime/2,:))
plot(distancearray, Xrna(numtime,:))
hold off
xlabel('Distance (um)'); ylabel('Concentration (mM)')
legend('X RNA first','X RNA middle','X RNA last')
title('A diffusing pair of mutually-inhibiting genes (con vs dis)')

subplot(4,2,6)
hold on
plot(distancearray, Xprot(1,:))
plot(distancearray, Xprot(numtime/2,:))
plot(distancearray, Xprot(numtime,:))
hold off
xlabel('Distance (um)'); ylabel('Concentration (mM)')
title('A diffusing pair of mutually-inhibiting genes(con vs dis)')
legend('X Protein first','X Protein middle','X Protein last')


subplot(4,2,7)
hold on
plot(distancearray, Yrna(1,:))
plot(distancearray, Yrna(numtime/2,:))
plot(distancearray, Yrna(numtime,:))

hold off
xlabel('Distance (um)'); ylabel('Concentration (mM)')
legend('Y RNA first','Y RNA middle','Y RNA last')
title('A diffusing pair of mutually-inhibiting genes (con vs dis)')

subplot(4,2,8)
hold on
plot(distancearray, Xprot(1,:))
plot(distancearray, Xprot(numtime/2,:))
plot(distancearray, Xprot(numtime,:))

hold off
xlabel('Distance (um)'); ylabel('Concentration (mM)')
title('A diffusing pair of mutually-inhibiting genes (con vs dis)')
legend('Y Protein first','Y Protein middle','Y Protein last')

%%%%%%%%%%%%%%%
%% 3
%%%%%%%%%%%%%%%

u=1; w=1; CXprot=1; CXrna=1; K=0.33; Dprot=1*10^(-4); Drna=1*10^(-4);
totaltime = 50;% I changed the value of totaltime to 100,200,400 to get each graphs
totaldistance = 3;
deltat = 0.01;
deltad = 0.02;
numtime = totaltime/deltat;
numposition = totaldistance/deltad;
timearray = zeros(numtime,1);
distancearray = zeros(numposition,1);
Xrna = zeros(numtime, numposition);
Xprot = zeros(numtime, numposition);
Yrna = zeros(numtime, numposition);
Yprot = zeros(numtime, numposition);
Xrna(1,1) = 1; Xrna(1,numposition) = 0;
Xprot(1,1) = 1; Xprot(1,numposition) = 0;
Yrna(1,1) = 0; Yrna(1,numposition) = 1;
Yprot(1,1) = 0; Yprot(1,numposition) = 1;
for t = 1:numtime-1
   for x = 1:numposition
       if x == 1
           laplacian_xp = (Xprot(t,x+1)-2*Xprot(t,x))/(deltad^2);
           laplacian_xr = (Xrna(t,x+1)-2*Xrna(t,x))/(deltad^2);
           laplacian_yp = (Yprot(t,x+1)-2*Yprot(t,x))/(deltad^2);
           laplacian_yr = (Yrna(t,x+1)-2*Yrna(t,x))/(deltad^2);             
       elseif x == numposition
           laplacian_xp = (Xprot(t,x-1)-2*Xprot(t,x))/(deltad^2);
           laplacian_xr = (Xrna(t,x-1)-2*Xrna(t,x))/(deltad^2);
           laplacian_yp = (Yprot(t,x-1)-2*Yprot(t,x))/(deltad^2);
           laplacian_yr = (Yrna(t,x-1)-2*Yrna(t,x))/(deltad^2);
       else
           laplacian_xp = (Xprot(t,x+1)+Xprot(t,x-1)-2*Xprot(t,x))/(deltad^2);
           laplacian_xr = (Xrna(t,x+1)+Xrna(t,x-1)-2*Xrna(t,x))/(deltad^2);
           laplacian_yp = (Yprot(t,x+1)+Yprot(t,x-1)-2*Yprot(t,x))/(deltad^2);
           laplacian_yr = (Yrna(t,x+1)+Yrna(t,x-1)-2*Yrna(t,x))/(deltad^2);             
       end
      
       dXrnadt = u*(1-((Yprot(t,x)^2)/(K^2+Yprot(t,x)^2)))-CXrna*Xrna(t,x)+Drna*laplacian_xr;
       dXprotdt = w*Xrna(t,x)-CXprot*Xprot(t,x)+Dprot*laplacian_xp;
       dYrnadt = u*(1-((Xprot(t,x)^2)/(K^2+Xprot(t,x)^2)))-CXrna*Yrna(t,x)+Drna*laplacian_yr;
       dYprotdt = w*Yrna(t,x)-CXprot*Yprot(t,x)+Dprot*laplacian_yp;
      
       Xrna(t+1,x) = Xrna(t,x) + dXrnadt*deltat;
       Xprot(t+1,x) = Xprot(t,x) + dXprotdt*deltat;
       Yrna(t+1,x) = Yrna(t,x) + dYrnadt*deltat;
       Yprot(t+1,x) = Yprot(t,x) + dYprotdt*deltat;
        if x ~= numposition
           distancearray(x+1) = distancearray(x) + deltad;
       end
   end
   timearray(t+1) = timearray(t) + deltat;   
end
figure(1)
subplot(4,2,1)
hold on
plot(timearray, Xrna(:,1))
plot(timearray, Xrna(:,numposition/2))
plot(timearray, Xrna(:,numposition))
hold off
xlabel('Time (s)'); ylabel('Concentration (mM)')
title('A diffusing pair of mutually-inhibiting genes (con vs time)')
legend('X RNA first','X RNA middle','X RNA last')
subplot(4,2,2)
hold on
plot(timearray, Xprot(:,1))
plot(timearray, Xprot(:,numposition/2))
plot(timearray, Xprot(:,numposition))
hold off
xlabel('Time (s)'); ylabel('Concentration (mM)')
title('A diffusing pair of mutually-inhibiting genes (con vs time)')
legend('X Protein first','X Protein middle','X Protein last')
subplot(4,2,3)
hold on
plot(timearray, Yrna(:,1))
plot(timearray, Yrna(:,numposition/2))
plot(timearray, Yrna(:,numposition))
hold off
xlabel('Time (s)'); ylabel('Concentration (mM)')
title('A diffusing pair of mutually-inhibiting genes (con vs time)')
legend('Y RNA first','Y RNA middle','Y RNA last')
subplot(4,2,4)
hold on
plot(timearray, Yprot(:,1))
plot(timearray, Yprot(:,numposition/2))
plot(timearray, Yprot(:,numposition))
hold off
xlabel('Time (s)'); ylabel('Concentration (mM)')
title('A diffusing pair of mutually-inhibiting genes (con vs time)')
legend('Y Protein first','Y Protein middle','Y Protein last')
subplot(4,2,5)
hold on
plot(distancearray, Xrna(1,:))
plot(distancearray, Xrna(numtime/2,:))
plot(distancearray, Xrna(numtime,:))
hold off
xlabel('Distance (um)'); ylabel('Concentration (mM)')
legend('X RNA first','X RNA middle','X RNA last')
title('A diffusing pair of mutually-inhibiting genes (con vs dis)')
subplot(4,2,6)
hold on
plot(distancearray, Xprot(1,:))
plot(distancearray, Xprot(numtime/2,:))
plot(distancearray, Xprot(numtime,:))
hold off
xlabel('Distance (um)'); ylabel('Concentration (mM)')
title('A diffusing pair of mutually-inhibiting genes(con vs dis)')
legend('X Protein first','X Protein middle','X Protein last')
subplot(4,2,7)
hold on
plot(distancearray, Yrna(1,:))
plot(distancearray, Yrna(numtime/2,:))
plot(distancearray, Yrna(numtime,:))
hold off
xlabel('Distance (um)'); ylabel('Concentration (mM)')
legend('Y RNA first','Y RNA middle','Y RNA last')
title('A diffusing pair of mutually-inhibiting genes (con vs dis)')
subplot(4,2,8)
hold on
plot(distancearray, Xprot(1,:))
plot(distancearray, Xprot(numtime/2,:))
plot(distancearray, Xprot(numtime,:))
hold off
xlabel('Distance (um)'); ylabel('Concentration (mM)')
title('A diffusing pair of mutually-inhibiting genes (con vs dis)')
legend('Y Protein first','Y Protein middle','Y Protein last')

%%%%%%%%%%%%%%%
%% 4
%%%%%%%%%%%%%%%

u=1; w=1; CXprot=1; CXrna=1; K=0.33; Dprot=1*10^(-4); Drna=1*10^(-4);
totaltime = 50;% I changed the value of totaltime to 100,200,400 to get each graphs
totaldistance = 3;
deltat = 0.01;
deltad = 0.02;
numtime = totaltime/deltat;
numposition = totaldistance/deltad;
timearray = zeros(numtime,1);
distancearray = zeros(numposition,1);
Xrna = zeros(numtime, numposition);
Xprot = zeros(numtime, numposition);
Yrna = zeros(numtime, numposition);
Yprot = zeros(numtime, numposition);

Xrna(1,1) = 1; Xrna(1,numposition) = 1;
Xprot(1,1) = 1; Xprot(1,numposition) = 1;
Yrna(1,1) = 0; Yrna(1,numposition) = 0;
Yprot(1,1) = 0; Yprot(1,numposition) = 0;

Xrna(1,numposition/2) = 0;
Xprot(1,numposition/2) = 0;
Yrna(1,numposition/2) = 1;
Yprot(1,numposition/2) = 1;
  
for t = 1:numtime-1
  for x = 1:numposition
      if x == 1
          laplacian_xp = (Xprot(t,x+1)-2*Xprot(t,x))/(deltad^2);
          laplacian_xr = (Xrna(t,x+1)-2*Xrna(t,x))/(deltad^2);
          laplacian_yp = (Yprot(t,x+1)-2*Yprot(t,x))/(deltad^2);
          laplacian_yr = (Yrna(t,x+1)-2*Yrna(t,x))/(deltad^2);            
      elseif x == numposition
          laplacian_xp = (Xprot(t,x-1)-2*Xprot(t,x))/(deltad^2);
          laplacian_xr = (Xrna(t,x-1)-2*Xrna(t,x))/(deltad^2);
          laplacian_yp = (Yprot(t,x-1)-2*Yprot(t,x))/(deltad^2);
          laplacian_yr = (Yrna(t,x-1)-2*Yrna(t,x))/(deltad^2);
      else
          laplacian_xp = (Xprot(t,x+1)+Xprot(t,x-1)-2*Xprot(t,x))/(deltad^2);
          laplacian_xr = (Xrna(t,x+1)+Xrna(t,x-1)-2*Xrna(t,x))/(deltad^2);
          laplacian_yp = (Yprot(t,x+1)+Yprot(t,x-1)-2*Yprot(t,x))/(deltad^2);
          laplacian_yr = (Yrna(t,x+1)+Yrna(t,x-1)-2*Yrna(t,x))/(deltad^2);            
      end
    
      dXrnadt = u*(1-((Yprot(t,x)^2)/(K^2+Yprot(t,x)^2)))-CXrna*Xrna(t,x)+Drna*laplacian_xr;
      dXprotdt = w*Xrna(t,x)-CXprot*Xprot(t,x)+Dprot*laplacian_xp;
      dYrnadt = u*(1-((Xprot(t,x)^2)/(K^2+Xprot(t,x)^2)))-CXrna*Yrna(t,x)+Drna*laplacian_yr;
      dYprotdt = w*Yrna(t,x)-CXprot*Yprot(t,x)+Dprot*laplacian_yp;
    
      Xrna(t+1,x) = Xrna(t,x) + dXrnadt*deltat;
      Xprot(t+1,x) = Xprot(t,x) + dXprotdt*deltat;
      Yrna(t+1,x) = Yrna(t,x) + dYrnadt*deltat;
      Yprot(t+1,x) = Yprot(t,x) + dYprotdt*deltat;
       if x ~= numposition
          distancearray(x+1) = distancearray(x) + deltad;
      end
  end
  timearray(t+1) = timearray(t) + deltat;  
end
figure(1)
subplot(4,2,1)
hold on
plot(timearray, Xrna(:,1))
plot(timearray, Xrna(:,numposition/2))
plot(timearray, Xrna(:,numposition))
hold off
xlabel('Time (s)'); ylabel('Concentration (mM)')
title('A diffusing pair of mutually-inhibiting genes (con vs time)')
legend('X RNA first','X RNA middle','X RNA last')
subplot(4,2,2)
hold on
plot(timearray, Xprot(:,1))
plot(timearray, Xprot(:,numposition/2))
plot(timearray, Xprot(:,numposition))
hold off
xlabel('Time (s)'); ylabel('Concentration (mM)')
title('A diffusing pair of mutually-inhibiting genes (con vs time)')
legend('X Protein first','X Protein middle','X Protein last')
subplot(4,2,3)
hold on
plot(timearray, Yrna(:,1))
plot(timearray, Yrna(:,numposition/2))
plot(timearray, Yrna(:,numposition))
hold off
xlabel('Time (s)'); ylabel('Concentration (mM)')
title('A diffusing pair of mutually-inhibiting genes (con vs time)')
legend('Y RNA first','Y RNA middle','Y RNA last')
subplot(4,2,4)
hold on
plot(timearray, Yprot(:,1))
plot(timearray, Yprot(:,numposition/2))
plot(timearray, Yprot(:,numposition))
hold off
xlabel('Time (s)'); ylabel('Concentration (mM)')
title('A diffusing pair of mutually-inhibiting genes (con vs time)')
legend('Y Protein first','Y Protein middle','Y Protein last')
subplot(4,2,5)
hold on
plot(distancearray, Xrna(1,:))
plot(distancearray, Xrna(numtime/2,:))
plot(distancearray, Xrna(numtime,:))
hold off
xlabel('Distance (um)'); ylabel('Concentration (mM)')
legend('X RNA first','X RNA middle','X RNA last')
title('A diffusing pair of mutually-inhibiting genes (con vs dis)')
subplot(4,2,6)
hold on
plot(distancearray, Xprot(1,:))
plot(distancearray, Xprot(numtime/2,:))
plot(distancearray, Xprot(numtime,:))
hold off
xlabel('Distance (um)'); ylabel('Concentration (mM)')
title('A diffusing pair of mutually-inhibiting genes(con vs dis)')
legend('X Protein first','X Protein middle','X Protein last')
subplot(4,2,7)
hold on
plot(distancearray, Yrna(1,:))
plot(distancearray, Yrna(numtime/2,:))
plot(distancearray, Yrna(numtime,:))

hold off
xlabel('Distance (um)'); ylabel('Concentration (mM)')
legend('Y RNA first','Y RNA middle','Y RNA last')
title('A diffusing pair of mutually-inhibiting genes (con vs dis)')
subplot(4,2,8)
hold on
plot(distancearray, Xprot(1,:))
plot(distancearray, Xprot(numtime/2,:))
plot(distancearray, Xprot(numtime,:))
hold off
xlabel('Distance (um)'); ylabel('Concentration (mM)')
title('A diffusing pair of mutually-inhibiting genes (con vs dis)')
legend('Y Protein first','Y Protein middle','Y Protein last')





