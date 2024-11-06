%%%% Part A
%%%%%%%%%%%%%%%
%% 1 
%%%%%%%%%%%%%%%
% K = K(1/2)


%case 1 with all initial concentrations set to 0
Xcirna=1.2; Xciprot=1.2;
Xcrorna=0.8; Xcroprot=0.8;
Wci=50;Uci=50;Wcro=50;Ucro=50;
Kci=10; Kcro=10;

totaltime = 15;
stepsize = 0.01;
numsteps = totaltime/stepsize;

Ciprot = zeros(numsteps,1);
Cirna = zeros(numsteps,1);
Croprot = zeros(numsteps,1);
Crorna = zeros(numsteps,1);
t = zeros(numsteps,1);

Ciprot(1)=0; Cirna(1)=0;
Croprot(1)=0; Crorna(1)=0; t(1)=0;

for (i=1:numsteps-1)
    dCiprot_dt = Wci*Cirna(i) - Xciprot*Ciprot(i);
    dCirna_dt = Uci*(1-(Croprot(i)^2)/(Kcro^2+Croprot(i)^2))-Xcirna*Cirna(i);
    dCroprot_dt = Wcro*Crorna(i) - Xcroprot*Croprot(i);
    dCrorna_dt = Ucro*(1-(Ciprot(i)^2)/(Kci^2+Ciprot(i)^2))-Xcrorna*Crorna(i);     
    Ciprot(i+1) = Ciprot(i) + dCiprot_dt*stepsize;
    Cirna(i+1) = Cirna(i) + dCirna_dt*stepsize;
    Croprot(i+1) = Croprot(i) + dCroprot_dt*stepsize;
    Crorna(i+1) = Crorna(i) + dCrorna_dt*stepsize;        
    t(i+1) = t(i) + stepsize;
end
hold on
plot(t,Ciprot);plot(t,Cirna);plot(t,Croprot);plot(t,Crorna);
hold off
xlabel('Time (s)'); ylabel('Number of Molecules')
title('Forward Euler algorithm about Ci & Cro')
legend('Ciprot','Cirna','Croprot','Crorna')

%%%%%%%%%%%%%%%
%% case 2 with starting only with 20 molecules of cro RNA present
%%%%%%%%%%%%%%%
Xcirna=1.2; Xciprot=1.2;
Xcrorna=0.8; Xcroprot=0.8;
Wci=50;Uci=50;Wcro=50;Ucro=50;
Kci=10; Kcro=10;
totaltime = 15;
stepsize = 0.01;
numsteps = totaltime/stepsize;

Ciprot = zeros(numsteps,1);
Cirna = zeros(numsteps,1);
Croprot = zeros(numsteps,1);
Crorna = zeros(numsteps,1);
t = zeros(numsteps,1);

Ciprot(1)=0; Cirna(1)=0;
Croprot(1)=0; Crorna(1)=20; t(1)=0;

for (i=1:numsteps-1)
    dCiprot_dt = Wci*Cirna(i) - Xciprot*Ciprot(i);
    dCirna_dt = Uci*(1-(Croprot(i)^2)/(Kcro^2+Croprot(i)^2))-Xcirna*Cirna(i);
    dCroprot_dt = Wcro*Crorna(i) - Xcroprot*Croprot(i);
    dCrorna_dt = Ucro*(1-(Ciprot(i)^2)/(Kci^2+Ciprot(i)^2))-Xcrorna*Crorna(i);     
    Ciprot(i+1) = Ciprot(i) + dCiprot_dt*stepsize;
    Cirna(i+1) = Cirna(i) + dCirna_dt*stepsize;
    Croprot(i+1) = Croprot(i) + dCroprot_dt*stepsize;
    Crorna(i+1) = Crorna(i) + dCrorna_dt*stepsize;        
    t(i+1) = t(i) + stepsize;
end
hold on
plot(t,Ciprot);plot(t,Cirna);plot(t,Croprot);plot(t,Crorna);
hold off
xlabel('Time (s)'); ylabel('Number of Molecules')
title('Forward Euler algorithm about Ci & Cro')
legend('Ciprot','Cirna','Croprot','Crorna')

%%%%%%%%%%%%%%%
%% case 3 with starting only with 50 molecules of cI RNA present
%%%%%%%%%%%%%%%
Xcirna=1.2; Xciprot=1.2;
Xcrorna=0.8; Xcroprot=0.8;
Wci=50;Uci=50;Wcro=50;Ucro=50;
Kci=10; Kcro=10;
totaltime = 15;
stepsize = 0.01;
numsteps = totaltime/stepsize;

Ciprot = zeros(numsteps,1);
Cirna = zeros(numsteps,1);
Croprot = zeros(numsteps,1);
Crorna = zeros(numsteps,1);
t = zeros(numsteps,1);

Ciprot(1)=0; Cirna(1)=50;
Croprot(1)=0; Crorna(1)=0; t(1)=0;

for (i=1:numsteps-1)
    dCiprot_dt = Wci*Cirna(i) - Xciprot*Ciprot(i);
    dCirna_dt = Uci*(1-(Croprot(i)^2)/(Kcro^2+Croprot(i)^2))-Xcirna*Cirna(i);
    dCroprot_dt = Wcro*Crorna(i) - Xcroprot*Croprot(i);
    dCrorna_dt = Ucro*(1-(Ciprot(i)^2)/(Kci^2+Ciprot(i)^2))-Xcrorna*Crorna(i);     
    Ciprot(i+1) = Ciprot(i) + dCiprot_dt*stepsize;
    Cirna(i+1) = Cirna(i) + dCirna_dt*stepsize;
    Croprot(i+1) = Croprot(i) + dCroprot_dt*stepsize;
    Crorna(i+1) = Crorna(i) + dCrorna_dt*stepsize;        
    t(i+1) = t(i) + stepsize;
end
hold on
plot(t,Ciprot);plot(t,Cirna);plot(t,Croprot);plot(t,Crorna);
hold off
xlabel('Time (s)'); ylabel('Number of Molecules')
title('Forward Euler algorithm about Ci & Cro')
legend('Ciprot','Cirna','Croprot','Crorna')

%%%% 2
%%%%%%%%%%%%%%%
%% 2-(1)
%%%%%%%%%%%%%%%

Xcirna=1.2; Xciprot=1.2;
Xcrorna=0.8; Xcroprot=0.8;
Wci=50;Uci=50;Wcro=50;Ucro=50;
Kci=10; Kcro=10;
totaltime = 15;
stepsize = 0.01;
numsteps = totaltime/stepsize;

Ciprot = zeros(numsteps,1);
Cirna = zeros(numsteps,1);
Croprot = zeros(numsteps,1);
Crorna = zeros(numsteps,1);
t = zeros(numsteps,1);

Ciprot(1)=0; Cirna(1)=0;
Croprot(1)=0; Crorna(1)=0; t(1)=0;

for (k=0:1:20)
   Crorna(1)=k;
   for (j=0:1:20)
       Cirna(1)=j;
       for (i=1:numsteps-1)
           dCiprot_dt = Wci*Cirna(i) - Xciprot*Ciprot(i);
           dCirna_dt = Uci*(1-(Croprot(i)^2)/(Kcro^2+Croprot(i)^2))-Xcirna*Cirna(i);
           dCroprot_dt = Wcro*Crorna(i) - Xcroprot*Croprot(i);
           dCrorna_dt = Ucro*(1-(Ciprot(i)^2)/(Kci^2+Ciprot(i)^2))-Xcrorna*Crorna(i);     
           Ciprot(i+1) = Ciprot(i) + dCiprot_dt*stepsize;
           Cirna(i+1) = Cirna(i) + dCirna_dt*stepsize;
           Croprot(i+1) = Croprot(i) + dCroprot_dt*stepsize;
           Crorna(i+1) = Crorna(i) + dCrorna_dt*stepsize;        
       end
       hold on     
       plot(Ciprot,Croprot)
       plot(Ciprot(end),Croprot(end),'rx')
       xlabel('[CI prot]'); ylabel('[Cro prot]')
       title('Simulation of Ci VS Cro')
   end
end
hold off

%%%%%%%%%%%%%%%
%% 2-(2)
%%%%%%%%%%%%%%%

Xcirna=1.2; Xciprot=1.2;
Xcrorna=0.8; Xcroprot=0.8;
Wci=50;Uci=50;Wcro=50;Ucro=50;
Kci=10; Kcro=10;
totaltime = 15;
stepsize = 0.01;
numsteps = totaltime/stepsize;

Ciprot = zeros(numsteps,1);
Cirna = zeros(numsteps,1);
Croprot = zeros(numsteps,1);
Crorna = zeros(numsteps,1);
t = zeros(numsteps,1);

Ciprot(1)=0; Cirna(1)=0;
Croprot(1)=0; Crorna(1)=0; t(1)=0;

for (k=0:500:2000)
   Crorna(1)=k;
   for (j=0:500:2000)
       Cirna(1)=j;
       for (i=1:numsteps-1)
           dCiprot_dt = Wci*Cirna(i) - Xciprot*Ciprot(i);
           dCirna_dt = Uci*(1-(Croprot(i)^2)/(Kcro^2+Croprot(i)^2))-Xcirna*Cirna(i);
           dCroprot_dt = Wcro*Crorna(i) - Xcroprot*Croprot(i);
           dCrorna_dt = Ucro*(1-(Ciprot(i)^2)/(Kci^2+Ciprot(i)^2))-Xcrorna*Crorna(i);     
           Ciprot(i+1) = Ciprot(i) + dCiprot_dt*stepsize;
           Cirna(i+1) = Cirna(i) + dCirna_dt*stepsize;
           Croprot(i+1) = Croprot(i) + dCroprot_dt*stepsize;
           Crorna(i+1) = Crorna(i) + dCrorna_dt*stepsize;        
       end        
       hold on
       plot(Ciprot,Croprot)
       plot(Ciprot(end),Croprot(end),'rx')
       xlabel('[CI prot]'); ylabel('[Cro prot]')
       title('Simulation of Ci VS Cro')
   end
end
hold off

%%%%%%%%%%%%%%%
%% 4
%%%%%%%%%%%%%%%

Xcirna=1.2; Xciprot=1.2;
Xcrorna=0.8; Xcroprot=0.8;
Wci=50;Uci=50;Wcro=50;Ucro=50;
Kci=10; Kcro=10;
totaltime = 30;
stepsize = 0.01;
numsteps = totaltime/stepsize;
Ciprot = zeros(numsteps,1);
Cirna = zeros(numsteps,1);
Croprot = zeros(numsteps,1);
Crorna = zeros(numsteps,1);
t = zeros(numsteps,1);
Ciprot(1)=0; Cirna(1)=50;
Croprot(1)=0; Crorna(1)=0; t(1)=0;
for (i=1:numsteps-1)
   if i>(numsteps/2)
       Xciprot=10;
   end
   dCiprot_dt = Wci*Cirna(i) - Xciprot*Ciprot(i);
   dCirna_dt = Uci*(1-(Croprot(i)^2)/(Kcro^2+Croprot(i)^2))-Xcirna*Cirna(i);
   dCroprot_dt = Wcro*Crorna(i) - Xcroprot*Croprot(i);
   dCrorna_dt = Ucro*(1-(Ciprot(i)^2)/(Kci^2+Ciprot(i)^2))-Xcrorna*Crorna(i);    
   Ciprot(i+1) = Ciprot(i) + dCiprot_dt*stepsize;
   Cirna(i+1) = Cirna(i) + dCirna_dt*stepsize;
   Croprot(i+1) = Croprot(i) + dCroprot_dt*stepsize;
   Crorna(i+1) = Crorna(i) + dCrorna_dt*stepsize;       
   t(i+1) = t(i) + stepsize;
end
hold on
plot(t,Ciprot);plot(t,Cirna);plot(t,Croprot);plot(t,Crorna);
hold off
xlabel('Time (s)'); ylabel('Number of Molecules')
title('Forward Euler algorithm about Ci & Cro')
legend('Ciprot','Cirna','Croprot','Crorna')

%%%% Part B
%%%%%%%%%%%%%%%
%% 1.
%%%%%%%%%%%%%%%

Xcirna=1.2; Xciprot=1.2;
Xcrorna=0.8; Xcroprot=0.8;
Wci=50;Uci=50;Wcro=50;Ucro=50;
Kci=10; Kcro=10;

initial_num = 0;
final_num = 50000;
num_reactions=8;

Ciprot = zeros(final_num,1);
Cirna = zeros(final_num,1);
Croprot = zeros(final_num,1);
Crorna = zeros(final_num,1);
time = zeros(final_num,1);

Ciprot(1)=0; Cirna(1)=0;
Croprot(1)=0; Crorna(1)=0;
v = zeros(1,num_reactions);

for i=1:(final_num-1)
    v(1)=Wci*Cirna(i);
    v(2)=Xciprot*Ciprot(i);
    v(3)=Uci*(1-(Croprot(i)^2)/(Kcro^2+Croprot(i)^2));
    v(4)=Xcirna*Cirna(i);
    v(5)=Wcro*Crorna(i);
    v(6)=Xcroprot*Croprot(i);
    v(7)=Ucro*(1-(Ciprot(i)^2)/(Kci^2+Ciprot(i)^2));
    v(8)=Xcrorna*Crorna(i);
    alpha0=sum(v);
    y = rand();
    y_mod = (alpha0)*y;
    tau = -log(y)/(alpha0);
    time(i+1) = time(i)+tau;    
    next_v = 0;
    for j=1:num_reactions
        if (y_mod < sum(v(1:j)))
            next_v = j;
            break;
        end
    end 
    if (next_v == 1)
        Ciprot(i+1:final_num) = Ciprot(i) + 1;
    elseif (next_v == 2)
        Ciprot(i+1:final_num) = Ciprot(i) - 1;
    elseif (next_v == 3)
        Cirna(i+1:final_num) = Cirna(i) + 1;
    elseif (next_v == 4)
        Cirna(i+1:final_num) = Cirna(i) - 1;
    elseif (next_v == 5)
        Croprot(i+1:final_num) = Croprot(i) + 1;
    elseif (next_v == 6)
        Croprot(i+1:final_num) = Croprot(i) - 1;
    elseif (next_v == 7)
        Crorna(i+1:final_num) = Crorna(i) + 1;
    elseif (next_v == 8)
        Crorna(i+1:final_num) = Crorna(i) - 1;
    else
        disp('something went wrong');
    end
end

hold on
plot(time,Ciprot,'red');plot(time,Cirna);
plot(time,Croprot);plot(time,Crorna);
hold off
xlabel('Time (s)'); ylabel('Number of Molecules')
title('Stochastic model')
legend('Ciprot','Cirna','Croprot','Crorna')

%%%%%%%%%%%%%%%
%% 2 -1
%%%%%%%%%%%%%%%

Xcirna=1.2; Xciprot=1.2;
Xcrorna=0.8; Xcroprot=0.8;
Wci=50;Uci=50;Wcro=50;Ucro=50;
Kci=10; Kcro=10;

initial_num = 0;
final_num = 50000;
num_reactions=8;

Ciprot = zeros(final_num,1);
Cirna = zeros(final_num,1);
Croprot = zeros(final_num,1);
Crorna = zeros(final_num,1);
time = zeros(final_num,1);

Ciprot(1)=0; Cirna(1)=0;
Croprot(1)=0; Crorna(1)=0;
v = zeros(1,num_reactions);
n=0;
for k=1:20
    for i=1:(final_num-1)
        v(1)=Wci*Cirna(i);
        v(2)=Xciprot*Ciprot(i);
        v(3)=Uci*(1-(Croprot(i)^2)/(Kcro^2+Croprot(i)^2));
        v(4)=Xcirna*Cirna(i);
        v(5)=Wcro*Crorna(i);
        v(6)=Xcroprot*Croprot(i);
        v(7)=Ucro*(1-(Ciprot(i)^2)/(Kci^2+Ciprot(i)^2));
        v(8)=Xcrorna*Crorna(i);
        alpha0=sum(v);
        y = rand();
        y_mod = (alpha0)*y;
        tau = -log(y)/(alpha0);
        time(i+1) = time(i)+tau;    
        next_v = 0;
        for j=1:num_reactions
            if (y_mod < sum(v(1:j)))
                next_v = j;
                break;
            end
        end 
        if (next_v == 1)
            Ciprot(i+1:final_num) = Ciprot(i) + 1;
        elseif (next_v == 2)
            Ciprot(i+1:final_num) = Ciprot(i) - 1;
        elseif (next_v == 3)
            Cirna(i+1:final_num) = Cirna(i) + 1;
        elseif (next_v == 4)
            Cirna(i+1:final_num) = Cirna(i) - 1;
        elseif (next_v == 5)
            Croprot(i+1:final_num) = Croprot(i) + 1;
        elseif (next_v == 6)
            Croprot(i+1:final_num) = Croprot(i) - 1;
        elseif (next_v == 7)
            Crorna(i+1:final_num) = Crorna(i) + 1;
        elseif (next_v == 8)
            Crorna(i+1:final_num) = Crorna(i) - 1;
        else
            disp('something went wrong');
        end
    end
    n=n+1;
    subplot(4,5,n);
    hold on
    plot(time,Ciprot,'red');plot(time,Cirna);
    plot(time,Croprot,'blue');plot(time,Crorna);
    hold off
    xlabel('Time (s)'); ylabel('Number of Molecules')
    title('Stochastic model')
    legend('Ciprot','Cirna','Croprot','Crorna')
    
end

%%%%%%%%%%%%%%%
%% 2-2
%%%%%%%%%%%%%%%

Xcirna=1.2; Xciprot=1.2;
Xcrorna=0.8; Xcroprot=0.8;
Wci=50;Uci=50;Wcro=50;Ucro=50;
Kci=10; Kcro=10;

initial_num = 0;
final_num = 50000;
num_reactions=8;

Ciprot = zeros(final_num,1);
Cirna = zeros(final_num,1);
Croprot = zeros(final_num,1);
Crorna = zeros(final_num,1);
time = zeros(final_num,1);

Ciprot(1)=0; Cirna(1)=0;
Croprot(1)=0; Crorna(1)=0;
v = zeros(1,num_reactions);
n=0;
for k=1:20
    for i=1:(final_num-1)
        v(1)=Wci*Cirna(i);
        v(2)=Xciprot*Ciprot(i);
        v(3)=Uci*(1-(Croprot(i)^2)/(Kcro^2+Croprot(i)^2));
        v(4)=Xcirna*Cirna(i);
        v(5)=Wcro*Crorna(i);
        v(6)=Xcroprot*Croprot(i);
        v(7)=Ucro*(1-(Ciprot(i)^2)/(Kci^2+Ciprot(i)^2));
        v(8)=Xcrorna*Crorna(i);
        alpha0=sum(v);
        y = rand();
        y_mod = (alpha0)*y;
        tau = -log(y)/(alpha0);
        time(i+1) = time(i)+tau;    
        next_v = 0;
        for j=1:num_reactions
            if (y_mod < sum(v(1:j)))
                next_v = j;
                break;
            end
        end
        if (next_v == 1)
            Ciprot(i+1:final_num) = Ciprot(i) + 1;
        elseif (next_v == 2)
            Ciprot(i+1:final_num) = Ciprot(i) - 1;
        elseif (next_v == 3)
            Cirna(i+1:final_num) = Cirna(i) + 1;
        elseif (next_v == 4)
            Cirna(i+1:final_num) = Cirna(i) - 1;
        elseif (next_v == 5)
            Croprot(i+1:final_num) = Croprot(i) + 1;
        elseif (next_v == 6)
            Croprot(i+1:final_num) = Croprot(i) - 1;
        elseif (next_v == 7)
            Crorna(i+1:final_num) = Crorna(i) + 1;
        elseif (next_v == 8)
            Crorna(i+1:final_num) = Crorna(i) - 1;
        else
            disp('something went wrong');
        end
    end    
    n=n+1;
    subplot(4,5,n);
    hold on
    plot(Ciprot,Croprot);
    plot(Ciprot(end),Croprot(end),'o')
    xlabel('[Ci port]'); ylabel('[Cro prot]')
    title('Ciprot VS Croprot')
end

%%%%%%%%%%%%%%%
%% 3
%%%%%%%%%%%%%%%

Xcirna=1.2; Xciprot=1.2;
Xcrorna=0.8; Xcroprot=0.8;
Wci=50;Uci=50;Wcro=50;Ucro=50;
Kci=10; Kcro=10;

initial_num = 0;
final_num = 50000;
num_reactions=8;

Ciprot = zeros(final_num,1);
Cirna = zeros(final_num,1);
Croprot = zeros(final_num,1);
Crorna = zeros(final_num,1);
time = zeros(final_num,1);

Ciprot(1)=0; Cirna(1)=0;
Croprot(1)=0; Crorna(1)=20;
v = zeros(1,num_reactions);
n=0;
for k=1:20
    for i=1:(final_num-1)
        v(1)=Wci*Cirna(i);
        v(2)=Xciprot*Ciprot(i);
        v(3)=Uci*(1-(Croprot(i)^2)/(Kcro^2+Croprot(i)^2));
        v(4)=Xcirna*Cirna(i);
        v(5)=Wcro*Crorna(i);
        v(6)=Xcroprot*Croprot(i);
        v(7)=Ucro*(1-(Ciprot(i)^2)/(Kci^2+Ciprot(i)^2));
        v(8)=Xcrorna*Crorna(i);
        alpha0=sum(v);
        y = rand();
        y_mod = (alpha0)*y;
        tau = -log(y)/(alpha0);
        time(i+1) = time(i)+tau;    
        next_v = 0;
        for j=1:num_reactions
            if (y_mod < sum(v(1:j)))
                next_v = j;
                break;
            end
        end 
        if (next_v == 1)
            Ciprot(i+1:final_num) = Ciprot(i) + 1;
        elseif (next_v == 2)
            Ciprot(i+1:final_num) = Ciprot(i) - 1;
        elseif (next_v == 3)
            Cirna(i+1:final_num) = Cirna(i) + 1;
        elseif (next_v == 4)
            Cirna(i+1:final_num) = Cirna(i) - 1;
        elseif (next_v == 5)
            Croprot(i+1:final_num) = Croprot(i) + 1;
        elseif (next_v == 6)
            Croprot(i+1:final_num) = Croprot(i) - 1;
        elseif (next_v == 7)
            Crorna(i+1:final_num) = Crorna(i) + 1;
        elseif (next_v == 8)
            Crorna(i+1:final_num) = Crorna(i) - 1;
        else
            disp('something went wrong');
        end
    end
    n=n+1;
    subplot(4,5,n);
    hold on
    plot(time,Ciprot,'red');plot(time,Cirna);
    plot(time,Croprot,'blue');plot(time,Crorna);
    hold off
    xlabel('Time (s)'); ylabel('Number of Molecules')
    title('Stochastic model')
    legend('Ciprot','Cirna','Croprot','Crorna')
    
end

%%%%%%%%%%%%%%%
%% 4.
%%%%%%%%%%%%%%%

Xcirna=1.2; Xciprot=10;
Xcrorna=0.8; Xcroprot=0.8;
Wci=50;Uci=50;Wcro=50;Ucro=50;
Kci=10; Kcro=10;

initial_num = 0;
final_num = 100000;
num_reactions=8;

Ciprot = zeros(final_num,1);
Cirna = zeros(final_num,1);
Croprot = zeros(final_num,1);
Crorna = zeros(final_num,1);
time = zeros(final_num,1);

Ciprot(1)=4000; Cirna(1)=0;
Croprot(1)=0; Crorna(1)=0;
v = zeros(1,num_reactions);
n=0;
for k=1:20
    for i=1:(final_num-1)
        v(1)=Wci*Cirna(i);
        v(2)=Xciprot*Ciprot(i);
        v(3)=Uci*(1-(Croprot(i)^2)/(Kcro^2+Croprot(i)^2));
        v(4)=Xcirna*Cirna(i);
        v(5)=Wcro*Crorna(i);
        v(6)=Xcroprot*Croprot(i);
        v(7)=Ucro*(1-(Ciprot(i)^2)/(Kci^2+Ciprot(i)^2));
        v(8)=Xcrorna*Crorna(i);
        alpha0=sum(v);
        y = rand();
        y_mod = (alpha0)*y;
        tau = -log(y)/(alpha0);
        time(i+1) = time(i)+tau;    
        next_v = 0;
        for j=1:num_reactions
            if (y_mod < sum(v(1:j)))
                next_v = j;
                break;
            end
        end 
        if (next_v == 1)
            Ciprot(i+1:final_num) = Ciprot(i) + 1;
        elseif (next_v == 2)
            Ciprot(i+1:final_num) = Ciprot(i) - 1;
        elseif (next_v == 3)
            Cirna(i+1:final_num) = Cirna(i) + 1;
        elseif (next_v == 4)
            Cirna(i+1:final_num) = Cirna(i) - 1;
        elseif (next_v == 5)
            Croprot(i+1:final_num) = Croprot(i) + 1;
        elseif (next_v == 6)
            Croprot(i+1:final_num) = Croprot(i) - 1;
        elseif (next_v == 7)
            Crorna(i+1:final_num) = Crorna(i) + 1;
        elseif (next_v == 8)
            Crorna(i+1:final_num) = Crorna(i) - 1;
        else
            disp('something went wrong');
        end
    end
    n=n+1;
    subplot(4,5,n);
    hold on
    plot(time,Ciprot,'red');plot(time,Cirna);
    plot(time,Croprot,'blue');plot(time,Crorna);
    hold off
    xlabel('Time (s)'); ylabel('Number of Molecules')
    title('Stochastic model')
    legend('Ciprot','Cirna','Croprot','Crorna')
    
end

