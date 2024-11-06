%%%% Part A
%%%% 1 (a)
Vmax = 5;
K = 20;
S1 = []; S2=[]; S3=[];
V1 = []; V2=[]; V3=[];
% h = 1
h = 1;
for (i=0:100)
    S1(i+1) = i; 
    V1(i+1) = (Vmax*(S1(i+1)^h))/(K^h + S1(i+1)^h);
end
plot(S1,V1)
% h = 2
h = 2;
for (i=0:100)
    S2(i+1) = i; 
    V2(i+1) = (Vmax*(S2(i+1)^h))/(K^h + S2(i+1)^h);
end
hold on
plot(S2,V2)
% h = 10
h = 10;
for (i=0:100)
    S3(i+1) = i; 
    V3(i+1) = (Vmax*(S3(i+1)^h))/(K^h + S3(i+1)^h);
end
plot(S3,V3)
hold off
xlabel('Substrate concentration(mM)')
ylabel('Rate of change(mM/s)')
legend('h=1','h=2','h=10')
title('Hill equation 1(a)')

%%%%% 1 (b)
Vmax = 5;
h = 2;
S1 = []; S2=[]; S3=[];
V1 = []; V2=[]; V3=[];
%K =10
K = 10;
for (i=0:100)
    S1(i+1) = i; 
    V1(i+1) = (Vmax*(S1(i+1)^h))/(K^h + S1(i+1)^h);
end
plot(S1,V1)
%K =20
K = 20;
for (i=0:100)
    S2(i+1) = i; 
    V2(i+1) = (Vmax*(S2(i+1)^h))/(K^h + S2(i+1)^h);
end
hold on
plot(S2,V2)
%K =40
K = 40;
for (i=0:100)
    S3(i+1) = i; 
    V3(i+1) = (Vmax*(S3(i+1)^h))/(K^h + S3(i+1)^h);
end
plot(S3,V3)
hold off
xlabel('Substrate concentration(mM)')
ylabel('Rate of change(mM/s)')
legend('K=10','K=20','K=40')
title('Hill equation 1(b)')

%%%% 1(c)
K = 20;
h = 2;
S1 = []; S2=[]; S3=[];
V1 = []; V2=[]; V3=[];
%Vmax=2
Vmax = 2;
for (i=0:100)
    S1(i+1) = i; 
    V1(i+1) = (Vmax*(S1(i+1)^h))/(K^h + S1(i+1)^h);
end
plot(S1,V1)
%Vmax=5
Vmax =5;
for (i=0:100)
    S2(i+1) = i; 
    V2(i+1) = (Vmax*(S2(i+1)^h))/(K^h + S2(i+1)^h);
end
hold on
plot(S2,V2)
%Vmax=10
Vmax = 10;
for (i=0:100)
    S3(i+1) = i; 
    V3(i+1) = (Vmax*(S3(i+1)^h))/(K^h + S3(i+1)^h);
end
plot(S3,V3)
hold off
xlabel('Substrate concentration(mM)')
ylabel('Rate of change(mM/s)')
legend('Vmax=2','Vmax=5','Vmax=10')
title('Hill equation 1(c)')

%%%% 2(b)
K = 50; h = 4;
S1 = []; V1 = [];
Vmax = 20;
for (i=0:100)
    S1(i+1) = i; 
    V1(i+1) = (Vmax*(S1(i+1)^h))/(K^h + S1(i+1)^h);
end
plot(S1,V1)
xlabel('Substrate concentration(mM)')
ylabel('Rate of change(mM/s)')
title('Hill equation 2(b)')

%%%% 3 (a)
h=4; Vmax = 20; K=50;
totaltime = 10;
stepsize = 0.01;
numsteps = totaltime/stepsize;

% let S = [S] : the concetration of substrate
S = zeros(numsteps,1);
t = zeros(numsteps,1);

S(1)=100; t(1)=0;

for (i=1:numsteps-1)
    dS_dt = -Vmax*(S(i)^h)/(K^h + S(i)^h);
    S(i+1) = S(i) + dS_dt*stepsize;
    t(i+1) = t(i) + stepsize;
end

plot(t,S)
xlabel('time (s)')
ylabel('[S] (mM)')
title('Forward Euler algorithm 3(a)')



%%%% Part B
%%%% 1 (a)

% Xprot means [Xprot] Xrna means [Xrna] from the ODE
% CXprot means Xprot  CXrna means Xrna from the ODE
% K = K(1/2)
w=1; u=1; CXprot=1; CXrna=1; K=0.33;
totaltime = 20;
stepsize = 0.01;
numsteps = totaltime/stepsize;

Xprot = zeros(numsteps,1);
Xrna = zeros(numsteps,1);
t = zeros(numsteps,1);

Xprot(1)=0.5; Xrna(1)=0.5; t(1)=0;

for (i=1:numsteps-1)
    dXrna_dt = u*(Xprot(i)^2)/(K^2+Xprot(i)^2)-CXrna*Xrna(i);
    dXprot_dt = w*Xrna(i) - CXprot*Xprot(i);
    Xrna(i+1) = Xrna(i) + dXrna_dt*stepsize;
    Xprot(i+1) = Xprot(i) + dXprot_dt*stepsize;
    t(i+1) = t(i) + stepsize;
end
figure(1)
subplot(2,1,1)
plot(t,Xrna)
xlabel('time (s)')
ylabel('[Xrna] (mM)')
title('Forward Euler algorithm 1(a) about [Xrna]')
subplot(2,1,2)
plot(t,Xprot)
ylabel('[Xprot] (mM)')
xlabel('time (s)')
title('Forward Euler algorithm 1(a) about [Xprot]')

%%%% 1 (b)

w=1; u=1; CXprot=1; CXrna=1; K=0.33;
totaltime = 20;
stepsize = 0.01;
numsteps = totaltime/stepsize;

Xprot = zeros(numsteps,1);
Xrna = zeros(numsteps,1);
t = zeros(numsteps,1);

Xprot(1)=0.2; Xrna(1)=0; t(1)=0;

for (i=1:numsteps-1)
    dXrna_dt = u*(Xprot(i)^2)/(K^2+Xprot(i)^2)-CXrna*Xrna(i);
    dXprot_dt = w*Xrna(i) - CXprot*Xprot(i);
    Xrna(i+1) = Xrna(i) + dXrna_dt*stepsize;
    Xprot(i+1) = Xprot(i) + dXprot_dt*stepsize;
    t(i+1) = t(i) + stepsize;
end
figure(1)
subplot(2,1,1)
plot(t,Xrna)
xlabel('time (s)')
ylabel('[Xrna] (mM)')
title('Forward Euler algorithm 1(a) about [Xrna]')
subplot(2,1,2)
plot(t,Xprot)
ylabel('[Xprot] (mM)')
xlabel('time (s)')
title('Forward Euler algorithm 1(a) about [Xprot]')

%%%% 1 (c)
w=1; u=1; CXprot=1; CXrna=1; K=0.33;
totaltime = 20;
stepsize = 0.01;
numsteps = totaltime/stepsize;

Xprot = zeros(numsteps,1);
Xrna = zeros(numsteps,1);
t = zeros(numsteps,1);

Xprot(1)=0.5; Xrna(1)=0; t(1)=0;

for (i=1:numsteps-1)
    dXrna_dt = u*(Xprot(i)^2)/(K^2+Xprot(i)^2)-CXrna*Xrna(i);
    dXprot_dt = w*Xrna(i) - CXprot*Xprot(i);
    Xrna(i+1) = Xrna(i) + dXrna_dt*stepsize;
    Xprot(i+1) = Xprot(i) + dXprot_dt*stepsize;
    t(i+1) = t(i) + stepsize;
end
figure(1)
subplot(2,1,1)
plot(t,Xrna)
xlabel('time (s)')
ylabel('[Xrna] (mM)')
title('Forward Euler algorithm 1(a) about [Xrna]')
subplot(2,1,2)
plot(t,Xprot)
ylabel('[Xprot] (mM)')
xlabel('time (s)')
title('Forward Euler algorithm 1(a) about [Xprot]')

%%%% 1 (d)

w=1; u=1; CXprot=1; CXrna=1; K=0.33;
totaltime = 20;
stepsize = 0.01;
numsteps = totaltime/stepsize;

Xprot = zeros(numsteps,1);
Xrna = zeros(numsteps,1);
t = zeros(numsteps,1);

Xprot(1)=0; Xrna(1)=0.2; t(1)=0;

for (i=1:numsteps-1)
    dXrna_dt = u*(Xprot(i)^2)/(K^2+Xprot(i)^2)-CXrna*Xrna(i);
    dXprot_dt = w*Xrna(i) - CXprot*Xprot(i);
    Xrna(i+1) = Xrna(i) + dXrna_dt*stepsize;
    Xprot(i+1) = Xprot(i) + dXprot_dt*stepsize;
    t(i+1) = t(i) + stepsize;
end
figure(1)
subplot(2,1,1)
plot(t,Xrna)
xlabel('time (s)')
ylabel('[Xrna] (mM)')
title('Forward Euler algorithm 1(a) about [Xrna]')
subplot(2,1,2)
plot(t,Xprot)
ylabel('[Xprot] (mM)')
xlabel('time (s)')
title('Forward Euler algorithm 1(a) about [Xprot]')

%%%% 1 (e)
w=1; u=1; CXprot=1; CXrna=1; K=0.33;
totaltime = 20;
stepsize = 0.01;
numsteps = totaltime/stepsize;

Xprot = zeros(numsteps,1);
Xrna = zeros(numsteps,1);
t = zeros(numsteps,1);

Xprot(1)=0; Xrna(1)=0.5; t(1)=0;

for (i=1:numsteps-1)
    dXrna_dt = u*(Xprot(i)^2)/(K^2+Xprot(i)^2)-CXrna*Xrna(i);
    dXprot_dt = w*Xrna(i) - CXprot*Xprot(i);
    Xrna(i+1) = Xrna(i) + dXrna_dt*stepsize;
    Xprot(i+1) = Xprot(i) + dXprot_dt*stepsize;
    t(i+1) = t(i) + stepsize;
end
figure(1)
subplot(2,1,1)
plot(t,Xrna)
xlabel('time (s)')
ylabel('[Xrna] (mM)')
title('Forward Euler algorithm 1(a) about [Xrna]')
subplot(2,1,2)
plot(t,Xprot)
ylabel('[Xprot] (mM)')
xlabel('time (s)')
title('Forward Euler algorithm 1(a) about [Xprot]')

%%%% 2 (a)
w=1; u=1; CXprot=1; CXrna=1; K=0.33;
totaltime = 20;
stepsize = 0.01;
numsteps = totaltime/stepsize;

Xprot = zeros(numsteps,1);
Xrna = zeros(numsteps,1);
t = zeros(numsteps,1);
for (i=0:0.2:1.4)
    Xrna(1)=i;
    for (j=0:0.2:1.4)
        Xprot(1)=j;
        t(1)=0;
        for (i=1:numsteps-1)
            dXrna_dt = u*(Xprot(i)^2)/(K^2+Xprot(i)^2)-CXrna*Xrna(i);
            dXprot_dt = w*Xrna(i) - CXprot*Xprot(i);
            Xrna(i+1) = Xrna(i) + dXrna_dt*stepsize;
            Xprot(i+1) = Xprot(i) + dXprot_dt*stepsize;
            t(i+1) = t(i) + stepsize;
        end
    end
end

%%%% 2 (b)
w=1; u=1; CXprot=1; CXrna=1; K=0.33;
totaltime = 20;
stepsize = 0.01;
numsteps = totaltime/stepsize;

Xprot = zeros(numsteps,1);
Xrna = zeros(numsteps,1);
t = zeros(numsteps,1);
n=0;
for (i=0:0.2:1.4)
    Xrna(1)=i;
    for (j=0:0.2:1.4)
        Xprot(1)=j;
        t(1)=0; n=n+1;
        for (i=1:numsteps-1)
            dXrna_dt = u*(Xprot(i)^2)/(K^2+Xprot(i)^2)-CXrna*Xrna(i);
            dXprot_dt = w*Xrna(i) - CXprot*Xprot(i);
            Xrna(i+1) = Xrna(i) + dXrna_dt*stepsize;
            Xprot(i+1) = Xprot(i) + dXprot_dt*stepsize;
            t(i+1) = t(i) + stepsize;
        end
        hold on
        plot(Xrna,Xprot)
        plot(Xrna(end),Xprot(end),'o')
        xlabel('[Xrna] (mM)')
        ylabel('[Xprot] (mM)')
        title('[Xrna] vs [Xprot]')
    end
end
hold off

