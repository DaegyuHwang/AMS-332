%%%% Part A

%%%%%%%%%%%%%%
%% 1
%%%%%%%%%%%%%%
% let λ = lambda;

N = 50;
lambda = 10;
T = 5;
n = 2*lambda*T;
S = zeros(N,2*lambda*T);

for i=1:N
    u = rand(n,1);
    ISI = -log(u)/lambda;
    t = cumsum(ISI);
    for j=1:n
        S(i,j) = t(j);
    end
end

S(S>T) = NaN;
plot(S,1:N,'.k');
xlabel('time (s)'); ylabel('Number of trials')
title('Neural spike trains')

%%%%%%%%%%%%%%
%% 2
%%%%%%%%%%%%%%

N = 50;
lambda = 10;
T = 5;
n = 2*lambda*T;
S = zeros(N,2*lambda*T);

for i=1:N
    u = rand(n,1);
    ISI = -log(u)/lambda;
    t = cumsum(ISI);
    for j=1:n
        S(i,j) = t(j);
    end
end

S(S>T) = NaN;
plot(S,1:N,'.k');
xlabel('time (s)'); ylabel('Number of trials #')
title('Raster plot')

%%%%%%%%%%%%%%
%% 3
%%%%%%%%%%%%%%

% let λ = lambda;
N = 50;
lambda = 10;
T = 5;
n = 2*lambda*T;
S = zeros(N,2*lambda*T);

% number of total spikes across trials
sum_f = 0;
% average firing rate
a_f = 0;

for i=1:N
    u = rand(n,1);
    ISI = -log(u)/lambda;
    t = cumsum(ISI);
    f = 0;
    for j=1:n
        S(i,j) = t(j);
        if S(i,j)>T
            if f == 0
                f = j-1;
                sum_f = sum_f+f;
            end
        end
    end
end

S(S>T) = NaN;
a_f = sum_f/(N*T);

%%%%%%%%%%%%%%
%% 4
%%%%%%%%%%%%%%

% let λ = lambda;
N = 50;
lambda = 10;
T = 5;
n = 2*lambda*T;
S = zeros(N,2*lambda*T);
% let delta_t = dt
dt = 0.2;
n_b = T/dt;
timearray = zeros(n_b,1);
timearray(1) = 0.1;

% P = firing rates
P = zeros(n_b,1);


for i=1:N
    u = rand(n,1);
    ISI = -log(u)/lambda;
    t = cumsum(ISI);
    f = 0;
    for j=1:n
        S(i,j) = t(j);
    end
end
S(S>T) = NaN;

for tb=0:(n_b-1)
    if tb ~= n_b-1
        timearray(tb+2) = timearray(tb+1)+0.2;
    end
    n_s = 0;
    for j=1:n
        for i=1:N
            if (tb*0.2<=S(i,j)) && (S(i,j)<=tb*0.2+dt)
                n_s = n_s +1;
            end
        end
    end
    P(1+tb) = n_s/(N*dt);
end

plot(timearray, P)
xlabel('time (s)'); ylabel('firing rate (spk/s)')
title('PSTH')

%%%%%%%%%%%%%%
%% 5
%%%%%%%%%%%%%%

% let λ = lambda;
N = 50;
lambda = 10;
T = 5;
n = 2*lambda*T;
S = zeros(N,2*lambda*T);
% CV
CV = zeros(N,1);
trial = zeros(N,1);
a_CV = 0;


for i=1:N
    trial(i) = i;
    u = rand(n,1);
    ISI = -log(u)/lambda;
    t = cumsum(ISI);
    f = 0;
    for j=1:n
        S(i,j) = t(j);
    end
    t(t>T) = NaN;
    nn = ~isnan(t);
    n_s = length(find(nn));
    t1 = t(1:n_s);
    CV(i) = std(diff(t1))/mean(diff(t1));
end

S(S>T) = NaN;
a_CV = 1/N*(sum(CV));

plot(trial, CV)
xlabel('Number of trials #'); ylabel('CV of ISIs')
title('CV of ISIs in each trial')

%%%%%%%%%%%%%%
%% 6
%%%%%%%%%%%%%%

% let λ = lambda;
N = 50;
lambda = 10;
T = 5;
n = 2*lambda*T;
S = zeros(N,2*lambda*T);

% number of total spikes across trials
sum_f = 0;
% firing factor of all trials
all_f = zeros(N,1);

for i=1:N
    u = rand(n,1);
    ISI = -log(u)/lambda;
    t = cumsum(ISI);
    f = 0;
    for j=1:n
        S(i,j) = t(j);
    end
    t(t>T) = NaN;
    nn = ~isnan(t);
    n_s = length(find(nn));
    all_f(i) = n_s;
end
S(S>T) = NaN;

FF = var(all_f)/mean(all_f);

%%%%%%%%%%%%%%
%% 7
%%%%%%%%%%%%%%

% let λ = lambda;
N = 50;
lambda = 10;


T = [10, 100, 200, 400, 800, 1600, 3200, 6400];
% CV
CV_t = zeros(length(T),1);

% firing factor of all trials
all_f = zeros(N,1);
FF_t = zeros(length(T),1);

for k=1:length(T)
    n = 2*lambda*T(k);
    S = zeros(N,2*lambda*T(k));
    for i=1:N
        u = rand(n,1);
        ISI = -log(u)/lambda;
        t = cumsum(ISI);
        f = 0;
        for j=1:n
            S(i,j) = t(j);
        end
        t(t>T(k)) = NaN;
        nn = ~isnan(t);
        n_s = length(find(nn));
        all_f(i) = n_s;        
    end
    t_last = t(1:n_s);
    CV_t(k) = std(diff(t_last))/mean(diff(t_last));
    FF_t(k) = var(all_f)/mean(all_f);
end

% I took the last spike train for CV

hold on
plot(T,CV_t);
plot(T,FF_t);
hold off
xlabel('Time (s)'); ylabel('value')
title('FF & CV according to T')
legend('CV','FF')


%%%% Part B

%%%%1

%%%%%%%%%%%%%%
%% AND
%%%%%%%%%%%%%%

M=4;
N=2+1;
% η=eta : learning rate
eta = 1;
T = 500;

x = [1 1 -1; 1 -1 -1; -1 1 -1; -1 -1 -1]';
y_t = [1 -1 -1 -1]; % target 'AND'

w = zeros(N,1);

% performance =1 if correnct otherwise 0
error = zeros(T,1); 

presentation = zeros(T,1);


for t=1:T
    presentation(t) = t;
    index = randi(M);
    y = sign(w'*x(:,index));
    error(t) = sign(y_t(index)*y); 
    w = w+eta*(y_t(index)-y)*x(:,index);
end

error(error==-1)=0;

plot(presentation, error,'-o')
xlabel('presentation number #'); ylabel('performance');
title('performance vs. presentation number')
ylim([-0.2 1.2])

%%%%%%%%%%%%%%
%% OR
%%%%%%%%%%%%%%

M=4;
N=2+1;
% η=eta : learning rate
eta = 1;
T = 500;

x = [1 1 -1; 1 -1 -1; -1 1 -1; -1 -1 -1]';
y_t = [1 1 1 -1]; % target 'OR'

w = zeros(N,1);

% performance =1 or 0
error = zeros(T,1); 

presentation = zeros(T,1);


for t=1:T
    presentation(t) = t;
    index = randi(M);
    y = sign(w'*x(:,index));
    error(t) = sign(y_t(index)*y); 
    w = w+eta*(y_t(index)-y)*x(:,index);
end

error(error==-1)=0; % change -1 to 0

plot(presentation, error,'-o')
xlabel('presentation number #'); ylabel('performance');
title('[OR]   performance vs. presentation number')
ylim([-0.2 1.2])


%%%%%%%%%%%%%%
%% XOR
%%%%%%%%%%%%%%

M=4;
N=2+1;
% η=eta : learning rate
eta = 1;
T = 1000;

x = [1 1 -1; 1 -1 -1; -1 1 -1; -1 -1 -1]';
y_t = [-1 1 1 -1]; % target 'XOR'

w = zeros(N,1);

% performance =1 or 0
error = zeros(T,1); 

presentation = zeros(T,1);


for t=1:T
    presentation(t) = t;
    index = randi(M);
    y = sign(w'*x(:,index));
    error(t) = sign(y_t(index)*y); 
    w = w+eta*(y_t(index)-y)*x(:,index);
end

error(error==-1)=0; % change -1 to 0

plot(presentation, error,'-o')
xlabel('presentation number #'); ylabel('performance');
title('[XOR]   performance vs. presentation number')
ylim([-0.2 1.2])

%%%%% 2
%%%%%%%%%%%%%%
%% 2-1
%%%%%%%%%%%%%%

N = 50+1;
M = 40;
T = [100,1000,5000];%stepas

% η=eta : learning rate
eta = 1;

X = randi(2,N,M); % 2 or 1 random value
X(X==2) = -1;
X(end,:) = -1;
y_t = randi(2,M,1); % M x 1 vector
y_t(y_t==2) = -1;

per_t1 = zeros(100,1);
per_t2 = zeros(1000,1);
per_t3 = zeros(5000,1);


for i=T
    performace = zeros(i,1); 
    w=zeros(N,1);
    for t=1:i
        index = randi(M);
        y=sign(w'*X(:, index)); 
        performace(t)=sign(y_t(index)*y);
        w = w + eta*(y_t(index)-y)*X(:,index);
    end
    performace(performace==-1)=0;
    if i==T(1)
        per_t1 = performace;
    elseif i==T(2)
        per_t2 = performace;
    elseif i==T(3)
        per_t3 = performace;
    end
end


figure(1)
sgtitle('performance vs. presentation number')
subplot(3,1,1)
plot(per_t1, '.')
xlabel('presentation number #'); ylabel('performance');
title('<T=100>')
subplot(3,1,2)
plot(per_t2, '.')
xlabel('presentation number #'); ylabel('performance');
title('<T=1000>')
subplot(3,1,3)
plot(per_t3, '.')
xlabel('presentation number #'); ylabel('performance');
title('<T=5000>')

%%%%%%%%%%%%%%
%% 2-3
%%%%%%%%%%%%%%

N = 50+1;
M = 40;
T = 1000;%stepas
% η=eta : learning rate
eta = [0.1, 1, 10];
X = randi(2,N,M); % 2 or 1 random value
X(X==2) = -1;
X(end,:) = -1;
y_t = randi(2,M,1); % M x 1 vector
y_t(y_t==2) = -1;
per_t1 = zeros(T,1);
per_t2 = zeros(T,1);
per_t3 = zeros(T,1);
for i=eta
   performace = zeros(T,1);
   w=zeros(N,1);
   for t=1:T
       index = randi(M);
       y=sign(w'*X(:, index));
       performace(t)=sign(y_t(index)*y);
       w = w + i*(y_t(index)-y)*X(:,index);
   end
   performace(performace==-1)=0;
   if i==eta(1)
       per_t1 = performace;
   elseif i==eta(2)
       per_t2 = performace;
   elseif i==eta(3)
       per_t3 = performace;
   end
end
figure(1)
sgtitle('performance vs. presentation number')
subplot(3,1,1)
plot(per_t1, '.')
xlabel('presentation number #'); ylabel('performance');
title('<η=0.1>')
subplot(3,1,2)
plot(per_t2, '.')
xlabel('presentation number #'); ylabel('performance');
title('<η=1>')
subplot(3,1,3)
plot(per_t3, '.')
xlabel('presentation number #'); ylabel('performance');
title('<η=10>')

%%%%%%%%%%%%%%
%% 2-4
%%%%%%%%%%%%%%

%  using 15 runs
all_count = zeros(15,1);



figure(1)
sgtitle('performance vs. presentation number')

for k=1:15
    count = 0;
    N = 50+1;
    M = 40;
    T = 1000;%stepas
    
    % η=eta : learning rate
    eta = 1;
    
    X = randi(2,N,M); % 2 or 1 random value
    X(X==2) = -1;
    X(end,:) = -1;
    y_t = randi(2,M,1); % M x 1 vector
    y_t(y_t==2) = -1;
    
    w=zeros(N,1);
    performance = zeros(T,1); 
    presentation = zeros(T,1);
    
    for t=1:T
        presentation(t) = t;
        index = randi(M);
        y=sign(w'*X(:, index)); 
        performance(t)=sign(y_t(index)*y);
        w = w + eta*(y_t(index)-y)*X(:,index);

        if performance(t)==1
            count = count +1;
            if count <= 200
                all_count(k)=t-199;                         
            end
        elseif performance(t)==-1
            count = 0;
        end
    end    
    performance(performance==-1)=0;

    subplot(5,3,k)
    plot(performance, '.')
    xlabel('presentation number #'); ylabel('performance');
end

average = sum(all_count)/15

%%%%%%%%%%%%%%
%% 2-5
%%%%%%%%%%%%%%


%  using 15 runs
all_count = zeros(15,1);

figure(1)
sgtitle('performance vs. presentation number')

for k=1:15
    count = 0;
    N = 50+1;
    M = 90;
    T = M*1000;%stepas
    
    % η=eta : learning rate
    eta = 1;
    
    X = randi(2,N,M); % 2 or 1 random value
    X(X==2) = -1;
    X(end,:) = -1;
    y_t = randi(2,M,1); % M x 1 vector
    y_t(y_t==2) = -1;
    
    w=zeros(N,1);
    performance = zeros(T,1); 
    presentation = zeros(T,1);
    
    for t=1:T
        presentation(t) = t;
        index = randi(M);
        y=sign(w'*X(:, index)); 
        performance(t)=sign(y_t(index)*y);
        w = w + eta*(y_t(index)-y)*X(:,index);

        if performance(t)==1
            count = count +1;
            if count <= 200
                all_count(k)=t-199;                         
            end
        elseif performance(t)==-1
            count = 0;
        end
    end    
    performance(performance==-1)=0;

    subplot(5,3,k)
    plot(performance, '.')
    xlabel('presentation number #'); ylabel('performance');
end

average = sum(all_count)/15
mean_pattern = average/M




























