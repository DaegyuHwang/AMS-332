%%%% Part A

%%%%%%%%%%%%%%%
%% 1
%%%%%%%%%%%%%%%

GNa = 400;
GK= 200;
GL = 2;
ENa = 99;
EK = -85;
VL = -65;
C = 2;
% time
totalt = 200;
deltat = 0.01;
numtime = totalt/deltat+1;
timearray = zeros(numtime,1);
V = zeros(numtime,1);
m = zeros(numtime,1);
h = zeros(numtime,1);
n = zeros(numtime,1);
% Input current
Ie = zeros(numtime,1);
Im = zeros(numtime,1);
% Initial conditions
V(1) = VL;
for t=1 : numtime-1
   am_V = 0.1*(V(t)+40)/(1-exp(-0.1*(V(t)+40)));
   ah_V = 0.07*exp(-0.05*(V(t)+65));
   an_V = 0.01*(V(t)+55)/(1-exp(-0.1*(V(t)+55)));
   Bm_V = 4*exp(-0.0556*(V(t)+65));
   Bh_V = 1/(1+exp(-0.1*(V(t)+35)));
   Bn_V = 0.125*exp(-0.0125*(V(t)+65));
   m(1) = am_V/(am_V+Bm_V);
   h(1) = ah_V/(ah_V+Bh_V);
   n(1) = an_V/(an_V+Bn_V);
   Im(1) = GL*(V(1)-VL) + GNa*(m(t)^3)*h(t)*(V(t)-ENa) + GK*(n(t)^4)*(V(t)-EK);
   dm_dt = am_V*(1-m(t))-Bm_V*m(t);
   dh_dt = ah_V*(1-h(t))-Bh_V*h(t);
   dn_dt = an_V*(1-n(t))-Bn_V*n(t);
   m(t+1) = m(t)+ dm_dt*deltat;
   h(t+1) = h(t)+ dh_dt*deltat;
   n(t+1) = n(t)+ dn_dt*deltat;
   if (t-1)*deltat > 40
       Ie(t) = 200;
   end
   dV_dt = 1/C*(-GL*(V(t)-VL)-GNa*(m(t)^3)*h(t)*(V(t)-ENa)-GK*(n(t)^4)*(V(t)-EK)+Ie(t));
   V(t+1) = V(t) + dV_dt*deltat;
   Im(t+1) = GL*(V(t+1)-VL) + GNa*(m(t+1)^3)*h(t+1)*(V(t+1)-ENa) + GK*(n(t+1)^4)*(V(t+1)-EK);
   timearray(t+1) = timearray(t)+deltat;
end
figure(1)
sgtitle('Hodgkin Huxley model')
subplot(5,1,1)
plot(timearray,V)
xlabel('time (ms)'); ylabel('membrane potential (mV)')
title('membrane potential vs time')
subplot(5,1,2)
plot(timearray,Im)
xlabel('time (ms)'); ylabel('total membrane current (pA)')
title('total membrane current vs time')
subplot(5,1,3)
plot(timearray,GNa.*(m.^3).*h)
xlabel('time (ms)'); ylabel('sodium conductance (nS)')
title('sodium conductance vs time')
subplot(5,1,4)
plot(timearray,GK.*(n.^4))
xlabel('time (ms)'); ylabel('potassium conductance (nS)')
title('potassium conductance vs time')
subplot(5,1,5)
plot(timearray, Ie)
xlabel('time (ms)'); ylabel('current Ie (pA)')
title('current Ie vs time')


%%%%%%%%%%%%%%%
%% 2
%%%%%%%%%%%%%%%

GNa = 400;
GK= 200;
GL = 2;
ENa = 99;
EK = -85;
VL = -65;
C = 2;
% about spike
V_spk = -20;
n_spk = 0;
n_t = zeros(numtime,1);
% time
totalt = 200;
deltat = 0.01;
numtime = totalt/deltat+1;
timearray = zeros(numtime,1);
V = zeros(numtime,1);
m = zeros(numtime,1);
h = zeros(numtime,1);
n = zeros(numtime,1);
% Input current
Ie = zeros(numtime,1);
Im = zeros(numtime,1);
% Initial conditions
V(1) = VL;
for t=1 : numtime-1
  am_V = 0.1*(V(t)+40)/(1-exp(-0.1*(V(t)+40)));
  ah_V = 0.07*exp(-0.05*(V(t)+65));
  an_V = 0.01*(V(t)+55)/(1-exp(-0.1*(V(t)+55)));
  Bm_V = 4*exp(-0.0556*(V(t)+65));
  Bh_V = 1/(1+exp(-0.1*(V(t)+35)));
  Bn_V = 0.125*exp(-0.0125*(V(t)+65));
  m(1) = am_V/(am_V+Bm_V);
  h(1) = ah_V/(ah_V+Bh_V);
  n(1) = an_V/(an_V+Bn_V);
  Im(1) = GL*(V(1)-VL) + GNa*(m(t)^3)*h(t)*(V(t)-ENa) + GK*(n(t)^4)*(V(t)-EK);
  dm_dt = am_V*(1-m(t))-Bm_V*m(t);
  dh_dt = ah_V*(1-h(t))-Bh_V*h(t);
  dn_dt = an_V*(1-n(t))-Bn_V*n(t);
  m(t+1) = m(t)+ dm_dt*deltat;
  h(t+1) = h(t)+ dh_dt*deltat;
  n(t+1) = n(t)+ dn_dt*deltat;
  if (t-1)*deltat > 40
      Ie(t) = 200;
  end
  dV_dt = 1/C*(-GL*(V(t)-VL)-GNa*(m(t)^3)*h(t)*(V(t)-ENa)-GK*(n(t)^4)*(V(t)-EK)+Ie(t));
  V(t+1) = V(t) + dV_dt*deltat;
  Im(t+1) = GL*(V(t+1)-VL) + GNa*(m(t+1)^3)*h(t+1)*(V(t+1)-ENa) + GK*(n(t+1)^4)*(V(t+1)-EK);
  timearray(t+1) = timearray(t)+deltat;
  if V(t+1) >= V_spk
      if V(t) < V_spk
          n_t(n_spk+1) = t+1;
          n_spk = n_spk+1;
      end      
  end
end
% the number of total spikes in a given time
n_spk
hold on
plot(timearray,V)
plot(timearray(n_t(1)),V(n_t(1)),'-o')
hold off
xlim([timearray(n_t(1)-500),timearray(n_t(1)+500)])
xlabel('time (ms)'); ylabel('membrane potential (mV)')
title('membrane potential vs time')

%%%%%%%%%%%%%%%
%% 3
%%%%%%%%%%%%%%%

GNa = 400;
GK= 200;
GL = 2;
ENa = 99;
EK = -85;
VL = -65;
C = 2;
% about spike
V_spk = -20;
n_spk = 0;
n_t = zeros(numtime,1);
% time
totalt = 200;
deltat = 0.01;
numtime = totalt/deltat+1;
timearray = zeros(numtime,1);
V = zeros(numtime,1);
m = zeros(numtime,1);
h = zeros(numtime,1);
n = zeros(numtime,1);
% Input current
Ie = zeros(numtime,1);
Im = zeros(numtime,1);
I_0 = 2;
% Initial conditions
V(1) = VL;
% at least value of input current that arises a spike
I_1 = 0;
for i=I_0 : 200
   for t=1 : numtime-1
       am_V = 0.1*(V(t)+40)/(1-exp(-0.1*(V(t)+40)));
       ah_V = 0.07*exp(-0.05*(V(t)+65));
       an_V = 0.01*(V(t)+55)/(1-exp(-0.1*(V(t)+55)));
       Bm_V = 4*exp(-0.0556*(V(t)+65));
       Bh_V = 1/(1+exp(-0.1*(V(t)+35)));
       Bn_V = 0.125*exp(-0.0125*(V(t)+65));
  
       m(1) = am_V/(am_V+Bm_V);
       h(1) = ah_V/(ah_V+Bh_V);
       n(1) = an_V/(an_V+Bn_V);
       Im(1) = GL*(V(1)-VL) + GNa*(m(t)^3)*h(t)*(V(t)-ENa) + GK*(n(t)^4)*(V(t)-EK);
  
       dm_dt = am_V*(1-m(t))-Bm_V*m(t);
       dh_dt = ah_V*(1-h(t))-Bh_V*h(t);
       dn_dt = an_V*(1-n(t))-Bn_V*n(t);
  
       m(t+1) = m(t)+ dm_dt*deltat;
       h(t+1) = h(t)+ dh_dt*deltat;
       n(t+1) = n(t)+ dn_dt*deltat;
  
       if (t-1)*deltat > 40
           Ie(t) = i;
       end
       dV_dt = 1/C*(-GL*(V(t)-VL)-GNa*(m(t)^3)*h(t)*(V(t)-ENa)-GK*(n(t)^4)*(V(t)-EK)+Ie(t));
       V(t+1) = V(t) + dV_dt*deltat;
  
       Im(t+1) = GL*(V(t+1)-VL) + GNa*(m(t+1)^3)*h(t+1)*(V(t+1)-ENa) + GK*(n(t+1)^4)*(V(t+1)-EK);
  
       timearray(t+1) = timearray(t)+deltat;
  
       if V(t+1) >= V_spk
           if V(t) < V_spk
               n_t(n_spk+1) = t+1;
               n_spk = n_spk+1;
               I_1 = i;               
           end       
       end       
   end
   if I_1 ~= 0
       break
   end
end
% at least value of input current that arises a spike
I_1

figure(1)
sgtitle('Hodgkin Huxley model')
subplot(2,1,1)
plot(timearray,V)
xlabel('time (ms)'); ylabel('membrane potential (mV)')
title('membrane potential vs time')
subplot(2,1,2)
plot(timearray, Ie)
xlabel('time (ms)'); ylabel('current Ie (pA)')
title('current Ie vs time')

%%%%%%%%%%%%%%%
%% 4
%%%%%%%%%%%%%%%

GNa = 400;
GK= 200;
GL = 2;
ENa = 99;
EK = -85;
VL = -65;
C = 2;
% about spike  
V_spk = -20;
n_t = zeros(numtime,1);
% time
totalt = 400;
deltat = 0.01;
numtime = totalt/deltat+1;
timearray = zeros(numtime,1);
V = zeros(numtime,1);
m = zeros(numtime,1);
h = zeros(numtime,1);
n = zeros(numtime,1);
% Input current
Ie = zeros(numtime,1);
Im = zeros(numtime,1);
I_0 = 2;
% Initial conditions
V(1) = VL;
% Let I_rh : at least value of input current that arises Rheobase
% Let c_spk(i) : the number of spike according to 'i', input current.
% Let V_rh : the membrane potential with I_rh per time
I_rh = 0;
c_spk = zeros(300,1);
V_rh = zeros(numtime,1);
I__rh = zeros(numtime,1);
for i=200 : -1 : I_0
   n_spk = 0;
   for t=1 : numtime-1
       am_V = 0.1*(V(t)+40)/(1-exp(-0.1*(V(t)+40)));
       ah_V = 0.07*exp(-0.05*(V(t)+65));
       an_V = 0.01*(V(t)+55)/(1-exp(-0.1*(V(t)+55)));
       Bm_V = 4*exp(-0.0556*(V(t)+65));
       Bh_V = 1/(1+exp(-0.1*(V(t)+35)));
       Bn_V = 0.125*exp(-0.0125*(V(t)+65));
  
       m(1) = am_V/(am_V+Bm_V);
       h(1) = ah_V/(ah_V+Bh_V);
       n(1) = an_V/(an_V+Bn_V);
       Im(1) = GL*(V(1)-VL) + GNa*(m(t)^3)*h(t)*(V(t)-ENa) + GK*(n(t)^4)*(V(t)-EK);
  
       dm_dt = am_V*(1-m(t))-Bm_V*m(t);
       dh_dt = ah_V*(1-h(t))-Bh_V*h(t);
       dn_dt = an_V*(1-n(t))-Bn_V*n(t);
  
       m(t+1) = m(t)+ dm_dt*deltat;
       h(t+1) = h(t)+ dh_dt*deltat;
       n(t+1) = n(t)+ dn_dt*deltat;
  
       if (t-1)*deltat > 40
           Ie(t) = i;
       end
       dV_dt = 1/C*(-GL*(V(t)-VL)-GNa*(m(t)^3)*h(t)*(V(t)-ENa)-GK*(n(t)^4)*(V(t)-EK)+Ie(t));
       V(t+1) = V(t) + dV_dt*deltat;
  
       Im(t+1) = GL*(V(t+1)-VL) + GNa*(m(t+1)^3)*h(t+1)*(V(t+1)-ENa) + GK*(n(t+1)^4)*(V(t+1)-EK);
  
       timearray(t+1) = timearray(t)+deltat;
  
       if V(t+1) >= V_spk
           if V(t) < V_spk
               n_t(n_spk+1) = t+1;
               n_spk = n_spk+1;             
           end       
       end
   end
   if i == 109
       V_rh = V;
       I__rh = Ie;
   end
   c_spk(i) = n_spk;
   if c_spk(i)+3 < c_spk(i+1)
%If the number of spikes dramatically changes(I set  3), it implies that
% we passed the Rhehobase
       I_rh = i+1;      
       break
   end
end
figure(1)
sgtitle('Hodgkin Huxley model')
subplot(2,2,1)
plot(timearray,V)
xlabel('time (ms)'); ylabel('membrane potential (mV)')
title('I = 108 (pA)')
subplot(2,2,3)
plot(timearray, Ie)
xlabel('time (ms)'); ylabel('current Ie (pA)')
title('I = 108 (pA)')
subplot(2,2,2)
plot(timearray,V_rh)
xlabel('time (ms)'); ylabel('membrane potential (mV)')
title('I = 109 (pA)')
subplot(2,2,4)
plot(timearray, I__rh)
xlabel('time (ms)'); ylabel('current Ie (pA)')
title('I = 109 (pA)')
% I got I_rh = 109 when the total time is 200
% Finally, I confirmed the I_rh by increasing the total time 200=>400.
% However, the value of I_rh is the same.
% Thus I_rh = 109   /   c_spk(109)=27  when the totaltime:400  / c_spk(108)=5


%%%%%%%%%%%%%%%
%% 5
%%%%%%%%%%%%%%%

GNa = 400;
GK= 200;
GL = 2;
ENa = 99;
EK = -85;
VL = -65;
C = 2;
% about spike  
V_spk = -20;
% time
totalt = 2000;
deltat = 0.01;
numtime = totalt/deltat+1;
timearray = zeros(numtime,1);
V = zeros(numtime,1);
m = zeros(numtime,1);
h = zeros(numtime,1);
n = zeros(numtime,1);
% Input current
Ie = zeros(numtime,1);
Im = zeros(numtime,1);
I_0 = 2;
% Initial conditions
V(1) = VL;
% Let I_rh : at least value of input current that arises Rheobase
I_rh = 0;
c_spk = zeros(300,1);
% firing rate
fr_m1 = zeros(101,1);
fr_m2 = zeros(101,1);
I_fr = zeros(101,1);
for i=180 : -1 : 80
   I_fr(i-79) = i;
   n_spk = 0;
   n_t = zeros(numtime,1);
   for t=1 : numtime-1
       am_V = 0.1*(V(t)+40)/(1-exp(-0.1*(V(t)+40)));
       ah_V = 0.07*exp(-0.05*(V(t)+65));
       an_V = 0.01*(V(t)+55)/(1-exp(-0.1*(V(t)+55)));
       Bm_V = 4*exp(-0.0556*(V(t)+65));
       Bh_V = 1/(1+exp(-0.1*(V(t)+35)));
       Bn_V = 0.125*exp(-0.0125*(V(t)+65));
  
       m(1) = am_V/(am_V+Bm_V);
       h(1) = ah_V/(ah_V+Bh_V);
       n(1) = an_V/(an_V+Bn_V);
       Im(1) = GL*(V(1)-VL) + GNa*(m(t)^3)*h(t)*(V(t)-ENa) + GK*(n(t)^4)*(V(t)-EK);
  
       dm_dt = am_V*(1-m(t))-Bm_V*m(t);
       dh_dt = ah_V*(1-h(t))-Bh_V*h(t);
       dn_dt = an_V*(1-n(t))-Bn_V*n(t);
  
       m(t+1) = m(t)+ dm_dt*deltat;
       h(t+1) = h(t)+ dh_dt*deltat;
       n(t+1) = n(t)+ dn_dt*deltat;
  
       if (t-1)*deltat > 40
           Ie(t) = i;
       end
       dV_dt = 1/C*(-GL*(V(t)-VL)-GNa*(m(t)^3)*h(t)*(V(t)-ENa)-GK*(n(t)^4)*(V(t)-EK)+Ie(t));
       V(t+1) = V(t) + dV_dt*deltat;
  
       Im(t+1) = GL*(V(t+1)-VL) + GNa*(m(t+1)^3)*h(t+1)*(V(t+1)-ENa) + GK*(n(t+1)^4)*(V(t+1)-EK);
  
       timearray(t+1) = timearray(t)+deltat;
  
       if V(t+1) >= V_spk
           if V(t) < V_spk
               n_t(n_spk+1) = t+1;
               n_spk = n_spk+1;             
           end
       end
   end
   c_spk(i) = n_spk;
   if c_spk(i)+3 < c_spk(i+1) && I_rh==0
       I_rh = i+1;
   end
   if I_rh > i
       fr_m1(i-79) = 0;
       fr_m2(i-79) = 0;
   else
       fr_m1(i-79) = n_spk/(totalt-40);
       fr_m2(i-79) = 1/((n_t(2)-n_t(1))*deltat);
   end
end
figure(1)
sgtitle('Hodgkin Huxley model')
subplot(2,1,1)
plot(I_fr,fr_m1,'-o' )
xlabel('Input current (pA)'); ylabel('Firing frequency (Hz)')
title('f-I curve using firing method 1')
subplot(2,1,2)
plot(I_fr, fr_m2,'-o')
xlabel('Input current (pA)'); ylabel('Firing frequency (Hz)')
title('f-I curve using firing method 2')

% fr_m1 and fr_m2 contain firing rates for each input current
% which is from 80 to 130


%%%% Part B
%%%%%%%%%%%%%%%
%% 1
%%%%%%%%%%%%%%%



% C = 1 nF = 1000 pF
GL = 50;
VL = -65;
C = 1000;
% time
totalt = 200;
deltat = 0.1;
numtime = totalt/deltat+1;
timearray = zeros(numtime,1);
% about spike
V_spk = -45;
V_reset = -65;
t_arp = 2;
n_spk = 0;
n_t = zeros(numtime,1);
V = zeros(numtime,1);
% Input current     1.1nA = 1100pA
Ie = zeros(numtime,1);
% Initial conditions
V(1) = VL;
for t=1 : numtime-1
   if (t-1)*deltat > 40
       Ie(t) = 1100;
   end
   timearray(t+1) = timearray(t)+deltat;
   if n_spk ~= 0
       if n_t(n_spk)+t_arp/deltat > t
           V(t+1) = V_reset;
           continue;
       end
   end
   dV_dt = 1/C*(-GL*(V(t)-VL)+Ie(t));
   V(t+1) = V(t) + dV_dt*deltat;
   if V(t+1) >= V_spk
       if V(t) < V_spk
           n_t(n_spk+1) = t+1;
           n_spk = n_spk+1;
           V(t+1) = V_spk;
       end       	
   end
end
figure(1)
sgtitle('Leaky Integrate-and-Fire neuron model')
subplot(2,1,1)
plot(timearray,V,'-o')
xlabel('time (ms)'); ylabel('membrane potential (mV)')
title('membrane potential vs time')
subplot(2,1,2)
plot(timearray,Ie)
xlabel('time (ms)'); ylabel('external current (pA)')
title('external current vs time')

%%%%%%%%%%%%%%%
%% 2
%%%%%%%%%%%%%%%

% C = 1 nF = 1000 pF
GL = 50;
VL = -65;
C = 1000;
% time 
totalt = 10000;
deltat = 0.1;
numtime = totalt/deltat+1;
timearray = zeros(numtime,1);
% about spike
V_spk = -45;
V_reset = -65;
t_arp = 2;
% Let I_rh : at least value of input current that arises Rheobase
% Let c_spk(i) : the number of spike according to 'i', input current.
I_rh = 0;
c_spk = zeros(numtime,1);
V = zeros(numtime,1); % all membrane totential depends on time
% firing rate (using method 2)
fr_m2 = zeros(2100-950+1,1);
I_fr = zeros(2100-950+1,1);
% Input current     1.1nA = 1100pA
Ie = zeros(numtime,1);
% Initial conditions
V(1) = VL;
% i means Ie
for i=2100 : -1 : 950
   I_fr(i-949) = i;
   n_spk = 0;
   n_t = zeros(numtime,1); % contains particles mean the time of each spike arises
   for t=1 : numtime-1
       if (t-1)*deltat > 40
           Ie(t) = i;
       end
       timearray(t+1) = timearray(t)+deltat;
       if n_spk ~= 0
           if n_t(n_spk)+t_arp/deltat > t
               V(t+1) = V_reset;
               continue;
           end
       end
  
       dV_dt = 1/C*(-GL*(V(t)-VL)+Ie(t));
       V(t+1) = V(t) + dV_dt*deltat;
  
       if V(t+1) >= V_spk
           if V(t) < V_spk
               n_t(n_spk+1) = t+1;
               n_spk = n_spk+1;
               V(t+1) = V_spk;
           end       
       end
   end
   c_spk(i) = n_spk;
   if c_spk(i)+10 < c_spk(i+1) && I_rh ==0
       I_rh = i+1;      
   end
   if I_rh > i
       fr_m2(i-949) = 0;
   else
       fr_m2(i-949) = 1/((n_t(2)-n_t(1))*deltat);
   end
end
plot(I_fr, fr_m2,'-o')
xlabel('I (pA)'); ylabel('firing rate (spikes/ms)')
title('f-I curve using firing method 2')



%%%% 3
%%%%%%%%%%%%%%%
%% 3-1
%%%%%%%%%%%%%%%

GNa = 400;
GK= 200;
GL = 2;
ENa = 99;
EK = -85;
VL = -65;
C = 2;
% about spike  
V_spk = -20;
% time
totalt = 2000;
deltat = 0.05;
numtime = totalt/deltat+1;
% firing rate
  
fr_m2 = zeros(101,1);
I_fr = zeros(101,1);
  
I_fr_2 = zeros(101,1);
fr_m2_2 = zeros(101,1);
for deltat = 0.05 : 0.05 : 0.1
   timearray = zeros(numtime,1);
   V = zeros(numtime,1);
   m = zeros(numtime,1);
   h = zeros(numtime,1);
   n = zeros(numtime,1);
  
   % Input current
   Ie = zeros(numtime,1);
   Im = zeros(numtime,1);
   I_0 = 2;
  
   % Initial conditions
   V(1) = VL;
  
  
   % Let I_rh : at least value of input current that arises Rheobase
   I_rh = 0;
   c_spk = zeros(300,1);
  
   for i=180 : -1 : 80
       I_fr(i-79) = i;
       n_spk = 0;
       n_t = zeros(numtime,1);
       for t=1 : numtime-1
           am_V = 0.1*(V(t)+40)/(1-exp(-0.1*(V(t)+40)));
           ah_V = 0.07*exp(-0.05*(V(t)+65));
           an_V = 0.01*(V(t)+55)/(1-exp(-0.1*(V(t)+55)));
           Bm_V = 4*exp(-0.0556*(V(t)+65));
           Bh_V = 1/(1+exp(-0.1*(V(t)+35)));
           Bn_V = 0.125*exp(-0.0125*(V(t)+65));
      
           m(1) = am_V/(am_V+Bm_V);
           h(1) = ah_V/(ah_V+Bh_V);
           n(1) = an_V/(an_V+Bn_V);
           Im(1) = GL*(V(1)-VL) + GNa*(m(t)^3)*h(t)*(V(t)-ENa) + GK*(n(t)^4)*(V(t)-EK);
      
           dm_dt = am_V*(1-m(t))-Bm_V*m(t);
           dh_dt = ah_V*(1-h(t))-Bh_V*h(t);
           dn_dt = an_V*(1-n(t))-Bn_V*n(t);
      
           m(t+1) = m(t)+ dm_dt*deltat;
           h(t+1) = h(t)+ dh_dt*deltat;
           n(t+1) = n(t)+ dn_dt*deltat;
      
           if (t-1)*deltat > 40
               Ie(t) = i;
           end
           dV_dt = 1/C*(-GL*(V(t)-VL)-GNa*(m(t)^3)*h(t)*(V(t)-ENa)-GK*(n(t)^4)*(V(t)-EK)+Ie(t));
           V(t+1) = V(t) + dV_dt*deltat;
      
           Im(t+1) = GL*(V(t+1)-VL) + GNa*(m(t+1)^3)*h(t+1)*(V(t+1)-ENa) + GK*(n(t+1)^4)*(V(t+1)-EK);
      
           timearray(t+1) = timearray(t)+deltat;
      
           if V(t+1) >= V_spk
               if V(t) < V_spk
                   n_t(n_spk+1) = t+1;
                   n_spk = n_spk+1;             
               end
           end
       end
       c_spk(i) = n_spk;
       if c_spk(i)+3 < c_spk(i+1) && I_rh==0
           I_rh = i+1;
       end
  
       if I_rh > i
           fr_m1(i-79) = 0;
           fr_m2(i-79) = 0;
       else
           fr_m1(i-79) = n_spk/(totalt-40);
           fr_m2(i-79) = 1/((n_t(2)-n_t(1))*deltat);
       end
   end
   if deltat == 0.05
       I_fr_2 = I_fr;
       fr_m2_2 = fr_m2;
   end
 
end
figure(1)
sgtitle('Hodgkin Huxley model according to "dt" ')
subplot(2,1,1)
plot(I_fr_2,fr_m2_2,'-o' )
xlabel('I  (pA)'); ylabel('f  (spikes/ms)')
title('f-I curve with dt = 0.05ms')
subplot(2,1,2)
plot(I_fr, fr_m2,'-o')
xlabel('I  (pA)'); ylabel('f  (spikes/ms)')
title('f-I curve with dt = 0.1ms ')


%%%%%%%%%%%%%%%
%% 3-2
%%%%%%%%%%%%%%%

% C = 1 nF = 1000 pF
GL = 50;
VL = -65;
C = 1000;

% time 
totalt = 10000;
deltat = 0.2;
numtime = totalt/deltat+1;
timearray = zeros(numtime,1);

% about spike
V_spk = -45;
V_reset = -65;
t_arp = 2;

% Let I_rh : at least value of input current that arises Rheobase
% Let c_spk(i) : the number of spike according to 'i', input current.

I_rh = 0;
c_spk = zeros(numtime,1);
V = zeros(numtime,1); % all membrane totential depends on time

% firing rate (using method 2)
fr_m2 = zeros(2100-950+1,1);
I_fr = zeros(2100-950+1,1);


% Input current     1.1nA = 1100pA
Ie = zeros(numtime,1);

% Initial conditions
V(1) = VL;

% i means Ie
for i=2100 : -1 : 950
    I_fr(i-949) = i;
    n_spk = 0;
    n_t = zeros(numtime,1); % contains particles mean the time of each spike arises
    for t=1 : numtime-1
        if (t-1)*deltat > 40
            Ie(t) = i;
        end
        timearray(t+1) = timearray(t)+deltat;

        if n_spk ~= 0
            if n_t(n_spk)+t_arp/deltat > t
                V(t+1) = V_reset;
                continue;
            end
        end
    
        dV_dt = 1/C*(-GL*(V(t)-VL)+Ie(t));
        V(t+1) = V(t) + dV_dt*deltat;
    
        if V(t+1) >= V_spk
            if V(t) < V_spk
                n_t(n_spk+1) = t+1;
                n_spk = n_spk+1;
                V(t+1) = V_spk;
            end        
        end
    end
    c_spk(i) = n_spk;
    if c_spk(i)+5 < c_spk(i+1) && I_rh ==0
        I_rh = i+1;       
    end

    if I_rh > i 
        fr_m2(i-949) = 0;
    else 
        fr_m2(i-949) = 1/((n_t(2)-n_t(1))*deltat);
    end
end

plot(I_fr, fr_m2,'-o')
xlabel('I (pA)'); ylabel('f (spikes/ms)')
title('f-I curve  with  LIF ')



%%%%%%%%%%%%%%%
%% 3-3
%%%%%%%%%%%%%%%


tic

% C = 1 nF = 1000 pF
GL = 50;
VL = -65;
C = 1000;

% time 
totalt = 2000;
deltat = 0.002;
numtime = totalt/deltat+1;
timearray = zeros(numtime,1);

% about spike
V_spk = -45;
V_reset = -65;
t_arp = 2;

% Let I_rh : at least value of input current that arises Rheobase
% Let c_spk(i) : the number of spike according to 'i', input current.

I_rh = 0;
c_spk = zeros(numtime,1);
V = zeros(numtime,1); % all membrane totential depends on time

% firing rate (using method 2)
fr_m2 = zeros(2100-950+1,1);
I_fr = zeros(2100-950+1,1);


% Input current     1.1nA = 1100pA
Ie = zeros(numtime,1);

% Initial conditions
V(1) = VL;

% i means Ie
for i=2100 : -1 : 950
    I_fr(i-949) = i;
    n_spk = 0;
    n_t = zeros(numtime,1); % contains particles mean the time of each spike arises
    for t=1 : numtime-1
        if (t-1)*deltat > 40
            Ie(t) = i;
        end
        timearray(t+1) = timearray(t)+deltat;

        if n_spk ~= 0
            if n_t(n_spk)+t_arp/deltat > t
                V(t+1) = V_reset;
                continue;
            end
        end
    
        dV_dt = 1/C*(-GL*(V(t)-VL)+Ie(t));
        V(t+1) = V(t) + dV_dt*deltat;
    
        if V(t+1) >= V_spk
            if V(t) < V_spk
                n_t(n_spk+1) = t+1;
                n_spk = n_spk+1;
                V(t+1) = V_spk;
            end        
        end
    end
    c_spk(i) = n_spk;
    if c_spk(i)+5 < c_spk(i+1) && I_rh ==0
        I_rh = i+1;       
    end

    if I_rh > i 
        fr_m2(i-949) = 0;
    else 
        fr_m2(i-949) = 1/((n_t(2)-n_t(1))*deltat);
    end
end

plot(I_fr, fr_m2,'-o')
xlabel('I (pA)'); ylabel('f (spikes/ms)')
title('f-I curve  with  LIF ')

toc

%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%


tic


GNa = 400;
GK= 200;
GL = 2;
ENa = 99;
EK = -85;
VL = -65;
C = 2;

% about spike   
V_spk = -20;


% time
totalt = 2000;
deltat = 0.002;
numtime = totalt/deltat+1;

timearray = zeros(numtime,1);
V = zeros(numtime,1);
m = zeros(numtime,1);
h = zeros(numtime,1);
n = zeros(numtime,1);

% Input current
Ie = zeros(numtime,1);
Im = zeros(numtime,1);
I_0 = 2;

% Initial conditions
V(1) = VL;


% Let I_rh : at least value of input current that arises Rheobase
I_rh = 0;
c_spk = zeros(300,1);

% firing rate
fr_m1 = zeros(101,1);
I_fr = zeros(101,1);


for i=180 : -1 : 80
    I_fr(i-79) = i;
    n_spk = 0;
    n_t = zeros(numtime,1);
    for t=1 : numtime-1
        am_V = 0.1*(V(t)+40)/(1-exp(-0.1*(V(t)+40)));
        ah_V = 0.07*exp(-0.05*(V(t)+65));
        an_V = 0.01*(V(t)+55)/(1-exp(-0.1*(V(t)+55)));
        Bm_V = 4*exp(-0.0556*(V(t)+65));
        Bh_V = 1/(1+exp(-0.1*(V(t)+35)));
        Bn_V = 0.125*exp(-0.0125*(V(t)+65));
    
        m(1) = am_V/(am_V+Bm_V);
        h(1) = ah_V/(ah_V+Bh_V);
        n(1) = an_V/(an_V+Bn_V);
        Im(1) = GL*(V(1)-VL) + GNa*(m(t)^3)*h(t)*(V(t)-ENa) + GK*(n(t)^4)*(V(t)-EK);
    
        dm_dt = am_V*(1-m(t))-Bm_V*m(t);
        dh_dt = ah_V*(1-h(t))-Bh_V*h(t);
        dn_dt = an_V*(1-n(t))-Bn_V*n(t);
    
        m(t+1) = m(t)+ dm_dt*deltat;
        h(t+1) = h(t)+ dh_dt*deltat;
        n(t+1) = n(t)+ dn_dt*deltat;
    
        if (t-1)*deltat > 40
            Ie(t) = i;
        end
        dV_dt = 1/C*(-GL*(V(t)-VL)-GNa*(m(t)^3)*h(t)*(V(t)-ENa)-GK*(n(t)^4)*(V(t)-EK)+Ie(t));
        V(t+1) = V(t) + dV_dt*deltat;
    
        Im(t+1) = GL*(V(t+1)-VL) + GNa*(m(t+1)^3)*h(t+1)*(V(t+1)-ENa) + GK*(n(t+1)^4)*(V(t+1)-EK);
    
        timearray(t+1) = timearray(t)+deltat;
    
        if V(t+1) >= V_spk 
            if V(t) < V_spk
                n_t(n_spk+1) = t+1;
                n_spk = n_spk+1;              
            end
        end
    end
    c_spk(i) = n_spk;
    if c_spk(i)+3 < c_spk(i+1) && I_rh==0
        I_rh = i+1; 
    end

    if I_rh > i 
        fr_m2(i-79) = 0;
    else 
        fr_m2(i-79) = 1/((n_t(2)-n_t(1))*deltat);
    end
end


plot(I_fr, fr_m2,'-o')
xlabel('I'); ylabel('f')
title('f-I curve using firing method 2')

toc

