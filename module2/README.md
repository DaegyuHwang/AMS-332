### #Deterministic & stochastic model of the cro-cI genetic network.

When a bacteria is infected by the bacteriophage lambda (λ), two fates are possible. First, the virus may replicate many progeny within the bacteria, ultimately resulting in the bursting of the bacteria to release new phage; this is known as the lytic pathway, or lysis. However, in some cases, the phage DNA remains in the bacteria without making new phage; the phage DNA is replicated along with the bacterial genome during cell division, and thus all descendents of the infected cell also carry the phage DNA. This state is known as lysogeny; the lysogenic state is generally stable, meaning that all progeny of a lysogenic bacteria remain in the lysogenic state indefnitely. However, exposure to stress (such as radiation) can convert a cell in the lysogenic state into the lytic phase.

I constructed a model for cell-fate decisions regarding lysis and lysogeny in bacteriophage, considering their mutually inhibitory mechanism: a stochastic model of the Cro/CI genetic network. 
Unlike deterministic models, I was able to overcome the limitations of Mass-Action Kinetic models by introducing Gillespie's algorithm into this model and accounting for the randomness of molecular collisions that fundamentally drive chemical reactions. 


### #Used Computational model

Determinstic model of the cro cI genetic network:
![image](https://github.com/user-attachments/assets/5e7e058f-3f35-48e9-881c-1f422e3ea351)

Stochastic model of the cro cI genetic network:
![image](https://github.com/user-attachments/assets/00462655-f77f-47ac-879f-a1b6d6100565)


