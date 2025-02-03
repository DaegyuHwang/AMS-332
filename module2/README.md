### #Deterministic & stochastic model of the cro-cI genetic network

When a bacteria is infected by the bacteriophage lambda (λ), two fates are possible. First, the virus may replicate many progeny within the bacteria, ultimately resulting in the bursting of the bacteria to release new phage; this is known as the lytic pathway, or lysis. However, in some cases, the phage DNA remains in the bacteria without making new phage; the phage DNA is replicated along with the bacterial genome during cell division, and thus all descendents of the infected cell also carry the phage DNA. This state is known as lysogeny; the lysogenic state is generally stable, meaning that all progeny of a lysogenic bacteria remain in the lysogenic state indefnitely. However, exposure to stress (such as radiation) can convert a cell in the lysogenic state into the lytic phase.

The decision to enter lysis or lysogeny is based on a pair of mutually repressive transcription factors, cI and cro. When levels of cI are high and levels of cro are low, an infected bacteria will be in the lysogenic phase; when levels of cI are low and levels of cro are high, the lytic phase is prefered. Thus, we made our model based on the cI and cro.

I developed a model for cell-fate decisions regarding lysis and lysogeny in bacteriophage, considering their mutually inhibitory mechanism: a stochastic model of the Cro/CI genetic network. Unlike deterministic models, I was able to overcome the limitations of Mass-Action Kinetic models by introducing Gillespie's algorithm into this model and accounting for the randomness of molecular collisions that fundamentally drive chemical reactions. 


### #Used Computational model

Determinstic model of the cro cI genetic network:

![image](https://github.com/user-attachments/assets/be249935-5326-4d8b-b55e-a2fdf8d53d71)


Stochastic model of the cro cI genetic network:

![image](https://github.com/user-attachments/assets/35e619ea-c4a3-4538-8b19-6480ffdfb420)


