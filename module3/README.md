### #A diffusing autoregulatory gene & A diffusing pair of mutually inhibiting genes.

I made models for a diffusing pair of mutually inhibiting genes & A diffusing pair of mutually inhibiting genes considering Fick’s law of diffusion. To obtain meaningful results, mathematical models must be consistent with biological knowledge. Biological systems have spatial structures and most are mostly spatially heterogeneous. Since this affects the flux, flow rate, and concentration change rate for each spatial structure, I created a model that explains these aforementioned biological facts by introducing a Reaction-Diffusion model where concentrations vary in space and time using partial differential equations. 


### #Used Computational model
<A diffusing autoregulatory gene>
In Project 1, we implemented a ODE-based model for a simple autoregulatory gene; that is, a gene that encodes a protein whose function (or part of whose function) is to activate the transcription of itself. The reaction-diffusion model for this system is:

![image](https://github.com/user-attachments/assets/16f0ab63-6ec5-4a14-a0f7-43d01e52a242)

<A diffusing pair of mutually-inhibiting genes>
In Project 2, we implemented an ODE-based model of the cro/cI pair of mutually inhibitory genes. Here will consider a simplifed version of that system with symmetric parameters, in a reaction-diffusion model:

![image](https://github.com/user-attachments/assets/5612ea15-6589-44cd-816d-b5082142ba97)
