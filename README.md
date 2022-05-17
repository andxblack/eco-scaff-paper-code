# Ecological scaffolding paper code

Simulation code, written in MATLAB, to support the paper: 

Andrew J. Black, Pierrick Bourrat and Paul Rainey, 2019. 
Ecological scaffolding and the evolution of individuality during the transition from cells to multicellular life.
Nature Eco. Evo., 4, 426-436, [doi:10.1038/s41559-019-1086-9](https://dx.doi.org/10.1038/s41559-019-1086-9)


This is relatively raw research code with many parameters being hard coded. I have no plans to tidy this up, mainly as newer work has superseded these models.
If you use any of this code, please cite the above paper and note the requirements of the MIT license. 

- Tested on MATLAB R2019a.
- Requires the function 'randomsample.m' which is included in the Machine Learning and Statistics Toolboox.

An interactive visualisation of the single type model can be viewed at [Observable](https://observablehq.com/@andxblack/nested-darwinian-populations). The JavaScript for that can be viewed on that site. Due to the limitations of the ODEX solver, the JavaScript version uses an approximation for generating the mutation times that doesn't generalise well when extra variance is added to the dynamics.
