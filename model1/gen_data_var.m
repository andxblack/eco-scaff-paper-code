% versions of the single type model that adds extra variablity in the
% dispersal time or other stuff.

clear

gens = 200;
patches = 100;

% variance in the growth times before dispersal.
% as for results shown in the maintext.
sigma_T = 2;
[av_beta1,av_size1,~] = single_run(patches,gens,0,sigma_T);

% variance in the concentration.
sigma_E0 = 2*10^4;
[av_beta2,av_size2,~] = single_run(patches,gens,sigma_E0,0);

% variance in `size' of the patch. Concentration remains fixed so peak time
% does not vary. 
sigma_E0 = 2*10^4;
[av_beta3,av_size3,~] = single_run_Nscaling(patches,gens,sigma_E0,0);
