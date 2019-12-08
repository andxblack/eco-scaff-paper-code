function [av_beta,av_q,av_repo,av_sizeA,av_sizeB,av_patch_fit] = evo_dynamics2(gens,T,pfFun)
    % evo dynamics for model with G and S types.
    
    % parameters
    eps = 10^-6;
    p = 10^-2;
    lambda = 1/eps * p;
    gamma = 1;
    max_mutations = 3;
    c = 1;
    d = 2;

    cols = 10;
    rows = 10;

    patches = rows*cols;

    % initial values of beta for each patch.
    beta = 1.8*ones(patches,1);
    q = 0*ones(patches,1);

    clear area

    b = zeros(patches,1);
    TA = zeros(patches,1);

    %gens = 200;

    av_patch_fit = zeros(gens,1);
    av_beta = zeros(gens,1);
    av_q = zeros(gens,1);
    av_repo = zeros(gens,1);
    av_sizeA = zeros(gens,1);
    av_sizeB = zeros(gens,1);
    av_part_fitness = zeros(gens,1);

    for gen=1:gens


    % create the circle objects that will be updated in the animation.
    for jj=1:cols 
    for kk=1:rows

        % counts how many we've plotted.
        g = (jj-1)*rows + kk;
        [area{g},b(g),beta_vec{g},q_vec{g}] = pois_sim_mod2(beta(g),q(g),gamma,T,lambda,eps,max_mutations,c,d);

    end
    end

    % calculate the average beta and q values
    av_beta(gen) = mean(beta);
    av_q(gen) = mean(q);
    av_repo(gen) = mean(beta.*(1-q));

    % calculate the 'fitnes' of each patch
    for ii=1:patches
        TA(ii) = pfFun(sum(area{ii}),b(ii));

    end
    
    for ii=1:patches
        patchAsize(ii) = sum(area{ii});
    end

    % average number of A in a given generation.
    av_sizeA(gen) = mean(patchAsize);
    av_sizeB(gen) = mean(b);
    
    % average fitness of the cells within a patch....
    av_part_fitness(gen) = sum(area{ii}.*beta_vec{ii}'.*(1-q_vec{ii}'))/(sum(area{ii})+b(ii));

    av_patch_fit(gen) = mean(TA);

    sam_patches = randsample(patches,patches,true,TA/sum(TA));

    % sample within the patch just according to the proportions of A's
    for ii=1:patches
        foo = area{sam_patches(ii)};

        bar = randsample(length(foo),1,true,foo/sum(foo));

        beta(ii) = beta_vec{sam_patches(ii)}(bar);
        q(ii) = q_vec{sam_patches(ii)}(bar);
    end

    end

end
