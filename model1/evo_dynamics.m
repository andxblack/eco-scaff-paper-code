
function [av_beta,av_size] = evo_dynamics(sigma1,sigma2,patch_fun,patches,maxG,T)
% run the evolutionary dynamics over some number of generations.
% sigma1 controls dispersal time for all patches within a generation
% sigma2 is for randomness within the geneartion.

    % initial values of beta for each patch.
    beta = 1.8*ones(patches,1);

    av_beta = zeros(maxG,1);
    av_size = zeros(maxG,1);

    % generate solutions for each patch
    for gen=1:maxG

        av_beta(gen) = mean(beta);
        
        randT = T + sigma1*rand(1);

        for ii=1:patches
            dT = randT + sigma2*randn(1);
            [area{ii},beta_vec{ii}] = patch_fun(beta(ii),dT);
        end

        % total area of each patch.
        for ii=1:patches
            pd(ii) = sum(area{ii});
        end

        av_size(gen) = mean(pd);

        sam_patches = randsample(patches,patches,true,pd/sum(pd));

        for ii=1:patches
            foo = area{sam_patches(ii)};
            bar = randsample(length(foo),1,true,foo/sum(foo));
            beta(ii) = beta_vec{sam_patches(ii)}(bar);
        end


    end

end
