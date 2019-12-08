% Change the volume, not the concentration so that the peak time
% remains the same, but the size at that time is now variable. 
% These results are shown in the Supp.

function [av_beta,av_size,output] = single_run_Nscaling(patches,maxG,sigma_E,sigma_T)

    T = 30;
    mu = 0.02;  % size of mutation.
    p = 10^-2;  % prob of mutation.

    % initial values of beta for each patch.
    beta = 1.8*ones(patches,1);
    
    % output arrays.
    av_beta = zeros(maxG,1);
    av_size = zeros(maxG,1);
    
    % generate solutions for each patch
    for gen=1:maxG

        E0 = normrnd(10^6,sigma_E,patches,1);
        E0(E0<0) = 0;

        for ii=1:patches
            
            dT = T + sigma_T*randn(1);

            if E0(ii) == 0
                area{ii} = [0,0,0];
                beta_vec{ii} = [beta(ii),beta(ii),beta(ii)];
            else  
                [area{ii},beta_vec{ii}] = pois_sim_var(E0(ii),beta(ii),mu,p,dT,E0(ii));
            end
        end
        
        
        output(gen).area = area;
        output(gen).beta_vec = beta_vec;
        output(gen).beta = beta;

        % total size of each patch.
        for ii=1:patches
            pd(ii) = sum(area{ii});
        end

        av_beta(gen) = mean(beta);
        av_size(gen) = mean(pd);

        sam_patches = randsample(patches,patches,true,pd/sum(pd));

        for ii=1:patches
            foo = area{sam_patches(ii)};
            bar = randsample(length(foo),1,true,foo/sum(foo));
            beta(ii) = beta_vec{sam_patches(ii)}(bar);
        end


    end

end
