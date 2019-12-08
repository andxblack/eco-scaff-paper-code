% Runs the model for the figure showing the GS dynamics.

reps = 1;
gens = 200;

% slow dispersal, with B's providing advantage to dispersal.
patch_fit = @(a,b) a*(1+200*b);
T = 30;

for ii=1:reps
    [av_beta,av_q,av_repo,av_sizeA,av_sizeB,av_patch_fit] = evo_dynamics2(gens,T,patch_fit);
end


% plot the trajectory
figure(1)
clf

subplot(2,3,1)
plot(av_beta)

xlim([0,200])
ylabel('mean(\beta)')
xlabel('generations')

subplot(2,3,2)
plot(av_q)
ylabel('mean(q)')
xlim([0,200])
xlabel('generations')

subplot(2,3,3)
plot(av_repo)
ylabel('mean(\beta(1-q))')
xlim([0,200])
xlabel('generations')

subplot(2,3,4)
plot(av_sizeA);
ylabel('mean(A)')
xlim([0,200])
xlabel('generations')


subplot(2,3,5)
plot(av_sizeB);
ylabel('mean(B)')
xlim([0,200])
xlabel('generations')

subplot(2,3,6)
plot(av_patch_fit)
ylabel('mean(patch fitness)')
xlim([0,200])
xlabel('generations')
