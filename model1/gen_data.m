% run the single type model with no extra variance and plot the results.

clear

% within-patch parameters.
N = 10^6;
mu = 0.02;  % size of mutation.
p = 0.05;  % prob of mutation.

patches = 64;
reps = 4;
maxG = 200;

% a function of the inital beta and the dispersal time.
patch_fun = @(beta,T) pois_sim(N,beta,mu,p,T);

b_slow = zeros(maxG,reps);
s_slow = zeros(maxG,reps);

tic
for ii=1:reps
    [b_slow(:,ii),s_slow(:,ii)] = evo_dynamics(0,0,patch_fun,patches,maxG,30);
end
toc


f5 = figure(5);
clf
subplot(1,2,1)
hold on
plot(b_slow,'color',[0.7,0.7,0.7]);
plot(mean(b_slow,2),'color',[0,0,0],'linewidth',2);
xlabel('generations')
ylabel('mean cell growth rate')


subplot(1,2,2)
hold on
plot(s_slow,'color',[0.7,0.7,0.7])
plot(mean(s_slow,2),'color',[0,0,0],'linewidth',2)

xlabel('generations')
ylabel('mean patch population')



