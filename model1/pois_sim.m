% Piecewise-deterministic simulation of within-patch dynamics.

function [XT,beta_final,sol_stuff] = pois_sim(N,beta_init,mu,p,T)

    beta1 = beta_init/N;
    gamma = 1;

    % 1. Solve upto the time of the first mutation.
    u = rand;
    bum = @(t,y) E1(t,y,p,u);
    options = odeset('Events',bum);

    X0 = [0; N-1;1];
    sol = ode45(@(t,y) eq1(t,y,beta1,gamma),[0 T],X0,options);

    t = sol.xe;
    X = sol.ye;
    
    t1_list = [0,t];
    sol_list = [sol];
  
    
    % create mutant and integrate forward until next mutation.
    
    bum = @(t,y) E1(t,y,p/2,rand);
    options = odeset('Events',bum);

    % stop the mutations going bigger than 2.
    % both mutants are smaller.
    if beta_init >= 2.0-mu
        beta2 = beta1 -mu/N;
        beta3 = beta1 -mu/N;
    else

        if rand < 0.5
            beta2 = beta1 + mu/N;
            beta3 = beta1 - mu/N;
        else 
            beta2 = beta1 - mu/N;
            beta3 = beta1 + mu/N;
        end
    end
    
    X0 = [0; X(2:3); 1];
    sol = ode45(@(t,y) eq2(t,y,beta1,beta2,gamma),[t T],X0,options);

    % create final mutant and grow until the end.
    t = sol.xe;
    X = sol.ye;
    
    t1_list = [t1_list, t];
    sol_list = [sol_list,sol];
    
    
    X0 = [0; X(2:4); 1];
    % have to still pass options here so that struct has the same fields.
    sol = ode45(@(t,y) eq3(t,y,beta1,beta2,beta3,gamma),[t T],X0,options);
    
    t1_list = [t1_list, T];
    sol_list = [sol_list,sol];
    

    XT = sol.y(3:5,end);
    beta_final = [beta1;beta2;beta3]*N;
    
    % useful output for creating plots.
    sol_stuff.t1_list = t1_list;
    sol_stuff.sol_list = sol_list;
    sol_stuff.beta_vec = beta_final;
    
    
end

function dy = eq1(t,y,beta,gamma)
    % equations before the first mutation.
    dy = zeros(3,1);
    
    dy(1) = y(2)*beta*y(3);
    dy(2) = -y(2)*beta*y(3);            % E eq
    dy(3) = y(2)*beta*y(3)-gamma*y(3);  % A1
  
end

function dy = eq2(t,y,beta1,beta2,gamma)
    % equations after the first mutation.
    dy = zeros(4,1);
    
    dy(1) = y(2)*beta1*y(3);
    dy(2) = -y(2)*beta1*y(3) - beta2*y(2)*y(4);        % E eq
    dy(3) = y(2)*beta1*y(3) - gamma*y(3);              % A1
    dy(4) = beta2*y(2)*y(4) - gamma*y(4);              % A2
  
end

function dy = eq3(t,y,beta1,beta2,beta3,gamma)
    % equations after the first mutation.
    dy = zeros(4,1);
    
    dy(1) = 0;
    dy(2) = -y(2)*beta1*y(3) - beta2*y(2)*y(4) - beta3*y(2)*y(5);                % E eq
    dy(3) = y(2)*beta1*y(3) - gamma*y(3);                % A1
    dy(4) = beta2*y(2)*y(4) - gamma*y(4);  % A2
    dy(5) = beta3*y(2)*y(5) - gamma*y(5);  % A3
    
  
end


function [position,isterminal,direction] = E1(t,y,rate,u)
  
    position = rate*y(1) +log(u); % The value that we want to be zero
    isterminal = 1;  % Halt integration 
    direction = 0;   % The zero can be approached from either direction
    
end
