% Poission mutaion simulation, for model II

function [A_prop,B_prop,beta_vec,q_vec,t1_list] = pois_sim_mod2(beta,q,gamma,T,lambda,eps,maxS,c,d)

    m = 1; % number of spieces currently in the system.
    X0 = [1-eps,0,0,eps];
    beta_vec = beta;
    q_vec = q;
    t1_list = [];
    
   
    % mutation size parameters
    br = 0.05;
    qr = 0.02; 

    % initial growth of single type.
    y0 = X0(2);
    ie = 1; % event indicator.
    t = 0;
    
    while ~isempty(ie) && m < maxS
    
        % set event for next mutation time.
        u = rand;
        bum = @(t,y) Event2(t,y,y0,lambda,u);
        options = odeset('Events',bum);

        sol = ode45(@(t,y) dy_multi_pop2_cost(t,y,beta_vec,q_vec,gamma,c,d,m),[t T],X0,options);
        ie = sol.ie;
        Y1 = sol.ye;
      
        if numel(sol.xe) > 0
            t = sol.xe;
        else
            % if no event then this is empty, so must be at the end.
            t= T;
        end
        
        if ie == 1
            
            % only the initial type mutates.  
            mut_type = 1;
            
            [beta_new,q_new] = gen_mutant(beta_vec(mut_type),q_vec(mut_type),br,qr );
 
            m = m+1;
            beta_vec(m) = beta_new;
            q_vec(m) = q_new;
            Y1(1) = Y1(1) - eps;
            X0 = [Y1;eps];
            y0 = sol.ye(2);
            t1_list = [t1_list,t];
            

        end
            
            
    end
    
    if t < T
        % integrate forward without anymore mutations.
        bum = @(t,y) Event2(t,y,y0,0,0); % nothing should happen
        options = odeset('Events',bum);
        sol = ode45(@(t,y) dy_multi_pop2_cost(t,y,beta_vec,q_vec,gamma,c,d,m),[t T],X0,options);
    end
        
    A_prop = deval(sol,T,4:4+m-1);
    B_prop = deval(sol,T,3);

end

function [position,isterminal,direction] = Event2(t,y,y0,rate,u)
    
    position = y0 -log(u)/rate - y(2); % The value that we want to be zero
    isterminal = 1;  % Halt integration 
    direction = 0;   % The zero can be approached from either direction
    
end


function [bn,qn] = gen_mutant(beta_old,q_old,br,qr)

    % generate these independently on interval br*[-1,1]
    bn = beta_old + br*2*(rand-0.5);
    qn = q_old + qr*2*(rand-0.5);

    % limits on the process like this.
    if bn < 0
        bn=0;
    end
    
    if qn < 0
        qn = 0;
    end
    
    if bn > 2
        bn = 2;
    end
    

end


function dy = dy_multi_pop2_cost(t,y,beta,q,gamma,c,d,nMax)

    % Eq order: E, lambda, B, A_i
    
    dy = zeros(nMax+3,1);
        
    dy(1) = y(1)*((beta.*((1-c)*q-1))*y(4:end)) - y(3)*d*y(1);        % E eq
    dy(2) = y(1)*((beta.*(1-q))*y(4:end));                   % lambda eq    
    dy(3) = y(1)*(beta.*q)*y(4:end) -gamma*y(3);               % B eq 
    dy(4:end) = y(1)*((beta.*(1-q))'.*y(4:end))-gamma*y(4:end); 

end

