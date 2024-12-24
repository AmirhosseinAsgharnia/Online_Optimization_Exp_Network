clear; clc
%%

load("Four_node.mat")

b(:,1) = 1;
b(:,2) = 2;
b(:,3) = 0;
b(:,4) = 1;
%% Time

time_horizon = 20000; % (T)
T_M = 100; % Look back time. T_M is used to measure the importance of a reservation based on previous T_M repeat. 
K = 50; % Time window

%%
Expected_reserve = zeros(time_horizon , 1);
Expected_block   = zeros(time_horizon , 1);
%% Core Hyper Parameters

lambda = 141; % Lagrange multiplier (\lambda)
eta    = sqrt(1/T_M); % (\eta)
v      = 20; % Blocking cost threshold (v)
epsilon= 0.0; % Exploration rate (\varepsilon)
Beta_LR= .1; % Learning rate (\beta)

%% Estimator's Hyper Parameters (Router)

router.num                  = 4; % Number of gateway routers
router.jobs                 = 5;
router.num_MF               = 3;
router.num_rules            = router.num_MF ^ router.num;
router.input_bounds         = ones (router.num , 2);
router.input_bounds(: , 2)  = router.jobs;

router.rule_base            = zeros (router.num_rules , 1);

%% Estimator's Hyper Parameters (Server)

server.num                  = 11;
server.jobs                 = 4; % Possible reservations. Action 1: 0 reservations / Action 2: 1 resersavtion / ... / Action 5: 4 reservations
server.job_vector           = server.jobs * ones( 1 , server.num ); 
server.num_possible_actions = prod ( server.job_vector );
server.num_MF               = 2;
server.num_rules            = server.num_MF ^ server.num;
server.input_bounds         = ones (server.num , 2);
server.input_bounds(: , 2)  = server.jobs;

server.rule_base_R          = ones (server.num_rules , 1) / server.num_rules;
server.rule_base_C          = ones (server.num_rules , router.num_rules) / server.num_rules;
server.probability          = zeros (time_horizon , 1);
server.P                    = zeros (server.num_rules , time_horizon);
probability_aux_normalized_pre = ones(server.num_rules , 1) / server.num_rules;
%%

f_t = zeros( server.num_rules , K);
F_s = zeros( server.num_rules , T_M);
% 
% Repeat_b = zeros(num_rules , T_M);
Repeat_a = zeros(server.num_rules , T_M);

%% Do not touch this part.

n = 1;
m = server.jobs;

M = zeros(2*n + 1 , m);

for h = 1:m

    if h > n && h <= m - n
        counter_1 = 0;
        for j = -n:n
            counter_1 = counter_1 + 1;
            M(counter_1,h) = h + j;
        end
    end

end

for h = 1:n
    M(:,h)  = M(:,n+1);
    M(:,m-h+1) = M(:,m-n);
end

%%

tic;
rng(111)

for It = 1 : time_horizon

    if It ~= 1

        [C , R , evaluated_costs] = cost_function(b(It-1 , :) , Repeat_a , M , epsilon , n_k , server);

        counter_1 = 0;

        for evaluate_input = evaluated_costs'

            counter_1 = counter_1 + 1;
            [H{1:server.num}] = ind2sub( server.jobs * ones(1 , server.num) , evaluate_input);
            A_eval = cell2mat(H);

            fuzzy_ind_router  = fuzzy_engine_4 (b(It-1 , :) , zeros(router.num_rules , 1) , router.num_MF , router.input_bounds);

            fuzzy_R  = fuzzy_engine_11 (A_eval , server.rule_base_R , server.num_MF , server.input_bounds);
            Reward_cost  = R (counter_1);
            server.rule_base_R (fuzzy_R.act) = server.rule_base_R (fuzzy_R.act) + Beta_LR * (Reward_cost - fuzzy_R.res) * fuzzy_R.phi;

            counter_2 = 1;
            
            for h = fuzzy_ind_router.act'
                fuzzy_C = fuzzy_engine_11 (A_eval , server.rule_base_C(:,h) , server.num_MF , server.input_bounds);
                Reward_block = C (counter_1);
                server.rule_base_C(fuzzy_C.act , h) = server.rule_base_C (fuzzy_C.act,h) + Beta_LR * (Reward_block - fuzzy_C.res) * fuzzy_C.phi * fuzzy_ind_router.phi(counter_2);
                counter_2 = counter_2 + 1;
            end

        end

        f_t (: , 1 : end - 1) = f_t (: , 2 : end);
        
        counter_3 = 1;
        for h = fuzzy_ind_router.act'
            f_t (: , end) = f_t (: , end) + server.rule_base_C(: , h) * fuzzy_ind_router.phi(counter_3);
            counter_3 = counter_3 + 1;
        end
        
        FFS = max(0 , min(K , It + 1) ^ -1 * (sum(f_t(: , max(1 , K - It + 2) : K) , 2)) - v);
        
        F_s (: , end) = FFS;
        
        probability_aux = exp(- eta * sum(server.rule_base_R + lambda * F_s (: , 1:It - 1) , 2));
        probability_aux_normalized = probability_aux ./ sum(probability_aux);
        
        %
        Expected_block   (It - 1) = sum(f_t(:,end) .* probability_aux_normalized);
        Expected_reserve (It - 1) = sum(server.rule_base_R .* probability_aux_normalized);
        
        server.probability(It) = sqrt( sum((probability_aux_normalized_pre - probability_aux_normalized) .^2) );

        probability_aux_normalized_pre = probability_aux_normalized;



        [~ , d] = max(probability_aux_normalized);
        [X{1:11}] = ind2sub( server.jobs*ones(1,11) , d);
        Y = cell2mat (X);
        d = sub2ind( server.jobs*ones(1,11) , Y(11) , Y(10) ,Y(9) ,Y(8) ,Y(7) ,Y(6) ,Y(5) ,Y(4) ,Y(3) ,Y(2) ,Y(1));
        Repeat_a( : , 1:end-1 ) = Repeat_a( : , 2:end ); Repeat_a( : , end ) = 0;
        Repeat_a( d , end )     = 1;


    end

    %%
    
    clc
    fprintf("The process has been %.3f %% completed! \n",It * 100 / time_horizon)

end
Time = toc;
save(sprintf('Result_lambda_%.0f.mat',lambda))