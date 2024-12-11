clear; clc

lambda = 141;

%%
MaxIt = 20000;
T_M = 100;

possible_actions = 5 * ones(1,11);
num_possible_actions = prod (possible_actions);
epsilon = 1;
num_evaluate_actions = ceil(num_possible_actions * epsilon);

num_inputs = 11;
num_membership_functions = 3;
num_rules = num_membership_functions ^ num_inputs;

input_bounds = ones (num_inputs , 2);
input_bounds(: , 2) = possible_actions';

rule_base_cost  = ones (num_rules , 1)/num_rules;
rule_base_block = ones (num_rules , num_rules)/num_rules;
rule_base       = ones (num_rules , 1)/num_rules;

f_t = zeros( num_rules , MaxIt);
F_s = zeros( num_rules , MaxIt);

Repeat_b = zeros(num_rules , T_M);
Repeat_a = zeros(num_rules , T_M);
%%

nu = 0.1;
Beta_LR = .01;
K = 141;

%%
P = zeros(num_rules , MaxIt);
P(:,1) = 1 / numel(P(:,1));
%%

n = 2;
m = 5;

M = zeros(2*n + 1 , m);

for h = 1:m

    if h > n && h <= m - n
        count = 0;
        for j = -n:n
            count = count + 1;
            M(count,h) = h + j;
        end
    end

end

for h = 1:n
    M(:,h)  = M(:,n+1);
    M(:,m-h+1) = M(:,m-n);
end
%%

tic;
eta = sqrt(1/MaxIt);
A_i = [1 1 1];

%%
Expected_reserve = zeros( num_rules , MaxIt);
Expected_block   = zeros( num_rules , MaxIt);
rng(111)
for It = 1 : MaxIt

    if It ~= 1

        [f , c , evaluated_costs] = cost_function(b(It-1 , :) , Cost , Repeat_a ,M);

        count = 0;

        for evaluate_input = evaluated_costs'

            count = count + 1;
            [H{1:numel([11 11 11])}] = ind2sub([11 11 11] , evaluate_input);
            A_eval = cell2mat(H);

            fuzzy_ind_b  = fuzzy_engine_3 (b(It-1 , :) , zeros(num_rules,1) , num_membership_functions , input_bounds);

            fuzzy_cost  = fuzzy_engine_3 (A_eval , rule_base_cost , num_membership_functions , input_bounds);
            Reward_cost  = c (count);
            rule_base_cost (fuzzy_cost.act) = rule_base_cost (fuzzy_cost.act) + Beta_LR * (Reward_cost - fuzzy_cost.res) * fuzzy_cost.phi;

            counter_2 = 1;
            for h = fuzzy_ind_b.act'
                fuzzy_block = fuzzy_engine_3 (A_eval , rule_base_block(:,h) , num_membership_functions , input_bounds);
                Reward_block = f (count);
                rule_base_block(fuzzy_block.act,h) = rule_base_block (fuzzy_block.act,h) + Beta_LR * (Reward_block - fuzzy_block.res) * fuzzy_block.phi * fuzzy_ind_b.phi(counter_2);
                counter_2 = counter_2 + 1;
            end

        end

        counter_2 = 1;
        for h = fuzzy_ind_b.act'
            f_t (: , It - 1) =  f_t (: , It - 1) + rule_base_block(:,h)* fuzzy_ind_b.phi(counter_2);
            counter_2 = counter_2 + 1;
        end

        FFS = max(0,min(K,It)^-1 * ( sum (f_t (:,It - min(It,K) + 1: It - 1) , 2)) - nu);


        F_s (: , It - 1) = FFS;
        rule_base = exp(-eta * sum( rule_base_cost + lambda * F_s (: , 1:It - 1) , 2));
        rule_base_norm = rule_base ./ sum(rule_base);

        %
        Expected_block   (:,It - 1) = f_t(: , It - 1)  .* rule_base_norm;
        Expected_reserve (:,It - 1) = rule_base_cost .* rule_base_norm;
        P(:,It) = rule_base_norm;

    end
    %%

    Y = sub2ind([11 11 11] , A(It,1)+1 , A(It,2)+1 , A(It,3)+1);
    Repeat_a( : , 1:end-1 ) = Repeat_a( : , 2:end ); Repeat_a( : , end ) = 0;
    Repeat_a( Y , end )     = 1;

    %%

    clc
    fprintf("The process has been %.3f %% completed! \n",It * 100 / MaxIt)

end
Time = toc;
save(sprintf('Result_lambda_%.0f.mat',lambda))