function [C , R , evaluated_costs] = cost_function ( b , Repeat_a , M , epsilon , n_k ,server)


%%

% b            = [2 2 2 2]; 
W            = [20 20 20 20];
job_request  = b';

%%

% Repeat_a = zeros(10000,1); Repeat_a (3250) = 1;
% epsilon = 0.0;

J = size(M,1);
r_a = sum (Repeat_a , 2);
[~ , sort_r_a] = sort (r_a,'descend');
a_ind = sort_r_a(1:1);

[H{1 : 11}] = ind2sub (server.jobs * ones(1 , 11) , a_ind);
A_n = cell2mat(H);

B_n = [];

for i = 1:11
    B_n = [B_n , M(:,A_n(i))];
end

C_n = [];
for i = 1:11
    C_n = [C_n , B_n(n_k( : , i) , i)];
end

a_max_ind = sub2ind(server.jobs * ones(1 , 11) , C_n(:,11) , C_n(:,10) , C_n(:,9) , C_n(:,8) , C_n(:,7) , C_n(:,6) , C_n(:,5) , C_n(:,4) , C_n(:,3) , C_n(:,2) , C_n(:,1) );

evaluated_costs = a_max_ind;

perm_vec = randperm(numel(evaluated_costs));
random_select = perm_vec(1 : ceil(epsilon*numel(evaluated_costs)));
evaluated_costs(random_select) = randi([1 , 362797056],numel(random_select),1);

perm_final = randperm(numel(evaluated_costs));
evaluated_costs = evaluated_costs(perm_final);

%%
[distance_vector , distance_index , ~] = network_analysis ();
[x1, x2, x3, x4] = ndgrid(0:job_request(1), 0:job_request(2), 0:job_request(3), 0:job_request(4));

job_request_aux = [x1(:), x2(:), x3(:), x4(:)];


evaluated_costs (randperm(ceil(numel(evaluated_costs)*0.999)))=[];


C = zeros(numel (evaluated_costs) , 1);
R = zeros(numel (evaluated_costs) , 1);



counter_1 = 0;
for instance = evaluated_costs'
    
    counter_1 = counter_1 + 1;
    [H{1 : 11}] = ind2sub (server.jobs * ones(1 , 11) , instance);
    reservations = cell2mat(H)';
    
    %%

    resource_capacity = [zeros(4,1) ; zeros(4,1) ; reservations];

    %%



    cost_pre = inf;

    for job_count = 1 : size(job_request_aux,1)

        cost_transfer    = 0;
        cost_violation   = 0;

        job_request_now = job_request_aux(job_count , :);
        resource_capacity_now = resource_capacity;

        for count = 1 : numel (distance_vector)

            if job_request_now(distance_index (count , 1)) ~= 0 && resource_capacity_now(distance_index (count , 2)) ~= 0

                D = min(job_request_now(distance_index (count , 1))  , resource_capacity_now(distance_index (count , 2)));

                job_request_now(distance_index (count , 1)) = job_request_now(distance_index (count , 1)) - D;
                resource_capacity_now(distance_index (count , 2)) = resource_capacity_now(distance_index (count , 2)) - D;

                cost_transfer = cost_transfer + D * distance_vector(count);

            end

            if sum(job_request_now) == 0 || sum(resource_capacity_now) == 0
                break;
            end

        end

        for count = 1:numel(job_request_now)

            job_request_diff = job_request(count) - (job_request_aux(job_count , count) - job_request_now(count) );

            cost_violation = cost_violation + W(count) * job_request_diff;

        end

        cost_reservation = sum(reservations);

        if cost_transfer + cost_violation < cost_pre
            cost_pre = cost_transfer + cost_violation;
            Cost_T = cost_transfer;
            Cost_V = cost_violation;
            Cost_R = cost_reservation;
        end
    end

    C(counter_1) = Cost_T + Cost_V;
    R(counter_1) = Cost_R;
end
