function Cost = cost_function (reservations_index,job_request_index)

reservations = [1 1 1 1 1 1 1 1 1 1 1]';
job_request  = [2 2 2 2]';
W            = [20 20 20 20];
%%

resource_capacity = [zeros(4,1) ; zeros(4,1) ; reservations];

%%

[distance_vector , distance_index , G] = network_analysis ();

%%

[x1, x2, x3, x4] = ndgrid(0:job_request(1), 0:job_request(2), 0:job_request(3), 0:job_request(4));

job_request_aux = [x1(:), x2(:), x3(:), x4(:)];

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

Cost.T = Cost_T; 
Cost.V = Cost_V; 
Cost.R = Cost_R; 
