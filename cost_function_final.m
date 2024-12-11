function [f_t , C] = cost_function_final (b , Cost)
%%

% J = size(M,1);
% 
% epsilon = 0.2;
% 
% r_a = sum (Repeat_a , 2);
% 
% [~ , sort_r_a] = sort (r_a,'descend');
% a_ind = sort_r_a(1:1);
% [ H{1:numel([11 11 11])} ] = ind2sub( [11 11 11] , a_ind );
% A_n = cell2mat(H);
% B_n = [M(:,A_n(1)),M(:,A_n(2)),M(:,A_n(3))];
% 
% n_k = nchoosek([1:J 1:J 1:J],3);
% n_k = unique(n_k,'rows');
% 
% C_n = [B_n(n_k(:,1),1) B_n(n_k(:,2),2) B_n(n_k(:,3),3)];
% a_max_ind = sub2ind([11 11 11],C_n(:,3),C_n(:,2),C_n(:,1));
%%
possible_actions = [11 11 11 11 11 11];

Sub_1 = zeros(1331 , 3);
Sub_2 = zeros(1331 , 3);

% evaluated_costs = a_max_ind;
% 
% perm_vec = randperm(numel(evaluated_costs));
% random_select = perm_vec(1 : ceil(epsilon*numel(evaluated_costs)));
% evaluated_costs(random_select) = randi([1,1331],numel(random_select),1);

evaluated_costs = [1:1331]';
count = 0;

for It = evaluated_costs'
    
    count = count + 1;

    [H{1:numel([11 11 11])}] = ind2sub([11 11 11],It);
    A = cell2mat(H);

    Sub_1( count , :) = A;
    Sub_2( count , :) = b+1;

end

Sub_total = [Sub_1 , Sub_2];

input_index = sub2ind (possible_actions , Sub_total (:,1) , Sub_total (:,2) , Sub_total (:,3) , Sub_total (:,4) , Sub_total (:,5) , Sub_total (:,6));

f_t = sum(Cost(input_index,4:6) , 2) + sum(Cost(input_index,7:9) , 2);
C   = sum(Cost(input_index,1:3) , 2);