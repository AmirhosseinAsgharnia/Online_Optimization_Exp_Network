function Cost = Cost_Calculator(A, B , F)

Cost.C_R  = 0;
Cost.C_V  = 0;
Cost.C_T  = 0;

D = A - B;
Res_2 = sum(D(D>0));
Res_3 = sum(D(D<0));
Cost_local_pre = 100000;
%%
Delta_available = zeros (numel(B) , numel(B));

for n = 1 : numel(B)
    for m = 1 : numel(B)
        
        if n == m
            
            Delta_available(n,m) = 0;
            
        else
            
            if B(n) - A(n) < 0; neg2zer_n = 0; else; neg2zer_n = B(n) - A(n); end
            if A(m) - B(m) < 0; neg2zer_m = 0; else; neg2zer_m = A(m) - B(m); end
            
            Delta_available(n,m) = min ( neg2zer_n , neg2zer_m );
            
        end
    end
end

K = [];
N = [];
M = [];

for n = 1 : numel(B)
    
    if sum(Delta_available(n,:)) == 0
        g = 0;
        N = [N,n];
        M = [M,n];
    else
        f = Delta_available(n,:);
        m = 1:numel(B);
        M = [M,m(f~=0)];
        g = f(f~=0);
        N = [N,n*ones(size(g))];
    end
    K = [K,g];
end

Delta_permutate = permutate(K);
Delta_permutate_res = Delta_permutate;
% Delta_permutate_res(sum(Delta_permutate,2)>Res_1,:)=[];
Delta_permutate_res(sum(Delta_permutate,2)>Res_2,:)=[];
Delta_permutate_res(sum(Delta_permutate,2)>abs(Res_3),:)=[];
max_cost_numbers = size(Delta_permutate_res,1);

Delta_permutate_sur = 0*Delta_permutate_res;

for cost_num = 1:max_cost_numbers

    f = find(Delta_permutate_res(cost_num,:) ~= 0);
    
    for j = 1 : numel(f)
        Delta_permutate_sur(cost_num , M(f(j))) = Delta_permutate_sur(cost_num , M(f(j))) - Delta_permutate_res(cost_num , f(j)); 
    end

end


for cost_num = 1:max_cost_numbers
    
    Cost_function_transfer = 0;
    
    for j = 1:numel(N)
        Cost_function_transfer = Cost_function_transfer+F.Ft{N(j),M(j)}(Delta_permutate_res(cost_num,j));
    end
    
    A_changed = A;
    
    Cost_function_violation = 0;
    
    for j = 1:numel(N)
        A_changed(N(j)) = A_changed(N(j)) + Delta_permutate_res(cost_num,j) + Delta_permutate_sur(cost_num,j);
    end
    
    D_2 = A_changed - B;
    D_3 = zeros(size(D_2));
    
    for j = 1:numel(D_2)
        if D_2(j) >= 0
            D_3(j) = 0;
        else
            D_3(j) = abs(D_2(j));
        end
    end
    
    for j = 1:numel(B)
        Cost_function_violation = Cost_function_violation+F.Fv{j}(D_3(j));
    end
    
    Cost_function_reservation = 0;
    
    D_4 = zeros(size(D_2));
    
    for j = 1:numel(D_2)
        if D_2(j) >= 0
            D_4(j) = D_2(j);
        else
            D_4(j) = 0;
        end
    end
    
    for j = 1:numel(B)
%         Cost_function_reservation = Cost_function_reservation + F.Fr(D_4(j));
        Cost_function_reservation = Cost_function_reservation + F.Fr(A(j));
    end
    
    Cost_local = min(Cost_function_violation + Cost_function_transfer, Cost_local_pre);
    
    if Cost_local ~= Cost_local_pre

        C_R = Cost_function_reservation;
        C_V = Cost_function_violation;
        C_T = Cost_function_transfer;

    end

    Cost_local_pre = Cost_local;
    
end

Cost.C_R = C_R;
Cost.C_V = C_V;
Cost.C_T = C_T;