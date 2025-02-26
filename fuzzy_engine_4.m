function Output=fuzzy_engine_4(input,rule,MF_num,bound)
%% Mapping inside the boundary
N = numel(input);
%% Mapping

for k=1:N
    input(k)=  (input(k)-bound(k,1)) / (bound(k,2)-bound(k,1));
end

input = max (input,0);
input = min (input,1);
%%
MF=zeros(N,MF_num); % membership degree
R =zeros(2,N);

s = 1/(MF_num-1);
for i=1:N
    k=1;
    for j=1:MF_num
        if j==1
            if input(i)>=0 && input(i)<=s
                MF(i,1)=-(1/s)*input(i)+j;
                R(i,k)=j; k=k+1;
            end
        elseif j==MF_num
            if input(i)>=1-s && input(i)<=1
                MF(i,MF_num)=(1/s)*input(i)-(j-2);
                R(i,k)=j; 
            end
        else
            if input(i)>=(j-2)*s && input(i)<=(j-2)*s+2*s
                MF(i,j)=min(abs((1/s)*input(i)-(j-2)),abs((-1/s)*input(i)+j));
                R(i,k)=j; k=k+1;
            end
        end
    end
end
%%
M =[ 1     1     1     1
     1     1     1     2
     1     1     2     1
     1     1     2     2
     1     2     1     1
     1     2     1     2
     1     2     2     1
     1     2     2     2
     2     1     1     1
     2     1     1     2
     2     1     2     1
     2     1     2     2
     2     2     1     1
     2     2     1     2
     2     2     2     1
     2     2     2     2];

MR=zeros(2^N,N);
RN=zeros(2^N,1);
%%
for i=1:N
    MR(:,i)=R(i,M(:,i));
    RN(:,1)=RN(:,1)+(MR(:,i)-1)*MF_num^(N-i);
end
RN=RN+1;
idy=any(MR==0,1);
MR(:,idy)=[];
RN(:,idy)=[];
%%
MU=zeros(N,size(MR,1));
for i = 1:N
    MU(i,:) = MF(i,MR(:,i));
end

Phi_den=sum(prod(MU));
Phi_num_2=prod(MU);
Output.phi = Phi_num_2'/Phi_den;

Phi_num=sum(prod([MU;rule(RN,1)']));
Output.res = Phi_num/Phi_den;

Output.act=RN;
