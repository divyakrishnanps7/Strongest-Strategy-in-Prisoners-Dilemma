% param1 and param2 are arrays with,
% for agent 1 p1=param1(1), p2=param1(2), p3 = param1(3), p4 = param1(4)
% for agent 2 p1=param2(1), p2=param2(2), p3 = param2(3), p4 = param2(4)
%number of parameters
T=1.5;
R=1;
S=0;
P=0;
n_param = 4;
n_steps = 150;
n_iterations_avg = 50;
n_realisation_per_initial_condition =25;


% T1_values = [0,inf,0,inf,1,0,1,inf,1];
% T2_values = [0,0,inf,inf,0,1,inf,1,1];
T1_values = [0];
T2_values = [1];

s1_initial = [1,1,-1,-1];
s2_initial = [1,-1,1,-1];
% s1_initial = [1];
% s2_initial = [1];

total_T_pair = length(T1_values); %=length(T2_values)
total_initial_cond = length(s1_initial);

for k=1:total_T_pair
    
    p1_matrix = zeros(2^n_param);
    p2_matrix = p1_matrix ;
    T1 = T1_values(k);
    T2 = T2_values(k);
    
    
    for i=1:16
        %%initialisating parameter 1 here
        param1 = fliplr(de2bi((i-1),n_param));
        
        for j=1:16
            %% initialising parameter 2 her2
            param2 =  fliplr(de2bi((j-1),n_param));
            payoff1 = zeros([1, total_initial_cond*n_realisation_per_initial_condition]);
            payoff2 = zeros([1, total_initial_cond*n_realisation_per_initial_condition]);
            
            var1 = zeros([1, total_initial_cond*n_realisation_per_initial_condition]);
            var2 = zeros([1, total_initial_cond*n_realisation_per_initial_condition]);

            for i1=1:total_initial_cond
                %% doing initialisation here
                s1_init = s1_initial(i1);
                s2_init = s2_initial(i1);


                for i2 = 1:n_realisation_per_initial_condition
                    temp_payoff1 = zeros([1,n_steps]);
                    temp_payoff2 = zeros([1,n_steps]); 

                    Prob1= zeros([1, n_steps]);
                    Prob2 = zeros([1,n_steps]);
                        
                    s1 = zeros([1, n_steps]);
                    s2 = zeros([1, n_steps]);
                    s1(1) = s1_init;
                    s2(1) = s2_init;
                    
                    sa1 = zeros([1, n_steps]);
                    sa2 = zeros([1, n_steps]);
                    sm1 = zeros([1,n_steps]);
                    sm2 = zeros([1,n_steps]);
                    
                    i3=1;
                    sa1(i3)=s1(i3)*sign((1/(1+exp(-(1/T1))))-rand);
                    sa2(i3)=s2(i3)*sign((1/(1+exp(-(1/T1))))-rand);
                    sm2(i3)=sa2(i3)*sign((1/(1+exp(-(1/T2))))-rand); 
                    sm1(i3)=sa1(i3)*sign((1/(1+exp(-(1/T2))))-rand); 


                    for i3=2:n_steps

                        Prob1(i3)=((param1(1)+param1(4)-param1(2)-param1(3))/4)*(1+sa1(i3-1))*(1+sm2(i3-1))...
                            + ((param1(2)-param1(4))/2)*(1+sa1(i3-1)) + ...
                            ((param1(3)-param1(4))/2)*(1+sm2(i3-1)) + param1(4);
                        Prob2(i3)=((param2(1)+param2(4)-param2(2)-param2(3))/4)*(1+sm1(i3-1))*(1+sa2(i3-1)) +...
                            ((param2(2)-param2(4))/2)*(1+sa2(i3-1)) + ...
                            ((param2(3)-param2(4))/2)*(1+sm1(i3-1)) + param2(4);
                         s1(i3)=sign(Prob1(i3)-rand);
                         s2(i3)=sign(Prob2(i3)-rand); 
                         sa1(i3)=s1(i3)*sign((1/(1+exp(-(1/T1))))-rand);
                         sa2(i3)=s2(i3)*sign((1/(1+exp(-(1/T1))))-rand);
                         sm1(i3)=sa1(i3)*sign((1/(1+exp(-(1/T2))))-rand);
                         sm2(i3)=sa2(i3)*sign((1/(1+exp(-(1/T2))))-rand); 
                        %Payoff of agent 1
                        temp_payoff1(i3)=((R+P-S-T)/4)*(1+sa1(i3))*(1+sa2(i3)) + ((S-P)/2)*(1+sa1(i3)) + ((T-P)/2)*(1+sa2(i3)) + P;
                        %Payoff of agent 2
                        temp_payoff2(i3) = ((R+P-S-T)/4)*(1+sa1(i3))*(1+sa2(i3)) + ((S-P)/2)*(1+sa2(i3)) + ((T-P)/2)*(1+sa1(i3)) + P;

                    end

                payoff1(n_realisation_per_initial_condition*(i1-1)+i2) = mean(temp_payoff1(n_steps - n_iterations_avg: n_steps));
                var1((n_realisation_per_initial_condition*(i1-1)+i2)) = var((temp_payoff1(n_steps - n_iterations_avg: n_steps)));

                payoff2(n_realisation_per_initial_condition*(i1-1)+i2) = mean(temp_payoff2(n_steps - n_iterations_avg: n_steps));
                var2((n_realisation_per_initial_condition*(i1-1)+i2)) = var(temp_payoff2(n_steps - n_iterations_avg: n_steps));
                end
            end

            p1_matrix(i,j)=mean(payoff1);

            p2_matrix(i,j)=mean(payoff2);
        end
    end
    
    T1_matrix_cell{k}=p1_matrix;   
    T2_matrix_cell{k}=p2_matrix;

end

subplot(1,2,1)
imagesc(p1_matrix)
title('P1 Matrix')

subplot(1,2,2)
imagesc(p2_matrix)
title('P2 Matrix')