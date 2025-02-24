% Econ_PE1

clc, clear all, close all

%%
% m = the number of non-capital goods
% n = the number of processes
% T = the number of simulation times
% dt = the time step
% q = the variable to reduce the step size when there is an overflow

% pi = n*(m+1) matrix characterizing the model economy
% c = 1*n matrix characterizing the model economy

% t(index) = the time
% r(index,i) = the intensity of process i at time index
% p(index,j) = the price of good j at time index
% e0(index) = the excess demand of the capital good at time index
% e(index,j) = the excess demand of the non-capital good j at time index
% C(index) = the total amount of capital in the economy at time index
% G(index) = the relative growth rate of the economy at time index
% D(index) = the diversity index of the economy at time index

% alpha(i) = the fraction of capital held by process i
% S = the entropy of the economy
%%

% Parameters
m = 100;
n = 1000;
T = 100000;
dt = 0.0001;
q = 1;

% load the fixed pi and c
% we change the varance of c by loading different variables in generator.mat
pi = cell2mat(struct2cell(load("generator.mat","pi")));
c = cell2mat(struct2cell(load("generator.mat","c_1over3")));

% initialize all the varibles interested
t = zeros(1,1);
r = zeros(1, n);
p = zeros(1, m);
e0 = zeros(1,1);
e = zeros(1, m);
C = zeros(1,1);
G = zeros(1,1);
D = zeros(1,1);

% get the intensity at time 0
pi_new = pi(:,2:m+1);  % drop the column of capital good
u = ones(1,n);
H = eye(n);
f = zeros(1,n);
Aeq = pi_new';
Beq = -pi_new'*u';
lb = (zeros(1,n)-u)';
w = quadprog(H,f,[],[],Aeq,Beq,lb,[]);
r_t = w'+u;

% perform the simulation
for index = 1:T
    t(index) = dt*index;
    
    % clear the variables of the previous time
    A = zeros(m,m); % A is calculated for minimizing phi
    b = zeros(m,1); % b is calculated to calculate f for minimizing phi
    e_t = zeros(m,1); 
    f = zeros(m,1); % f is calculated for minimizing phi
    alpha = zeros(1,n); 

    for j = 1:m
        for k = 1:m
            for i = 1:n
                % calculate A
                A(j,k) = A(j,k) + r_t(i)/c(i)*pi(i,j+1)*pi(i,k+1);
            end
        end
        
        for i = 1:n
            % calculate b
            b(j)= b(j) + r_t(i)/c(i)*pi(i,1)*pi(i,j+1);
            % calculate e_t
            e_t(j) = e_t(j) - r_t(i)*pi(i,j+1);
        end
        % calculate f
        f(j) = b(j)-e_t(j)/dt;
    end   
    % phi = 1/2 * (p_t.') * A * p_t + (p_t.') * (b - 1/dt * e);
    % we use quadprog as followed to minimize phi
    p_t = quadprog(A,f,[],[],[],[],zeros(m,1),[]);

    % calculate e
    e(index,:) = e_t;

    % update p
    p(index,:) = p_t;

    % update r
    for i = 1:n
        r_temp = 0;
        for j = 1:m
            r_temp = r_temp + pi(i,j+1)*p_t(j);
        end
        r_t(i) = r_t(i) * (1 + dt/c(i) * (pi(i,1) + r_temp));
    end
    r(index,:) = r_t;

    % calculate C and G
    dCdt =0;
    C(index) = 0;
    for i = 1:n
        dCdt_temp = 0;
        for j = 1:m
            dCdt_temp = dCdt_temp + pi(i,j+1)*p_t(j);
        end
        dCdt = dCdt + r_t(i)*(dCdt_temp +pi(i,1));
        C(index) = C(index) + r_t(i)*c(i);
    end
    G(index) = dCdt/C(index);
    G(index) = G(index) * q;
    
    % check and update
    if index > 1
        G_diff = G(index) - G(index-1);
        if G_diff > 1   % the signal for overflow
            pi = pi/2;
            q = q*2;   % update q

            % replace the wrong results of this iteration with the previous ones
            r_t = r(index-1,:);
            p_t = p(index-1,:);
            r(index,:) = r_t;
            p(index,:) = p_t;
            G(index) = G(index-1);
        end
    end


    % calculate alpha
    for i = 1:n
        alpha(i) = r_t(i)*c(i)/C(index);
    end
    
    % calculate S
    S_max =log(n);
    S = 0;
    for i =1:n
        S =  S-alpha(i)*log(alpha(i));
    end

    % calculate D
    D(index) = S/S_max;
   
    % check e0
    e0(index) = 0;
    for i = 1:n
        e0_temp = 0;
        for k = 1:m
            e0_temp = e0_temp + pi(i,k+1)*p_t(k);
        end
        e0(index) = e0(index) + r_t(i) * e0_temp;
    end
    
    % check e
    for j=1:m
        e(index,j) = 0;
        for i=1:n
            e(index,j) = e(index,j) - r_t(i)*pi(i,j+1);
        end
    end


end

% save the data
save("PE1_c_var_1over3.mat",'r','c','pi','p','e0','D','G','C','e','t');
