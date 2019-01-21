%% simple model contrasting mental rotation vs. response substitution in M1
%% see McDougle & Taylor, Nat. Comm. (2019), Figure S6 for details
%% inspired by Georgopoulos et al. (1989) and Cisek & Scott (1999)
%% 9/13/17; NYC
clear;clc;close all;

num_neurons = 300; % number of simulated units
theta = linspace(0,pi*2,num_neurons); % prefered tuning directions
targets = 0; % target position(s)
x_rs = zeros(num_neurons,length(targets),500); % processed perceptual input x; response subtitution model
y_rs = zeros(num_neurons,length(targets),500); % activity y; response subtitution model
x_mr = zeros(num_neurons,length(targets),500); % processed perceptualinput x; mental rotation model
y_mr = zeros(num_neurons,length(targets),500); % activity y; mental rotation model
gamma = 1.3; % time constant on activity change
delay_init = linspace(50,150,num_neurons); % uniform dist. of neuron response onset delays (in ms)
delay_init = delay_init(randperm(length(delay_init))); % shuffle

a = [ones(1,150) linspace(1,0,200) zeros(1,150);... % stim input vector for rs (linear decrease of target-aimed vector)
    zeros(1,150) linspace(0,1,200) ones(1,150)];    % goal input vector (linear increase of goal-aimed vector)

drive = 1; % fixed gain parameter
a_mr = [ones(1,150) linspace(drive,1,200) ones(1,150)]; % steadily fixed input (angle changes uniformly; "mental rotation")

iterations = 30; 
for iter = 1:iterations
    for i = 1:num_neurons
        
        %% set neuron's params (values from Cisek & Scott 1999)
        s = normrnd(0,.3); % offset of cosine tuning function for neuron i
        Gamma = normrnd(.4,.2); % activity threshold parameter for neuron i
        gain = normrnd(1.2,.5); % gain on excitation for neuron i
        delay = delay_init(i); % onset delay for neuron i (from pre-specified uniform dist)
        
        %% loop over targets
        for T = 1:1 % single target case; can sim with n targets
            
            %% input
            phi = [targets(T) targets(T)+(pi/2)]; % target and goal
            e(1,:) = a(1,:) .* (cos(theta(i) - phi(1)) + s); % excitation for target
            e(2,:) = a(2,:) .* (cos(theta(i) - phi(2)) + s); % excitation for goal
            
            %% avoid negative excitation
            id_neg_act = find(e(1,:) < 0);
            e(1,id_neg_act) = 0;
            id_neg_act2 = find(e(1,:) < 0);
            e(2,id_neg_act2) = 0;
            
            %% sum excitation vectors for target and goal (response substitution)
            E_rs(i,T,:) = e(1,:) + e(2,:);
            
            %% excitation for mental rotation
            mr_input = [zeros(1,150), linspace(phi(1),phi(2),200), phi(2)*ones(1,150)];
            E_mr(i,T,:) =  a_mr .* (cos(theta(i) - mr_input) + s);
            id_neg_act = find(E_mr(i,T,:) < 0);
            E_mr(i,T,id_neg_act) = 0;
            
            %% loop over time points (ms)
            for t = 1:length(a(1,:))
                if t > delay
                    
                    %% response substitution
                    x_rs(i,T,t+1) = x_rs(i,T,t) + gamma * (-x_rs(i,T,t) + gain*E_rs(i,T,t)); % processed input to cell i
                    y_rs(i,T,t) = x_rs(i,T,t) - Gamma; % activity of cell i
                    %% mental rotation
                    x_mr(i,T,t+1) = x_mr(i,T,t) + gamma * (-x_mr(i,T,t) + gain*E_mr(i,T,t)); % processed input to cell i
                    y_mr(i,T,t) = x_mr(i,T,t) - Gamma; % activity of cell i
                    %% no negative activity
                    if y_rs(i,T,t) < 0
                        y_rs(i,T,t) = 0; 
                    end
                    %% no negative activity
                    if y_mr(i,T,t) < 0
                        y_mr(i,T,t) = 0;
                    end
                end
            end
        end
        %% store values for population vector
        % response substitution
        weight_rs(i,:) = squeeze(y_rs(i,1,:)); % activity for neuron i
        unit_vector = [cos(theta(i)) sin(theta(i))]; % unit vector for neuron i
        px_rs(i,:) = weight_rs(i,:)*unit_vector(1); % Cartesian x-value for neuron i
        py_rs(i,:) = weight_rs(i,:)*unit_vector(2); % Cartesian y-value for neuron i
        % mental rotation
        weight_mr(i,:) = squeeze(y_mr(i,1,:)); % activity for neuron i
        unit_vector = [cos(theta(i)) sin(theta(i))]; % "unit vector" for neuron i
        px_mr(i,:) = weight_mr(i,:)*unit_vector(1); % Cartesian x-value for neuron i
        py_mr(i,:) = weight_mr(i,:)*unit_vector(2); % Cartesian y-value for neuron i
    end
    %% compute population vector
    %% response substitution
    Px_rs = sum(px_rs); % pop vector X
    Py_rs = sum(py_rs); % pop vector Y
    PA_rs(iter,:) = atan2(Py_rs,Px_rs)*180/pi; % pop vector ANGLE
    PL_rs(iter,:) = sqrt(Px_rs.^2 + Py_rs.^2); % pop vector LENGTH
    % mental rotation
    Px_mr = sum(px_mr); % pop vector X
    Py_mr = sum(py_mr); % pop vector Y
    PA_mr(iter,:) = atan2(Py_mr,Px_mr)*180/pi; % pop vector ANGLE
    PL_mr(iter,:) = sqrt(Px_mr.^2 + Py_mr.^2); % pop vector LENGTH
end

%% PLOTTING %%
%% plot angle over t
figure;hold on;
colors = {[.7 .1 .5],[.1 .7 .2]};
%% all iterations
for g = 1:iterations
    plot(PA_rs(g,:),'linewidth',1,'color',colors{1});
    plot(PA_mr(g,:),'linewidth',1,'color',colors{2});
end
legend('Substitution','Rotation','location','southeast');
xlabel('Time (ms)');
ylabel('Population vector angle (deg)');
axis([150 400 -30 100]);
box off


%% plot velocity over t
figure;hold on;
eta = .76; % velocity constant (arbitrary scaling parameter)
velocity_rs = PL_rs .* eta;
velocity_mr = PL_mr .* eta;
interval = 1:400;

for g = 1:iterations
plot(interval,velocity_rs(g,interval),'color',colors{1});
plot(interval,velocity_mr(g,interval),'color',colors{2});
end

h=legend('Substitution','Rotation','location','northeast');
set(h,'linewidth',1);
xlabel('Time (ms)');
ylabel('Movement speed (cm/s)');
box off
axis([148 402 20 63]);



