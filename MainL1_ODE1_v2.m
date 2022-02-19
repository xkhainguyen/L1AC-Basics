clc
close all;
clear all;

params.dt = 0.0001; % plant interval
params.Ts = 0.01;   % adaptive sampling time
Tf = 0.1;           % simulation time
t = 0:params.dt:Tf;
ts = 0:params.Ts:Tf;

params.a = 1; params.am = params.a; params.b = params.a;    % system matrix
params.kr = params.am/params.b; params.r = 0;

% Initialization
u = zeros(1, length(ts));
x = zeros(1, length(t));
dx = zeros(1, length(t));
x_hat = zeros(1, length(ts));
dx_hat = zeros(1, length(ts));
x_tilde = zeros(1, length(ts));
d = zeros(1, length(t));
sigma_hat = zeros(1, length(ts));
x(1) = 0; x_hat(1) = 0; 

global u_prev
u_prev = 0;

% Piecewise-Constant Adaptive Laws Setting
params.Phi_Ts = (1 - exp(-params.am*params.Ts))/params.am;
params.tau = 0.1;

j = 1;

for i=1:length(t)   % Real time
    if (~mod(i-1,round(params.Ts/params.dt)))&&(i>0) %Sampling
        x_tilde(j) = x_hat(j) - x(i);

        sigma_hat(j) = estimator(x_tilde(j), params); % Adaptive Law

        if(j>0) 
            u(j) = controller(sigma_hat(j), params);   % Controller + Filter
        end
        
        if (j == length(ts))    % last time
            [~, dx_hat(j)] = predictor(x_hat(j), u(j), sigma_hat(j), params);    % Predicted system          
        else        
            [x_hat(j+1), dx_hat(j)] = predictor(x_hat(j), u(j), sigma_hat(j), params);    % Predicted system
            % Predictor uses different sampling time then x_tilde varies 
            j = j + 1;
        end
           
    end   
    
    if (i == length(t)) % last time
        [~, dx(i), d(i)] = plant(x(i), u(j), t(i), params);    % Real system
    else
        [x(i+1), dx(i), d(i)] = plant(x(i), u(j-1), t(i), params);    % Real system
    end
end


%%
plot(t,x); title('x');
figure
stairs(ts,sigma_hat); hold on ;plot(t,d,'--r'); title('Estimation');
figure
plot(ts,u); title('u control');
figure
plot(ts,x_tilde); title('x^~');
%%
figure
plot(ts,u);hold on;stairs(ts,sigma_hat); title('Filter');
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%% F U N C T I O N S %%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

function [x, dx, d] = plant(x_prev, u, t, params)
    d = 100*(step(t - 0.01) - step(t - 0.02))...
    + (300 + 500*sin(500*t))*(step(t - 0.04) - step(t - 0.05))...
    + (500*sin(1000*t) - 500)*(step(t - 0.07) - step(t - 0.08));

    dx = -params.a*x_prev + params.b*(u + d);
    
    x = x_prev + dx*params.dt;  
end

function [x_hat, dx_hat] = predictor(x_hat_prev, u, sigma_hat, params)
   dx_hat = -params.am*x_hat_prev + params.b*(u + sigma_hat);
   x_hat = x_hat_prev + dx_hat*params.Ts;
end

function u_f = controller(sigma_hat, params)
    global u_prev
    u = params.kr*params.r - sigma_hat;
    u_f = (params.Ts/params.tau)*u + (1-params.Ts/params.tau)*u_prev;
    u_prev = u_f;
end

function sigma_hat = estimator(x_tilde, params)
    muy = exp(-params.am*params.Ts)*x_tilde;
    sigma_hat = -params.Phi_Ts\muy/params.b;
end

function out = step(t)
    if t<0
        out = 0;
    else
        out = 1;
    end
end