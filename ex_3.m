%
% Purpose: Exercise 3 in Bogacz free energy tutorial
% https://www.sciencedirect.com/science/article/pii/S0022249615000759
%
% The task considers a simple perceptual problem in which a value of 
% a single variable has to be inferred from a single observation. 
% We consider a simple organism that tries to infer the size or diameter 
% of a food item, which we denote by v , on the basis of light intensity 
% it observes.  Please see the paper for details of the exercise.
%
% (c) 2022 Paul Moore - moorep@maths.ox.ac.uk 
%
% This software is provided 'as is' with no warranty or other guarantee of
% fitness for the user's purpose.  Please let the author know of any bugs
% or potential improvements.


function ex_3

    ex();
    
    % published solution for comparison
    %soln();
end


function ex()

    % time step    
    dt = 0.01; 

    % observed light intensity mean and variance
    u = 2; 
    sigma_u = 1;  
    
    % size expectation and variance
    v_p = 3; 
    sigma_p = 1;  
    
    iterations = floor(5/dt);
    
    phi = zeros(iterations,1);
    e_p = zeros(iterations,1);
    e_u = zeros(iterations,1);
    
    % initialise values
    phi(1,1) = v_p;    
    e_p(1,1) = 0;
    e_u(1,1) = 0;
    
    for i=2:iterations
        
        phi_last = phi(i-1,1);
        e_p_last = e_p(i-1,1);
        e_u_last = e_u(i-1,1);
        
        % eqns (13) and (14)
        de_p = phi_last - v_p - (sigma_p * e_p_last);
        de_u = u - phi_last^2 - (sigma_u * e_u_last);
        
        gradient = e_u_last * 2*phi_last - e_p_last;
        
        phi(i,1) = phi_last + dt*gradient;
        e_p(i,1) = e_p_last + dt*de_p;
        e_u(i,1) = e_u_last + dt*de_u;
    end
    hold on;
    plot(1:iterations,phi);
    plot(1:iterations,e_p,'r');
    plot(1:iterations,e_u,'g');
    
    xlabel ('Time step');
    ylabel ('Value');
    legend ("phi ", "eps_p ", "eps_u ");
end


% published solution for comparison
function soln()

    v_p = 3;
    Sigma_p = 1;
    Sigma_u = 1; % mean of prior distribution of food size
    % variance of prior distribution
    % variance of sensory noise
    u = 2; % observed light intensity
    DT = 0.01;
    MAXT = 5; % integration step
    % maximum time considered
    phi (1) = v_p;
    % initializing the best guess of food size
    error_p (1) = 0; % initializing the prediction error of food size
    error_u (1) = 0; % initializing the prediction error of sensory input
    for i = 2: MAXT/DT
    phi(i) = phi(i -1) + DT * (- error_p (i -1) + error_u (i -1) * (2* phi(i -1)));
    error_p (i) = error_p (i -1) + DT * (phi(i -1) - v_p - Sigma_p * error_p (i -1));
    error_u (i) = error_u (i -1) + DT * (u - phi(i -1)^2 - Sigma_u * error_u (i -1));
    end
    plot ([ DT:DT:MAXT], phi , "k ");
    hold on;
    plot ([ DT:DT:MAXT], error_p , "k--");
    plot ([ DT:DT:MAXT], error_u , "k: ");
    xlabel ("Time ");
    ylabel ("Activity ");
    legend ("\phi ", "\epsilon_p ", "\epsilon_u ");
    axis ([0 MAXT -2 3.5]);

end