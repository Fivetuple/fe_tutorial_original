%
% Purpose: Exercise 2 in Bogacz free energy tutorial
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


function ex_2

    ex();
   
    % published solution for comparison
    %soln_ex_2();
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
   
    % phi is the most likely size of food item
    phi = zeros(iterations,1);
    phi(1) = 3;
    
    for i=2:iterations
        
        phi_last = phi(i-1,1);
        term1 = (v_p - phi_last)/sigma_p;
        term2 = (u - phi_last^2) * (2*phi_last)/sigma_u;
        
        % see eqn. (8)
        gradient = term1 + term2;
        
        phi(i,1) = phi_last + dt*gradient;
    end
    plot(1:iterations,phi);
    xlabel ('Time step');
    ylabel ('phi ');
end

% published solution for comparison
function soln_ex_2()
    v_p = 3;
    Sigma_p = 1;
    Sigma_u = 1; % mean of prior distribution of food size
    % variance of prior distribution
    % variance of sensory noise
    u = 2; % observed light intensity
    DT = 0.01;
    MAXT = 5; % integration step
    % maximum time considered
    phi (1) = v_p; % initializing the best guess of food size
    for i = 2: MAXT/DT
    phi(i) = phi(i -1) + DT * (( v_p - phi(i -1))/ Sigma_p + ...
    (u-phi(i -1)^2)/ Sigma_u * (2* phi(i -1)));
    end
    plot ([ DT:DT:MAXT], phi , "k ");
    xlabel ("Time ");
    ylabel ("\phi ");
    axis ([0 MAXT -2 3.5]);
 
end