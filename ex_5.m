%
% Purpose: Exercise 5 in Bogacz free energy tutorial
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

function ex_5

    figure('Name','ex');
    ex();

    % published solution for comparison
    % figure('Name','soln');
    % soln();
end


function ex()

    % alpha is the learning rate
    alpha = 0.01;
    
    % define mean and variance for phi_i
    mean = 5;
    variance = 2;
        
    % g is g_i(phi_i+1) in eqn (59)
    g = 5;
    
    rng(0);
    
    DT = 0.01;
    iterations = floor(20/DT);
    n_trials = 2000;
    
    sigma = 1*ones(n_trials,1);
    
    for trial = 2:n_trials

        e = ones(iterations,1);
        eps = ones(iterations,1);
    
        phi = normrnd(mean,sqrt(variance));

        sigma_last = sigma(trial-1,1);

        for i=2:iterations
      
            e_last = e(i-1,1);
            eps_last = eps(i-1,1);
            
            % eqns (59) and (60)
            eps(i,1) = eps_last +  DT*(phi - g - e_last);
            e(i,1) = e_last + DT*(sigma_last * eps_last - e_last);
        end     
        
        sigma(trial,1) = sigma_last + alpha*(eps(i,1) * e(i,1) - 1);
    end
    plot(1:n_trials,sigma);
    
    xlabel ("Trial number ");
    ylabel ("Sigma ");
end



% published solution for comparison
function soln()

    mean_phi = 5;
    Sigma_phi = 2;
    phi_above = 5; % mean of input from the current level
    % variance of input from the current level
    % input from the level above
    DT = 0.01;      %  integration step
    MAXT = 20;      %  maximum time considered
    TRIALS = 2000;  %  number of simulated trials
    LRATE = 0.01;   %  learning rate
    
    Sigma (1) = 1; % initializing the value of weight
    
    for trial = 2: TRIALS
        error (1) = 0; %#ok<*AGROW>
        % initializing the prediction error
        e(1) = 0;

        % initializing the interneuron
        phi = mean_phi + sqrt( Sigma_phi ) * randn ;

        for i = 2: MAXT/DT
            error (i) = error (i -1) + DT * (phi - phi_above - e(i -1));
            e(i) = e(i -1) + DT * ( Sigma (trial -1) * error (i -1) - e(i -1));
        end
        Sigma( trial ) = Sigma (trial -1) + LRATE * ( error (end )*e(end) - 1); 
    end
    plot (Sigma , "k ");
    xlabel ("Trial ");
    ylabel ("\Sigma ");
   
end