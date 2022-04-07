%
% Purpose: Exercise 1 in Bogacz free energy tutorial
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


function ex_1

    close all;
    ex();
    
    % published solution for comparison
    % soln_ex_1();
end


function ex()

    % normal function
    gauss = @(x,mu,sigma) (1/(sqrt(2*pi*sigma)) *  exp( ((x-mu).^2)/(-2*sigma)) );

    % range of size of food item
    v_min = 0.01;
    v_step = 0.01;
    v_max =  5;
    v = v_min:v_step:v_max; 

    % likelihood p(u|v), u is the estimate of light intensity
    p_u_given_v = gauss(2,v.^2,1);

    % prior probability of observing this value of v
    p_v = gauss(v,3,1);

    % marginal likelihood
    p_u = sum(p_v .* p_u_given_v * 0.01);
    
    % apply Bayes Rule
    p_v_given_u = p_v.*p_u_given_v/p_u;

    hold on;
    plot(v,p_v_given_u);
    plot(v,p_v,'g');
    plot(v,p_u_given_v,'y');
    
    legend('p(v|u)','p(v)','p(u|v)');
    xlabel ("Size of food item, v");
        
end



% published solution for comparison
function soln_ex_1()

    v_p = 3;
    sigma_p = 1;
    sigma_u = 1; % mean of prior distribution of food size
    % standard deviation of prior = sqrt (1)
    % standard deviation of sensory noise = sqrt (1)
    u = 2; % observed light intensity
    MINV = 0.01;
    % minimum value of v for which posterior computed
    DV = 0.01;
    % interval between values of v for which posterior found
    MAXV = 5;
    % maximum value of v for which posterior computed
    vrange = [MINV:DV:MAXV ];
    numerator = normpdf (vrange ,v_p , sigma_p ) .* normpdf (u, vrange .^2 , sigma_u );
    normalization = sum ( numerator * DV );
    p = numerator / normalization ;
    plot (vrange , p, "k ");
    xlabel ("v ");
    ylabel ("p(v|u) ");
end