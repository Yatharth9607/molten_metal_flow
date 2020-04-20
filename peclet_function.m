% peclet function (A(|P|))

% parameters required:
% P: Peclet number (array)
% solver_type: type of solver ('central_diff', 'exponential', 'power_law')

% return output:
% P: A(|P|) depending upon the solver function

function [P] = peclet_function(P, solver_type)
    
    % central difference scheme
    if solver_type == "central_diff"
        P = 1 - 0.5 * abs(P);
    
    % exponential scheme
    elseif solver_type == "exponential"
        if abs(P) < 700
            P = exp(abs(P))./ (exp(abs(P)) - 1);
        else
            P = 1;
        end

    % power law
    elseif solver_type == "power_law"
        mult = (1 - 0.1 * abs(P));
        P = max(zeros(size(P, 1), 1), mult.*mult.*mult.*mult.*mult);

    else
        error("Please enter correct solver type")
    
    end

end
    