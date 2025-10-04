%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Supporting Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function description = get_exit_flag_description(exitflag)
    % Convert exit flag to human-readable description
    switch exitflag
        case 1
            description = 'Mesh tolerance reached';
        case 2
            description = 'Function tolerance reached';
        case 3
            description = 'Step tolerance reached';
        case 4
            description = 'Max iterations reached';
        case 0
            description = 'Max function evaluations reached';
        case -1
            description = 'Stopped by output/plot function';
        case -2
            description = 'No feasible point found';
        otherwise
            description = 'Unknown';
    end
end