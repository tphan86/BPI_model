%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Supporting Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function stop = optimization_output(x, optimValues, state, param_names)
    % Custom output function for monitoring optimization progress
    persistent best_fval iteration_count;
    
    stop = false;
    
    switch state
        case 'init'
            fprintf('Initializing optimization...\n');
            best_fval = Inf;
            iteration_count = 0;
            
        case 'iter'
            iteration_count = iteration_count + 1;
            
            % Track best solution
            if optimValues.fval < best_fval
                best_fval = optimValues.fval;
            end
            
            % Periodic detailed reporting
            if mod(iteration_count, 25) == 0
                fprintf('  Iteration %d: Current = %.4e, Best = %.4e\n', ...
                    iteration_count, optimValues.fval, best_fval);
            end
            
        case 'done'
            fprintf('Optimization run completed after %d iterations.\n', iteration_count);
    end
end