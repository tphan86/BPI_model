%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Supporting Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function obj_val = robust_objective_function(para_fit, data_mean, time)
    % Enhanced objective function with comprehensive error handling
    try
        % Check for valid parameter values
        if any(~isfinite(para_fit)) || any(imag(para_fit) ~= 0)
            obj_val = 1e12;
            return;
        end
        
        % Evaluate model
        model_output = total_host_cells(para_fit, time);
        
        % Check for valid model outputs
        if any(~isfinite(model_output)) || any(imag(model_output) ~= 0) || any(model_output < 0)
            obj_val = 1e12;
            return;
        end
        
        % Calculate objective (norm of residuals)
        residuals = data_mean - model_output;
        obj_val = norm(residuals);
        
        % Optional: Add small regularization to prevent parameter explosion
        reg_weight = 1e-8;
        reg_term = reg_weight * norm(para_fit);
        obj_val = obj_val + reg_term;
        
    catch ME
        % Handle any unexpected errors
        warning('Error in model evaluation: %s', ME.message);
        obj_val = 1e12;
    end
end