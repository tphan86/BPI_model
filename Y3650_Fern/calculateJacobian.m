function J = calculateJacobian(modelFunc, params, time)
    numParams = length(params);
    numTimePoints = length(time);
    perturbation = 1e-2;
    J = zeros(numTimePoints, numParams);
    
    for j = 1:numParams

        % Perturb the parameter
        paramsPerturbed = params;

        if (j ~= 12) || (params(12) == 0)
            paramsPerturbed(j) = paramsPerturbed(j) + perturbation;
        else
            paramsPerturbed(j) = paramsPerturbed(j) + perturbation*params(j);
        end
        
        % Calculate the model output with perturbed parameters
        [T_perturbed, ~, ~, ~] = modelFunc(paramsPerturbed, time);
        
        % Calculate the model output with original parameters
        [T, ~, ~, ~] = modelFunc(params, time);
        
        % Compute the Jacobian (sensitivity) for this parameter
        J(:, j) = (T_perturbed - T) / perturbation;
    end
end