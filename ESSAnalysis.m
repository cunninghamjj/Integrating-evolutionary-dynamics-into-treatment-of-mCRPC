%% Calculate ESS's for each competition matrix to initialize simulation parameters.

% Run MatrixPermutations.m to get possible coefficients.
MatrixPermutations;


for matrixIndex = 1:1:size(matrixCoefficients, 1)
    
    disp(matrixIndex)
    
    %% Solve ESS using ODE
    % Set growth rates fast so convergence happens quickly.
    r = [1, 1, 1];
    
    % Initialize populations at arbitrary start.
    y = [100, 100, 100, 0];
    
    % Create solutions matrix with initial values.
    allSolution(1,:) = y;
    
    % Set simulation time to very large number for convergence.
    maxSimulationTime = 100000;
    
    % Set sigma for PSA decay rate.
    sigmaPSA = 0.5;
    
    % Loop through simulation.
    for time = 2:maxSimulationTime
        
        % Update carrying capacities with current symbiotic T+.
        k = [y(2) * 1.5, 10000, 10000];
        
        % T+, TP, T-, and PSA ODE's
        dydt = zeros(1, 4);
        
        dydt(1) = y(1) * r(1) * (1 - ( ( y(1) + matrixCoefficients(matrixIndex,1) * y(2) + matrixCoefficients(matrixIndex,2) * y(3) ) / k(1) ) );
        dydt(2) = y(2) * r(2) * (1 - ( ( matrixCoefficients(matrixIndex,3) * y(1) + y(2) + matrixCoefficients(matrixIndex,4) * y(3) ) / k(2) ) );
        dydt(3) = y(3) * r(3) * (1 - ( ( matrixCoefficients(matrixIndex,5) * y(1) + matrixCoefficients(matrixIndex,6) * y(2) + y(3) ) / k(3) ) );
        dydt(4) = sum(y(1:3)) - sigmaPSA * y(4);
        
        y = y + dydt;
        
         % Add to solutions matrix.
        allSolution(time, :) = y;
        
    end
    
    % Save ESS population densities
    ESS(matrixIndex,:) = allSolution(end, :);
    
    % Save ESS population frequencies
    ESSFrequency(matrixIndex,:) = [allSolution(end,1)./( sum(allSolution(end,1:3), 2) ) .* 100, allSolution(end,2)./( sum(allSolution(end,1:3), 2) ) .* 100, allSolution(end,3)./( sum(allSolution(end,1:3), 2) ) .* 100];
    
end


%% For ESS population densities we'd like to keep a presence of cells in
%% so the minimum density is floored at 1E-9. 

ESS(ESS < 1E-9) = 1E-9;




