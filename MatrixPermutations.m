%% Create permutations of matrix that satisfy inequalities.

%% Competition matrix

%     | T+  |  TP |  T- |       
% T+  | a11 | a12 | a13 |
% TP  | a21 | a22 | a23 |
% T-  | a31 | a32 | a33 | 


%% Inequalities 

% a31 > a21
% a32 > a12
% a13 > a23
% a13 > a12
% a23 > a21
% a32 > a31


% Create coefficient vector
coefficients = {'a12', 'a13', 'a21', 'a23', 'a31', 'a32'};

% Create all permutations of coefficients
allPossiblePermutations = perms(coefficients);

% Create empty vector for permutations that satisfy the inequalities.
satisfiedPermutations = {};

% Loop through possible permutations to find ordering that satify all six
% inequalities and append them to the satisfied permutations. 
for i = 1:1:length(allPossiblePermutations)
    
    disp(i)
    
    if(find(ismember(allPossiblePermutations(i,:),'a31')) > find(ismember(allPossiblePermutations(i,:),'a21')))
        
        if(find(ismember(allPossiblePermutations(i,:),'a32')) > find(ismember(allPossiblePermutations(i,:),'a12')))
            
            if(find(ismember(allPossiblePermutations(i,:),'a13')) > find(ismember(allPossiblePermutations(i,:),'a23')))
                
                if(find(ismember(allPossiblePermutations(i,:),'a13')) > find(ismember(allPossiblePermutations(i,:),'a12')))
                    
                    if(find(ismember(allPossiblePermutations(i,:),'a23')) > find(ismember(allPossiblePermutations(i,:),'a21')))
                        
                        if(find(ismember(allPossiblePermutations(i,:),'a32')) > find(ismember(allPossiblePermutations(i,:),'a31')))
                            
                            satisfiedPermutations(size(satisfiedPermutations,1) + 1, :) = allPossiblePermutations(i,:);
                            
                        end
                    end
                end
            end
        end
    end
        
end

%% Assign values from 0.4 - 0.9 to permutations.

values = [0.4 0.5 0.6 0.7 0.8 0.9];

for permutation = 1:size(satisfiedPermutations, 1)
        
    % Assign values for current order index
    for i = 1:1:6
        assignin('base', cell2mat(satisfiedPermutations(permutation,i)), values(i));
    end
    
    matrixCoefficients(permutation,:) = [a12 a13 a21 a23 a31 a32];
    
end


%% For Integrating evolutionary dynamics into treatment of metastatic 
%% castrate-resistant prostate cancer we use matrixCoefficients 7 for
%% Representative patient #1 and matrixCoefficients 5 for Representative 
%% patient #2



