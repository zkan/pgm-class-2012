%CLIQUETREECALIBRATE Performs sum-product or max-product algorithm for 
%clique tree calibration.

%   P = CLIQUETREECALIBRATE(P, isMax) calibrates a given clique tree, P 
%   according to the value of isMax flag. If isMax is 1, it uses max-sum
%   message passing, otherwise uses sum-product. This function 
%   returns the clique tree where the .val for each clique in .cliqueList
%   is set to the final calibrated potentials.
%
% Copyright (C) Daphne Koller, Stanford University, 2012

function P = CliqueTreeCalibrate(P, isMax)

% Number of cliques in the tree.
N = length(P.cliqueList);

% Setting up the messages that will be passed.
% MESSAGES(i,j) represents the message going from clique i to clique j. 
MESSAGES = repmat(struct('var', [], 'card', [], 'val', []), N, N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% We have split the coding part for this function in two chunks with
% specific comments. This will make implementation much easier.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% YOUR CODE HERE
% While there are ready cliques to pass messages between, keep passing
% messages. Use GetNextCliques to find cliques to pass messages between.
% Once you have clique i that is ready to send message to clique
% j, compute the message and put it in MESSAGES(i,j).
% Remember that you only need an upward pass and a downward pass.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if isMax, work in log-space
if isMax
    for i = 1:N
        P.cliqueList(i).val = log(P.cliqueList(i).val);
    end
end

i = -1;
j = -1;
count = 1;
while i ~= 0 && j ~= 0
    [i, j] = GetNextCliques(P, MESSAGES);
    if i ~= 0 && j ~= 0

        MESSAGES(i, j).var = intersect(P.cliqueList(i).var, P.cliqueList(j).var);
        MESSAGES(i, j).card = zeros(1, length(MESSAGES(i, j).var));

        % find the cardinality
        for k = 1:length(MESSAGES(i, j).var)
            for q = 1:length(P.cliqueList)
                if(~isempty(find(P.cliqueList(q).var == MESSAGES(i, j).var(k))))
                    MESSAGES(i, j).card(k) = P.cliqueList(q).card(find(P.cliqueList(q).var == MESSAGES(i, j).var(k)));
                    break;
                end
            end
        end

        % check if it is a leaf node
        if numel(find(P.edges(i, :) == 1)) == 1
            if ~isMax
                v_to_sum_out = setdiff(P.cliqueList(i).var, MESSAGES(i, j).var);
                MESSAGES(i, j) = FactorMarginalization(P.cliqueList(i), v_to_sum_out);
            else
                v_to_max_out = setdiff(P.cliqueList(i).var, MESSAGES(i, j).var);
                MESSAGES(i, j) = FactorMaxMarginalization(P.cliqueList(i), v_to_max_out);
            end
        else
            % messages C_i receives except C_j
            factors_to_recv_msg = setdiff(find(P.edges(i, :) == 1), j);

            if ~isMax
                MESSAGES(i, j).val = ones(1, prod(MESSAGES(i, j).card));
                for v = 1:numel(factors_to_recv_msg)
                    MESSAGES(i, j) = FactorProduct(MESSAGES(i, j), MESSAGES(factors_to_recv_msg(v), i));
                end
                MESSAGES(i, j) = FactorProduct(P.cliqueList(i), MESSAGES(i, j));
            
                v_to_sum_out = setdiff(MESSAGES(i, j).var, P.cliqueList(j).var);
                MESSAGES(i, j) = FactorMarginalization(MESSAGES(i, j), v_to_sum_out);
            else
                MESSAGES(i, j).val = zeros(1, prod(MESSAGES(i, j).card));
                for v = 1:numel(factors_to_recv_msg)
                    MESSAGES(i, j) = FactorSum(MESSAGES(i, j), MESSAGES(factors_to_recv_msg(v), i));
                end
                MESSAGES(i, j) = FactorSum(P.cliqueList(i), MESSAGES(i, j));
            
                v_to_max_out = setdiff(MESSAGES(i, j).var, P.cliqueList(j).var);
                MESSAGES(i, j) = FactorMaxMarginalization(MESSAGES(i, j), v_to_max_out);
            end
        end

        if ~isMax
            % normalize the message
            MESSAGES(i, j).val = MESSAGES(i, j).val / sum(MESSAGES(i, j).val);
        end
    end

%    if i == 7 && j == 8
%        break;
%    end
end

%MESSAGES(2, 9)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% YOUR CODE HERE
%
% Now the clique tree has been calibrated. 
% Compute the final potentials for the cliques and place them in P.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:N
    for j = 1:N
        if P.edges(i, j) == 1
            if ~isMax
                P.cliqueList(i) = FactorProduct(P.cliqueList(i), MESSAGES(j, i));
            else
                P.cliqueList(i) = FactorSum(P.cliqueList(i), MESSAGES(j, i));
            end
        end
    end
end

return;
