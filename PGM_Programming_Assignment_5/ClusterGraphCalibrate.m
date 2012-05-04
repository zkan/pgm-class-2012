% CLUSTERGRAPHCALIBRATE Loopy belief propagation for cluster graph calibration.
%   P = CLUSTERGRAPHCALIBRATE(P, useSmart) calibrates a given cluster graph, G,
%   and set of of factors, F. The function returns the final potentials for
%   each cluster. 
%   The cluster graph data structure has the following fields:
%   - .clusterList: a list of the cluster beliefs in this graph. These entries
%                   have the following subfields:
%     - .var:  indices of variables in the specified cluster
%     - .card: cardinality of variables in the specified cluster
%     - .val:  the cluster's beliefs about these variables
%   - .edges: A cluster adjacency matrix where edges(i,j)=1 implies clusters i
%             and j share an edge.
%  
%   UseSmart is an indicator variable that tells us whether to use the Naive or Smart
%   implementation of GetNextClusters for our message ordering
%
%   See also FACTORPRODUCT, FACTORMARGINALIZATION
%
% Copyright (C) Daphne Koller, Stanford University, 2012

function [P MESSAGES] = ClusterGraphCalibrate(P,useSmartMP)

if(~exist('useSmartMP','var'))
  useSmartMP = 0;
end

N = length(P.clusterList);

MESSAGES = repmat(struct('var', [], 'card', [], 'val', []), N, N);
[edgeFromIndx, edgeToIndx] = find(P.edges ~= 0);

for m = 1:length(edgeFromIndx),
    i = edgeFromIndx(m);
    j = edgeToIndx(m);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % YOUR CODE HERE
    %
    %
    %
    % Set the initial message values
    % MESSAGES(i,j) should be set to the initial value for the
    % message from cluster i to cluster j
    %
    % The matlab/octave functions 'intersect' and 'find' may
    % be useful here (for making your code faster)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    MESSAGES(i, j).var = intersect(P.clusterList(i).var, P.clusterList(j).var);
    MESSAGES(i, j).card = zeros(1, length(MESSAGES(i, j).var));

    % find the cardinality
    for k = 1:length(MESSAGES(i, j).var)
        for q = 1:length(P.clusterList)
            if(~isempty(find(P.clusterList(q).var == MESSAGES(i, j).var(k))))
                MESSAGES(i, j).card(k) = P.clusterList(q).card(find(P.clusterList(q).var == MESSAGES(i, j).var(k)));
                break;
            end
        end
    end

    % messages C_i receives except C_j
    factors_to_recv_msg = setdiff(find(P.edges(i, :) == 1), j);

    MESSAGES(i, j).val = ones(1, prod(MESSAGES(i, j).card));
    for v = 1:numel(factors_to_recv_msg)
        MESSAGES(i, j) = FactorProduct(MESSAGES(i, j), MESSAGES(factors_to_recv_msg(v), i));
    end
    MESSAGES(i, j) = FactorProduct(P.clusterList(i), MESSAGES(i, j));

    v_to_sum_out = setdiff(MESSAGES(i, j).var, P.clusterList(j).var);
    MESSAGES(i, j) = FactorMarginalization(MESSAGES(i, j), v_to_sum_out);
    
    % normalize the message
    MESSAGES(i, j).val = MESSAGES(i, j).val / sum(MESSAGES(i, j).val);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end;

% perform loopy belief propagation
tic;
iteration = 0;

lastMESSAGES = MESSAGES;

% this variable used for Quiz #1
iter = 0

while (1),
    iteration = iteration + 1;
    [i, j] = GetNextClusters(P, MESSAGES,lastMESSAGES, iteration, useSmartMP); 
    prevMessage = MESSAGES(i,j);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % YOUR CODE HERE
    % We have already selected a message to pass, \delta_ij.
    % Compute the message from clique i to clique j and put it
    % in MESSAGES(i,j)
    % Finally, normalize the message to prevent overflow
    %
    % The function 'setdiff' may be useful to help you
    % obtain some speedup in this function
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % no cluster sends a message to itself
    if i == j
        continue;
    end
    
    % messages C_i receives except C_j
    factors_to_recv_msg = setdiff(find(P.edges(i, :) == 1), j);
    
    MESSAGES(i, j).val = ones(1, prod(MESSAGES(i, j).card));
    for v = 1:numel(factors_to_recv_msg)
        MESSAGES(i, j) = FactorProduct(MESSAGES(i, j), MESSAGES(factors_to_recv_msg(v), i));
    end
    MESSAGES(i, j) = FactorProduct(P.clusterList(i), MESSAGES(i, j));

    v_to_sum_out = setdiff(MESSAGES(i, j).var, P.clusterList(j).var);
    MESSAGES(i, j) = FactorMarginalization(MESSAGES(i, j), v_to_sum_out);
    
    % normalize the message
    MESSAGES(i, j).val = MESSAGES(i, j).val / sum(MESSAGES(i, j).val);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if(useSmartMP==1)
      lastMESSAGES(i,j)=prevMessage;
    end

    % these 6 lines of code are used for Quiz #1
    if i == 19 && j == 3
        delta = MessageDelta(MESSAGES(i, j), lastMESSAGES(i, j))
        if delta >= 1.0e-6
            iter = iter + 1
        end
    end
    
    % Check for convergence every m iterations
    if mod(iteration, length(edgeFromIndx)) == 0
        if (CheckConvergence(MESSAGES, lastMESSAGES))
            break;
        end
        disp(['LBP Messages Passed: ', int2str(iteration), '...']);
        if(useSmartMP~=1)
          lastMESSAGES=MESSAGES;
        end
    end
    
end;
toc;
disp(['Total number of messages passed: ', num2str(iteration)]);


% Compute final potentials and place them in P
for m = 1:length(edgeFromIndx),
    j = edgeFromIndx(m);
    i = edgeToIndx(m);
    P.clusterList(i) = FactorProduct(P.clusterList(i), MESSAGES(j, i));
end


% Get the max difference between the marginal entries of 2 messages -------
function delta = MessageDelta(Mes1, Mes2)
delta = max(abs(Mes1.val - Mes2.val));
return;


