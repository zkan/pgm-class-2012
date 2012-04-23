%GETNEXTCLIQUES Find a pair of cliques ready for message passing
%   [i, j] = GETNEXTCLIQUES(P, messages) finds ready cliques in a given
%   clique tree, P, and a matrix of current messages. Returns indices i and j
%   such that clique i is ready to transmit a message to clique j.
%
%   We are doing clique tree message passing, so
%   do not return (i,j) if clique i has already passed a message to clique j.
%
%	 messages is a n x n matrix of passed messages, where messages(i,j)
% 	 represents the message going from clique i to clique j. 
%   This matrix is initialized in CliqueTreeCalibrate as such:
%      MESSAGES = repmat(struct('var', [], 'card', [], 'val', []), N, N);
%
%   If more than one message is ready to be transmitted, return 
%   the pair (i,j) that is numerically smallest. If you use an outer
%   for loop over i and an inner for loop over j, breaking when you find a 
%   ready pair of cliques, you will get the right answer.
%
%   If no such cliques exist, returns i = j = 0.
%
%   See also CLIQUETREECALIBRATE
%
% Copyright (C) Daphne Koller, Stanford University, 2012


function [i, j] = GetNextCliques(P, messages)

% initialization
% you should set them to the correct values in your code
i = 0;
j = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% YOUR CODE HERE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = size(messages, 1);
%message_passing_matrix = zeros(n);

% check upstream
for i = 1:n
    for j = i:n
        if P.edges(i, j)
%            printf('i %d -> j %d: %d\n', i, j, ~isempty(messages(i, j).var));
            if isempty(messages(i, j).var)

                % check if it is already received messages from its neighbors?
%                neighbors = find(P.edges(i, :) == 1);
%                for k = 1:length(neightbors)
%                    if neighbors(k) ~= j
%                        if(~isempty(message(i, j).var))                            
%                        end
%                    end
%                end
%                printf('upstream\n');
                return;
            end
        end
    end
end

%printf('\n');

% check downstream (if we pass upstream, all messages are sent to the root already.)
for j = n:-1:1
    for i = 1:n
        if P.edges(j, i)
%            printf('j %d -> i %d: %d\n', j, i, ~isempty(messages(j, i).var));
            if isempty(messages(j, i).var)
%                message_passing_matrix(j, i) = 1;
                tmp = i;
                i = j;
                j = tmp;
%                printf('downstream\n');
                return;
            end
        end
    end
end

%P.edges
%message_passing_matrix

% no such cliques exist, returns i = j = 0.
i = 0;
j = 0;

return;
