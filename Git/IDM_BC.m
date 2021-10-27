% 
% Copyright (c) 2021, CIMAT (www.cimat.com)
% All rights reserved. Please read the "license.txt" for license terms.
% 
% Project Title: Improvement Direction Mapping Method
% Publisher: CIMAT (www.cimat.mx)
% 
% Implemented by: S. Mostapha Kalami Heris, PhD
% 
% Contact Info: salvador.botello@cimat.mx
% 
% Base Reference Paper:
%
% 
% Reference Papaer URL:
% 

clear all;
close all;
clc;

%%
global maxeval;
global eval;

%%
[PF, mop] = SetMOP('ZDT6'); % Set MOP
func = @(x) mop.func(x')';
npop = 50;                  % Population size
maxeval = 10000;             % Max evaluations

%% Basic IDM Structure
evo = struct('npop', npop, 'nvar', mop.pd, 'nobj', mop.od, 'nweights', npop, 'func', func, ...
             'min_var', mop.domain(:, 1)', 'max_var', mop.domain(:, 2)', ...
             'dom', 'Aggregation', 'funcType', 'TCH', 'tensor', ["Broyden", ""], 'eta', 0.0);

%%  IDM method: Chebychev-Broyden
eval = 0;
evo = IDM_Full_DO(evo, 'init');                 % Initialize IDM method
for i = 1:evo.npop                              % Evaluate initial population
    evo.f(i,:) = mop.func(evo.x(i, :)')';
    eval = eval + 1;
end
evo.q = [];
evo.g = [];
evo = IDM_Full_DO(evo, 'update');                % Update IDM structure
gen = 1;                                         % Initialize generation counter

% Show Iteration Information
fprintf("Algorithm: Chebychev-Broyden, Problem: %s, Eval: %d, Gen: %d\n", mop.name, eval, gen);

% IDM main loop
while (eval < maxeval)
    evo = IDM_Full_DO(evo, 'step');              % Offspring generation/Offspring evaluation
    evo = IDM_Full_DO(evo, 'update');            % Selection process/Update IDM structure
    gen = gen + 1;                               % Update generation counter
    
    % Show Iteration Information
    fprintf("Algorithm: Chebychev-Broyden, Problem: %s, Eval: %d, Gen: %d\n", mop.name, eval, gen);
end

fprintf('Optimization Terminated.\n');

% Plot final front
plot(PF(:, 1), PF(:, 2), '-r', evo.f(:, 1), evo.f(:, 2), '.');