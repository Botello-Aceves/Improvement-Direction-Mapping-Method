%% JacobianDirectionalBias Object Data
% evo.npop
% evo.nvar
% evo.nobj
% evo.sd
% evo.rho

% evo.x 
% evo.f

% evo.func
% evo.df

% evo.idx
% evo.namda
% evo.a

% evo.min_var;
% evo.max_var;

%%
function evo = DS_Do(evo, do_string)
    switch do_string
        case 'init'
            evo = Aggregation_IDM_init(evo);
        case 'update'
            evo = Historic_Update(evo);
        case 'step'
            evo = Historic_step(evo);
        otherwise
            fprintf('Unknown command.\n');
    end
end

%%
function evo = Aggregation_IDM_init(evo)
    evo.f = zeros(evo.npop, evo.nobj);
    evo.x = rand(evo.npop, evo.nvar);
    evo.x = evo.x.*repmat(evo.max_var - evo.min_var, evo.npop, 1) + repmat(evo.min_var, evo.npop, 1);
    
    
    if(evo.nobj > 2)
        evo = Aggregation_IDM_weight_init_MD(evo);
    else
        evo = Aggregation_IDM_weight_init_2D(evo);
    end
    
    evo.H  = rand(evo.nvar, evo.nobj, evo.npop);
    evo.dx = zeros(evo.npop, evo.nvar);
    evo.dz = zeros(evo.npop, evo.nobj);
    evo.a  = ones(1, evo.npop);
end

%%
function evo = Aggregation_IDM_weight_init_2D(evo)
    evo.w = zeros(evo.npop, 2);
    evo.w(:, 1) = 0.01:(0.98/(evo.npop-1)):0.99;
    evo.w(:, 2) = 1 - evo.w(:, 1);
end

%%
function evo = Aggregation_IDM_weight_init_MD(evo)
    for i = 1:evo.npop
        w = rand(1, evo.nobj); w = w/sum(w);
%         w = [w - [0, w(1:end-1)], 1-w(end)];
        evo.w(i, :) = w;
    end
    
    t = evo.npop;
    t_max = 1e4 + evo.npop;
    
    dist = zeros(1, evo.npop + 1);
%     hv = zeros(evo.npop + 1, evo.nobj);
    while(t < t_max)
        w = rand(1, evo.nobj); w = w/sum(w);
%         w = [w - [0, w(1:end-1)], 1-w(end)];
        evo.w = [evo.w; w];
        
        for i = 1:size(evo.w, 1)
            d = vecnorm(repmat(evo.w(i, :), size(evo.w, 1), 1) - evo.w, 2, 2);
            [~, idx] = sort(d);
            dist(i) = sum(d(idx(1:evo.nobj+1)));
        end
        
%         hv_all = approximate_hypervolume_ms(evo.w', ones(evo.nobj, 1)*2);
%         for i = 1:size(evo.w, 1)
%             W = evo.w; W(i, :) = [];
%             hv(i) = hv_all - approximate_hypervolume_ms(W', ones(evo.nobj, 1)*2);
%         end
        [~, id]= min(dist);
        evo.w(id(1), :) = [];
        t = t +1;
    end
    plot3(evo.w(:, 1), evo.w(:, 2), evo.w(:, 3), '.'); pause(0.1);
end

%%
function evo = Historic_SearchDirection(evo)
    evo.dz = -evo.w;
end

%%
function evo = Historic_step(evo)
    global eval;
    for i = 1:evo.npop
        dx_ = evo.H(:, :, i)*evo.dz(i, :)';
        dx = evo.eta*evo.dx(i, :) + (1 - evo.eta)*dx_';
        a = evo.a(i);
        
        if (a < 1e-6)
            a = 1.0;
        end
        
        a_ = 0;
        da = a;
        k  = 0;
%         while(da > 1e-6 && k < 3)
            evo.q(i, :) = evo.x(i, :) + (a_ + da)*dx;
            evo.q(i, :) = max([evo.q(i, :); evo.min_var]);
            evo.q(i, :) = min([evo.q(i, :); evo.max_var]);
            evo.g(i, :) = evo.func(evo.q(i, :)); eval = eval+1;
            
            if(dominates(evo.g(i, :), evo.f(i, :)))
                evo.a = 1.2*evo.a;
            else
                evo.a = 0.8*evo.a;
            end
            k = k + 1;
%         end
        evo.dx(i, :) = dx;
    end
end

%%
function evo = Historic_Update(evo)
    X = [evo.x; evo.q];
    F = [evo.f; evo.g];
    
   if(~isempty(evo.g))
        for i = 1:evo.npop
            if(dominates(evo.g(i, :), evo.f(i, :)))
                evo.x(i, :) = evo.q(i, :);
                evo.f(i, :) = evo.g(i, :);
            end
        end
   end
    
    evo = Tensor_IDM(evo);
    evo = Historic_SearchDirection(evo);
end

%%
function J = PDS_NumJacobian(x, evo)
    n = length(x); 
    s = x + 1e-6*rand(n+1, n);
    for i = 1:n+1
        z(i, :) = evo.func(s(i, :));
    end
    
    A = MultivarLinearRegression(s, z);
    [~, J] = LinearFunction(A, x);
end