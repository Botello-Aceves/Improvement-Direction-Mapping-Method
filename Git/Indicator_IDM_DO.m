function evo = Indicator_IDM_DO(evo, do_string)
    switch do_string
        case 'init'
            evo = Indicator_IDM_init(evo);
        case 'directions'
            evo = Indicator_IDM_dir(evo);
        case 'step'
            evo = Indicator_IDM_Step(evo);
        case 'update'
            evo = Indicator_IDM_Update(evo);
        otherwise
            
    end
end

%%
function evo = Indicator_IDM_init(evo)
    evo.x = rand(evo.npop, evo.nvar);
    evo.x = evo.x.*repmat(evo.max_var - evo.min_var, evo.npop, 1) + repmat(evo.min_var, evo.npop, 1);
    
    evo.q = zeros(evo.npop, evo.nvar);
    evo.f = zeros(evo.npop, evo.nobj);
    evo.g = zeros(evo.npop, evo.nobj);
    
    if(evo.nobj > 2)
        evo = Indicator_IDM_weight_init_MD(evo);
    else
        evo = Indicator_IDM_weight_init_2D(evo);
    end
    
    evo.H  = rand(evo.nvar, evo.nobj, evo.npop);
    evo.dx = zeros(evo.npop, evo.nvar);
    evo.dz = zeros(evo.npop, evo.nobj);
    evo.a  = ones(1, evo.npop);
end

%%
function evo = Indicator_IDM_weight_init_2D(evo)
    evo.w = zeros(evo.nweights, 2);
    evo.w(:, 1) = 0.0:(1.0/(evo.nweights-1)):1.0;
    evo.w(:, 2) = 1 - evo.w(:, 1);
end

%%
function evo = Indicator_IDM_weight_init_MD(evo)
    for i = 1:evo.nweights
        x = rand(1, evo.nobj-1);
        w = [x - [0, x(1:end-1)], 1-x(end)];
        evo.w(i, :) = w;
    end
    
    t = evo.nweights;
    t_max = 1e4 + evo.nweights;
    hv = zeros(evo.nweights + 1, evo.nobj);
    while(t < t_max)
        x = rand(1, evo.nobj-1);
        w = [x - [0, x(1:end-1)], 1-x(end)];
        evo.w = [evo.w; w];
        hv_all = lebesgue_measure(evo.w', ones(evo.nobj, 1)*2);
        for i = 1:size(evo.w, 1)
            W = evo.w; W(i, :) = [];
            hv(i) = hv_all - lebesgue_measure(W', ones(evo.nobj, 1)*2);
        end
        [~, id]= min(hv);
        evo.w(id(1), :) = [];
        t = t +1;
    end
end

%%
function evo = Indicator_IDM_dir(evo)
    for i = 1:evo.npop
        evo.dz(i, :) = -Indicator_IDM_Diff(evo, evo.f(i, :));
    end
end

%%
function evo = Indicator_IDM_Step(evo)
    global eval;
    
    f0 = Indicator_IDM_Function(evo, evo.f);
    for i = 1:evo.npop
        evo.dx(i, :) = evo.H(:, :, i)*evo.dz(i, :)';
        a = evo.a(i);
        if (a < 1e-6)
            a = 1.0;
        end

        a_ = 0;
        da = a;
        k  = 0;
        A = evo.f; A(i, :) = [];
        f1 = Indicator_IDM_Function(evo, evo.f(i, :)) - f0;
        while(da > 1e-6 && k < 3)
            evo.q(i, :) = evo.x(i, :) + (a_ + da)*evo.dx(i, :);
            evo.q(i, :) = max([evo.q(i, :); evo.min_var]);
            evo.q(i, :) = min([evo.q(i, :); evo.max_var]);
            evo.g(i, :) = evo.func(evo.q(i, :)); eval = eval+1;
            
            B = [A; evo.g(i, :)]; 
            f2 = Indicator_IDM_Function(evo, evo.g(i, :)) - Indicator_IDM_Function(evo, B);
            
            if(f1 > f2)
                a_ = a_ + da;
            else
                da = 0.5*da;
            end
            k = k + 1;
        end
        evo.a(i) = a_;
    end
end

%%
function evo = R2_IDM_Update(evo)
    global kappa;
    X = [evo.x; evo.q];
    F = [evo.f; evo.g];
    H = repmat(evo.H, 1, 1, 2);
    
    fitness =  zeros(1, size(evo.F, 1));
    evo.ideal = min(F) - max(max(F) - min(F));
    evo.nadir = max(F);
    
    evo = R2_IBEA_Update_Indicator(evo);
    for i = 1:size(F, 1)
        for j = 1:size(evo.R2, 2)
            fitness(i) = fitness(i) - exp(-evo.R2(j, i)/kappa);
        end
    end
    
    while(size(F, 1) > evo.npop)
        [~, id] = min(fitness);
        for i = 1:size(F, 1)
            fitness(i) = fitness(i) + exp(-evo.R2(id, i)/kappa);
        end
        X(id, :) = [];
        F(id, :) = [];
        H(:, :, d)= [];
        fitness = [evo.fitness(1:id-1), evo.fitness(id+1:end)];
        evo.R2(id, :) = [];
        evo.R2(:, id) = [];
    end
    evo.x = X;
    evo.f = F;
    evo.H = H; 
    
    evo = Tensor_IDM(evo);
    evo = Indicator_IDM_dir(evo);
end


%%
function f = Indicator_IDM_Function(evo, x)
    z  = evo.ideal;
    y = (repmat(z, size(x, 1), 1) - x);
    
    switch evo.funcType
        case 'R2'
            f = 0;
            t = zeros(1, size(x, 1));
            for i = 1:size(evo.w, 1)
                for j = 1:size(x, 1)
                    t(i) = max(evo.w(i, :).*y(j, :));
                end
                f = f + min(t);
            end
            f = f/size(evo.w, 1);
    end
end

%%
function df = Indicator_IDM_Diff(evo, x)
    z  = evo.ideal;
    y = (repmat(z, size(x, 1), 1) - x);
    
    switch evo.funcType
        case 'R2'
            df = zeros(1, evo.nobj);
            t  = zeros(1, size(x, 1));
            id = zeros(1, size(x, 1));
            for i = 1:size(evo.w, 1)
                for j = 1:size(x, 1)
                    [t(i), id(i)] = max(evo.w(i, :).*y(j, :));
                end
                [~, ix] = min(t);
                df(id(ix)) =  df(id(ix)) + 1;
            end
            df = df./size(evo.w, 1);
    end
end