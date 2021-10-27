function evo = Aggregation_IDM_DO(evo, do_string)
    switch do_string
        case 'init'
            evo = Aggregation_IDM_init(evo);
        case 'step'
            evo = Aggregation_IDM_Step(evo);
        case 'update'
            evo = Aggregation_IDM_Update(evo);
        otherwise
            
    end
end

%%
function evo = Aggregation_IDM_init(evo)
    evo.x = rand(evo.npop, evo.nvar);
    evo.x = evo.x.*repmat(evo.max_var - evo.min_var, evo.npop, 1) + repmat(evo.min_var, evo.npop, 1);
    
    evo.q = zeros(evo.npop, evo.nvar);
    evo.f = zeros(evo.npop, evo.nobj);
    evo.g = zeros(evo.npop, evo.nobj);
    
    if(evo.nobj > 2)
        evo = Aggregation_IDM_weight_init_MD(evo);
    else
        evo = Aggregation_IDM_weight_init_2D(evo);
    end
    
    evo.ideal  = realmax*ones(1, evo.nobj);
    evo.zenith = realmax*ones(1, evo.nobj);
    evo.nadir  = realmin*ones(1, evo.nobj);
    
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
function evo = Aggregation_IDM_Step(evo)
    global eval;
    for i = 1:evo.npop
        if (evo.a(i) < 1e-6)
            evo.a(i) = 1;
        end
        
        evo.dz(i, :) = -Aggregation_IDM_Diff(evo, evo.f(i, :), evo.w(i, :));
%         evo.dz(i, :) = evo.dz(i, :)/norm(evo.dz(i, :));
        evo.dx(i, :) = evo.H(:, :, i)*evo.dz(i, :)';
        
        a_ = 0;
        da = evo.a(i);
        evo.q(i, :) = evo.x(i, :) + (a_ + da)*evo.dx(i, :);
        evo.q(i, :) = max([evo.q(i, :); evo.min_var]);
        evo.q(i, :) = min([evo.q(i, :); evo.max_var]);    
        evo.g(i, :) = evo.func(evo.q(i, :)); eval = eval+1;
        
        k  = 1;
        f1 = Aggregation_IDM_Function(evo, evo.f(i, :), evo.w(i, :));
        f2 = Aggregation_IDM_Function(evo, evo.g(i, :), evo.w(i, :));
        while(da > 1e-6 && k < 3)
            evo.q(i, :) = evo.x(i, :) + (a_ + da)*evo.dx(i, :);
            evo.q(i, :) = max([evo.q(i, :); evo.min_var]);
            evo.q(i, :) = min([evo.q(i, :); evo.max_var]);
            evo.g(i, :) = evo.func(evo.q(i, :)); eval = eval+1;
            f2 = Aggregation_IDM_Function(evo, evo.g(i, :), evo.w(i, :));
            
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
function evo = Aggregation_IDM_Update(evo)
    X = [evo.x; evo.q];
    F = [evo.f; evo.g];
    H = repmat(evo.H, 1, 1, 2);

%     evo.ideal = min(F) - max(max(F) - min(F));
    evo.zenith = min([evo.zenith; F]);
    evo.nadir = max([evo.nadir; F]);
    evo.ideal = evo.zenith;% - 0.01*(evo.nadir - evo.zenith);
    
    W = zeros(evo.nweights, size(F, 1));
    for i = 1:evo.nweights
        for j = 1:size(F, 1)
           W(i, j) = Aggregation_IDM_Function(evo, F(j, :), evo.w(i, :));
        end
    end
    
    ix = []; iw = [];
    xi = 1:size(F);
    wi = 1:evo.nweights;
    while(length(ix) < evo.npop)
%         [z, ID] = min(W, [], 2);
%         [~, id] = min(z);
        [nW, nI] = find(W == min(min(W)));
        
        W(nW, :) = [];
        W(:, nI) = [];
        
        ix = [ix, xi(nI)];
        iw = [iw, wi(nW)];
        xi(nI) = [];
        wi(nW) = [];
    end

    evo.x(iw, :) = X(ix, :);
    evo.f(iw, :) = F(ix, :);
    evo.H(:, :, iw) = H(:, :, ix);
    evo = Tensor_IDM(evo);
end

%%
function f = Aggregation_IDM_Function(evo, x, w)
    m  = evo.nobj;
    z  = evo.ideal;
    mx = evo.nadir;

    switch evo.funcType
        case 'GSS'
            a = 1.0;
            b = 0.1;
            w_ = (1./w)/norm(1./w);
            SIGMA = w_'*w_ + b*eye(m);
            y = (x - z)./mx;
            f = 1.0 - exp(-a*y*((SIGMA)\y'));

        case 'PSS'
            a = 1.0;
            b = 0.1;
            w_ = (1./w)/norm(1./w);
            SIGMA = w_'*w_ + b*eye(m);
            y = (x - z)./mx;
            f = a*y*((SIGMA)\y');

        case 'WCP'
            p = 9.0;
            y = (x - z)./mx;
            f = sum((w.*y).^p);

        case 'WS'
            y = (x - z)./mx;
            f = sum(w.*y);

        case 'EWC'
            p = 1.0;
            y = (x - z)./mx;
            f = sum((exp(p*w) - 1).*exp(p*y));
        
        case 'WPO'
            p = 3.0;
            y = (x - z)./mx;
            f = sum(w.*y.^p);

        case 'WPR'
            y = (x - z)./mx;
            f = prod(y.^w);

        case 'WN'
            p = 10;
            y = (x - z)./mx;
            f = sum(w.*y.^p)^(1/p);

        case 'TCH'
            y = (x - z)./mx;
            f = max(w.*y);

        case 'ATCH'
            a = 0.1;
            y = (x - z)./mx;
            f = max(w.*y) + a*sum(y);

        case 'MTCH'
            a = 0.1;
            y = (x - z)./mx;
            f = max(w.*y + a*sum(y));

        case 'ASF'
            y = (x - z)./mx;
            f = max(y./w);

        case 'AASF'
            a = 0.01;
            y = (x - z)./mx;
            f = max(y./w) + a*sum(y./w);

        case 'PBI'
            t = 5.0;
            y = (x - z)./mx;
            d1 = abs(dot(y, w/norm(w)));
            d2 = norm(y-d1*w/norm(w));
            f = d1 + t*d2;

        case 'IPBI'
            t = 5.0;
            y = (z - x)./mx;
            d1 = abs(dot(y, w/norm(w)));
            d2 = norm(y-d1*w/norm(w));
            f = t*d2 - d1;

        case 'CS'
            a = min(w);
            y = (x - z)./mx;
            f = sum(w.*y) + a*sum(abs(y));

        case 'VADS'
            p = 2;
            y = (x - z)./mx;
            u = norm(y);
            l = (w/norm(w))*(y/norm(y))';
            f = u/(l)^p;
    end
end

%%
function df = Aggregation_IDM_Diff(evo, x, w)
    m  = evo.nobj;
    z  = evo.ideal;
    mx = evo.nadir;

    switch evo.funcType
        case 'GSS'
            a = 1.0;
            b = 0.1;
            w_ = (1./w)/norm(1./w);
            SIGMA = w_'*w_ + b*eye(m);
            y = (x - z)./mx;
            df = 2*(a*(SIGMA\y'))*exp(-a*y*(SIGMA\y'));

        case 'PSS'
            a = 1.0;
            b = 0.1;
            w_ = (1./w)/norm(1./w);
            SIGMA = w_'*w_ + b*eye(m);
            y = (x - z)./mx;
            df = 2*a*(SIGMA\y');

        case 'WCP'
            p = 9.0;
            y = (x - z)./mx;
            df = (w.^p).*(y.^(p-1));

        case 'WS'
            df = w;

        case 'EWC'
            p = 1.0;
            y = (x - z)./mx;
            df = p*(exp(p*w )- 1).*exp(p*y);
        
        case 'WPO'
            p = 3.0;
            y = (x - z)./mx;
            df = p*w.*y.^(p-1);

        case 'WPR'
            y = (x - z)./mx;
            df = y.^w;
            for j = 1:m
                df(j) = w(j)*df(j)/y(j);
            end

        case 'WN'
            p = 10;
            y = (x - z)./mx;
            df = w.*y.^(p-1)*sum(w.*y.^p)^((1-p)/p);

        case 'TCH'
            y = x; %(x - z)./mx;
            [~, id] = max(w.*y);
            df = zeros(1, m);
            df(id) = w(id);

        case 'ATCH'
            a = 0.1;
            y = x; %(x - z)./mx;
            [~, id] = max(w.*y);
            df = a*ones(1, m);
            df(id) = w(id)+a;

        case 'MTCH'
            a = 0.1;
            y = x; %(x - z)./mx;
            df = zeros(1, m);
            [~, id] = max(w.*y + a*sum(y));
            df(id) = w(id)+a;

        case 'ASF'
            y = x; %(x - z)./mx;
            df = zeros(1, m);
            [~, id] = max(w./y);
            df(id) = 1./w(id);

        case 'AASF'
            a  = 0.01;
            y = x; %(x - z)./mx;
            df = zeros(1, m);
            [~, id] = max(w./y);
            df(id) = 1./w(id) + a;

        case 'PBI'
            t  = 5.0;
            y = (x - z)./mx;
            h  = w/norm(w);
            d1 = abs(dot(y, h));
            d2 = norm(y-d1*h);
            df = h + (t/d2)*(y - h*d1).*(1.0 - h.^2);

        case 'IPBI'
            t  = 5.0;
            y = (z - x)./mx;
            h  = w/norm(w);
            d1 = abs(dot(y, h));
            d2 = norm(y-d1*h);
            df = h - (t/d2)*(y - h*d1).*(1.0 - h.^2);

        case 'CS'
            a  = min(w);
            df = w + a;

        case 'VADS'
            p = 2;
            y = (x - z)./mx;
            u = norm(y); 
            v = norm(w);
            l = dot(w/v, y/u);
            du = y/u; 
            dl = p*l^(p-1)*(w/v).*((u^2 - y)/u^3);
            df = (du - dl)/l^2;
    end
end