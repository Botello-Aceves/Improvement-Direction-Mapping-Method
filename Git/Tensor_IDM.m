function evo = Tensor_IDM(evo)
    global N;
    switch evo.tensor(1)
        case 'Broyden'
            for i = 1:evo.npop
%                 H = rand(evo.nvar, evo.nobj);
                H = evo.H(: , :, i);
                p = evo.dz(i, :);
                x  = evo.x(i, :);
                f  = evo.f(i, :);
                
                d = p*H';
                x_ = x + d;
                x_ = max(x_, evo.min_var(:, 1));
                x_ = min(x_, evo.max_var(:, 2));
                f_ = evo.func(x_);

                kk = 1;
                y = f_ - f;
                eta = y/norm(y)*(p/norm(p))';
                while(eta < 1-1e-6 && kk < 1e2 && norm(x - x_) > 1e-6)
                    H_ = H'*H;
                    H = H*(eye(length(f)) + (p'-y')/(p*H_*y')*p*H_);
                    H = H/norm(H, 'fro');
                    if(max(max(isnan(H))) || max(max(isinf(H))))
                        H = evo.H(:, :, i) + (0.2*rand(length(x), length(f)) - 0.1);
                        H = H/norm(H, 'fro');
                    end
                    
                    d = p*H';
                    x_ = x + d;
                    x_ = max(x_, evo.min_var(:, 1));
                    x_ = min(x_, evo.max_var(:, 2));
                    f_ = evo.func(x_);
                                       
                    y = f_ - f;
                    eta = y/norm(y)*(p/norm(p))';
                    kk = kk + 1;
                end
                evo.H(:, :, i) = H/norm(H, 'fro');
            end
            
        case 'Analityc'
            for i = 1:evo.npop
                n = evo.nvar;
                sz = min(evo.npop, n*(n+1)+1);
                s = evo.x(i, :) + 1e-9*(2*rand(sz, n) - 1);
                for j = 1:sz
                    s(j, :) = max(s(j, :), evo.min_var(:, 1));
                    s(j, :) = min(s(j, :), evo.max_var(:, 2));
                    z(j, :) = evo.func(s(j, :));
                end

                A = MultivarQuadraticLeastSquareRegression(s, z);
                [~, H] = QuadraticFunction(A, evo.x(i, :));
                evo.H(:, :, i) = Tensor_trans(evo, H);
            end
        case 'Finite'
            for l = 1:evo.npop
                A = evo.f;
                xi = 1:size(A, 1);
                
                [a, ix] = ismember(evo.f(l, :), A, 'rows');
                while(a == true)
                    A(ix, :) = [];
                    xi(ix) = [];
                    [a, ix] = ismember(evo.f(l, :), A, 'rows');
                end
                ss = min(length(xi), N);
                if(ss > 0)
                   ix = setNeighborhood(evo.x(l, :), evo.x(xi, :), ss);
                   dz = evo.f(xi(ix), :) - repmat(evo.f(l, :), ss, 1);
                   dx = evo.x(xi(ix), :) - repmat(evo.x(l, :), ss, 1);

                   H = zeros(evo.nobj, evo.nvar);
                    for k = 1:ss
                         for i = 1:evo.nobj
                            for j = 1:evo.nvar
                                H(i, j) = H(i, j) + ...
                                    dz(k, i)*dx(k,j)./norm(dx(k, :)); 
                            end
                         end
                    end
                    H = (1/min(length(xi), N))*H;
                    evo.H(:, :, l) = Tensor_trans(evo, H);
                 else
                    evo.H(:, :, l) = eye(evo.nvar, evo.nobj);
                 end
            end
            
        case 'Regression'
            for i = 1:evo.npop
                ix = setNeighborhood(evo.x(i, :), evo.x, min(evo.nvar + 1, N));
                A = MultivarLinearRegression(evo.x(ix, :), evo.f(ix, :));
                [~, H] = LinearFunction(A, evo.x(i, :));
                evo.H(:, :, i) = Tensor_trans(evo, H);
            end
        otherwise
            evo.H = NaN;
    end
end

%%
function H = Tensor_trans(evo, H)
    switch evo.tensor(2)
        case 'Inverse'
            H = pinv(H);
        case 'Transpose'
            H = H';
        case 'Normal'
            H = Normalized_PseudoInverse(H);
        otherwise
            H = NaN;
    end
end