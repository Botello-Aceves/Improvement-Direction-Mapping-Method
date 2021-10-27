function evo = IDM_Full_DO(evo, do_string)
    switch do_string
        case 'init'
            evo = All_IDM_init(evo);
        case 'step'
            evo = All_IDM_Step(evo);
        case 'update'
            evo = All_IDM_Update(evo);
        otherwise
            
    end
end

%%
function evo = All_IDM_init(evo)
    
    switch evo.dom
        case 'Oliver'
            evo = DS_Do(evo, 'init');
        case 'Pareto'
            evo = JacobianDirectionalBias_Do(evo, 'init');
        case 'Aggregation'
            evo = Aggregation_IDM_DO(evo, 'init');
        case 'Indicator'
            evo = R2_IDM_DO(evo, 'init');
        otherwise
            fprintf('Unknown metric');
            evo = NULL;
    end
end

%%
function evo = All_IDM_Step(evo)
    switch evo.dom
        case 'Oliver'
            evo = DS_Do(evo, 'step');
        case 'Pareto'
            evo = JacobianDirectionalBias_Do(evo, 'step');
        case 'Aggregation'
            evo = Aggregation_IDM_DO(evo, 'step');
        case 'Indicator'
            evo = R2_IDM_DO(evo, 'step');
        otherwise
            fprintf('Unknown metric');
            evo = NULL;
    end
end

%%
function evo = All_IDM_Update(evo)
    switch evo.dom
        case 'Oliver'
            evo = DS_Do(evo, 'update');
        case 'Pareto'
            evo = JacobianDirectionalBias_Do(evo, 'update');
        case 'Aggregation'
            evo = Aggregation_IDM_DO(evo, 'update');
        case 'Indicator'
            evo = R2_IDM_DO(evo, 'update');
        otherwise
            fprintf('Unknown metric');
            evo = NULL;
    end
end