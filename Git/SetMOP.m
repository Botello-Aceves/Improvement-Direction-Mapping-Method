function [PF, mop] = SetMOP(name)
    global M k l;
    
    switch name
        case 'sphere'
            mop = benchmark('sphere', 2);
            paretoName = 'PF/Schaffer.pf';
            PF = 2*dlmread(paretoName);
            mop.ref = max(PF)';
        case 'schaffer'
            mop = benchmark('schaffer', 0);
            paretoName = 'PF/Schaffer.pf';
            PF = dlmread(paretoName);
            mop.ref = max(PF)';
        case 'schaffer2'
            mop = benchmark('schaffer2', 0);
            paretoName = 'PF/Schaffer2.pf';
            PF = dlmread(paretoName);
            mop.ref = max(PF)';
        case 'fonseca'
            mop = benchmark('fonseca', 10);
            paretoName = 'PF/Fonseca.pf';
            PF = dlmread(paretoName);
            mop.ref = max(PF)';
        case 'kursawe'
            mop = benchmark('kursawe', 0);
            paretoName = 'PF/Kursawe.pf';
            PF = dlmread(paretoName);
            mop.ref = max(PF)';
        case 'poloni'
            mop = benchmark('poloni', 0);
            paretoName = 'PF/Poloni.pf';
            PF = dlmread(paretoName);
            mop.ref = max(PF)';
        case 'viennet'
            mop = benchmark('viennet', 0);
            paretoName = 'PF/Viennet.pf';
            PF = dlmread(paretoName);
            mop.ref = max(PF)';
        case 'viennet2'
            mop = benchmark('viennet2', 0);
            paretoName = 'PF/Viennet2.pf';
            PF = dlmread(paretoName);
            mop.ref = max(PF)';
        case 'viennet3'
            mop = benchmark('viennet3', 0);
            paretoName = 'PF/Viennet3.pf';
            PF = dlmread(paretoName);
            mop.ref = max(PF)';
        case 'ZDT1'
            mop = benchmark('zdt1', 30);
            paretoName = ['PF/ZDT/', mop.name,'.pf'];
            PF = dlmread(paretoName);
            mop.ref = [11; 11];
        case 'ZDT2'
            mop = benchmark('zdt2', 30);
            paretoName = ['PF/ZDT/', mop.name,'.pf'];
            PF = dlmread(paretoName);
            mop.ref = [11; 11];
        case 'ZDT3'
            mop = benchmark('zdt3', 30);
            paretoName = ['PF/ZDT/', mop.name,'.pf'];
            PF = dlmread(paretoName);
            mop.ref = [11; 11];
        case 'ZDT4'
            mop = benchmark('zdt4', 10);
            paretoName = ['PF/ZDT/', mop.name,'.pf'];
            PF = dlmread(paretoName);
            mop.ref = [11; 11];
        case 'ZDT6'
            mop = benchmark('zdt6', 10);
            paretoName = ['PF/ZDT/', mop.name,'.pf'];
            PF = dlmread(paretoName);
            mop.ref = [11; 11];
        case 'DTLZ1'
            M = 2; k = 5; l = 1;
            mop = benchmark('DTLZ1', M + k - 1);
            paretoName = ['PF/DTLZ.', num2str(mop.od), 'D/', mop.name, '.', num2str(mop.od), 'D.pf'];
            PF = dlmread(paretoName);
            mop.ref = [11; 11];
        case 'DTLZ2'
            M = 2; k = 10; l = 1;
            mop = benchmark('DTLZ2', M + k - 1);
            paretoName = ['PF/DTLZ.', num2str(mop.od), 'D/', mop.name, '.', num2str(mop.od), 'D.pf'];
            PF = dlmread(paretoName);
            mop.ref = [11; 11];
        case 'DTLZ3'
            M = 2; k = 10; l = 1;
            mop = benchmark('DTLZ3', M + k - 1);
            paretoName = ['PF/DTLZ.', num2str(mop.od), 'D/', mop.name, '.', num2str(mop.od), 'D.pf'];
            PF = dlmread(paretoName);
            mop.ref = [11; 11];
        case 'DTLZ4'
            M = 2; k = 10; l = 1;
            mop = benchmark('DTLZ4', M + k - 1);
            paretoName = ['PF/DTLZ.', num2str(mop.od), 'D/', mop.name, '.', num2str(mop.od), 'D.pf'];
            PF = dlmread(paretoName);
            mop.ref = [11; 11];
        case 'DTLZ5'
            M = 2; k = 10; l = 1;
            mop = benchmark('DTLZ5', M + k - 1);
            paretoName = ['PF/DTLZ.', num2str(mop.od), 'D/', mop.name, '.', num2str(mop.od), 'D.pf'];
            PF = dlmread(paretoName);
            mop.ref = [11; 11];
        case 'DTLZ6'
            M = 2; k = 10; l = 1;
            mop = benchmark('DTLZ6', M + k - 1);
            paretoName = ['PF/DTLZ.', num2str(mop.od), 'D/', mop.name, '.', num2str(mop.od), 'D.pf'];
            PF = dlmread(paretoName);
            mop.ref = [11; 11];
        case 'DTLZ7'
            M = 2; k = 10; l = 1;
            mop = benchmark('DTLZ7', M + k - 1);
            paretoName = ['PF/DTLZ.', num2str(mop.od), 'D/', mop.name, '.', num2str(mop.od), 'D.pf'];
            PF = dlmread(paretoName);
            mop.ref = [11; 11];
        case 'WFG1'
            M = 2; k = 20; l = 4;
            mop = benchmark('wfg1', k + l);
            paretoName = ['PF/WFG.', num2str(mop.od), 'D/', mop.name, '.', num2str(mop.od), 'D.pf'];
            PF = dlmread(paretoName);
            mop.ref = [11; 11];
        case 'WFG2'
            M = 2; k = 20; l = 4;
            mop = benchmark('wfg2', k + l);
            paretoName = ['PF/WFG.', num2str(mop.od), 'D/', mop.name, '.', num2str(mop.od), 'D.pf'];
            PF = dlmread(paretoName);
            mop.ref = [11; 11];
        case 'WFG3'
            M = 2; k = 20; l = 4;
            mop = benchmark('wfg3', k + l);
            paretoName = ['PF/WFG.', num2str(mop.od), 'D/', mop.name, '.', num2str(mop.od), 'D.pf'];
            PF = dlmread(paretoName);
            mop.ref = [11; 11];
        case 'WFG4'
            M = 2; k = 20; l = 4;
            mop = benchmark('wfg4', k + l);
            paretoName = ['PF/WFG.', num2str(mop.od), 'D/', mop.name, '.', num2str(mop.od), 'D.pf'];
            PF = dlmread(paretoName);
            mop.ref = [11; 11];
        case 'WFG5'
            M = 2; k = 20; l = 4;
            mop = benchmark('wfg5', k + l);
            paretoName = ['PF/WFG.', num2str(mop.od), 'D/', mop.name, '.', num2str(mop.od), 'D.pf'];
            PF = dlmread(paretoName);
            mop.ref = [11; 11];
        case 'WFG6'
            M = 2; k = 20; l = 4;
            mop = benchmark('wfg6', k + l);
            paretoName = ['PF/WFG.', num2str(mop.od), 'D/', mop.name, '.', num2str(mop.od), 'D.pf'];
            PF = dlmread(paretoName);
            mop.ref = [11; 11];
        case 'WFG7'
            M = 2; k = 20; l = 4;
            mop = benchmark('wfg7', k + l);
            paretoName = ['PF/WFG.', num2str(mop.od), 'D/', mop.name, '.', num2str(mop.od), 'D.pf'];
            PF = dlmread(paretoName);
            mop.ref = [11; 11];
        case 'WFG8'
            M = 2; k = 20; l = 4;
            mop = benchmark('wfg8', k + l);
            paretoName = ['PF/WFG.', num2str(mop.od), 'D/', mop.name, '.', num2str(mop.od), 'D.pf'];
            PF = dlmread(paretoName);
            mop.ref = [11; 11];
        case 'WFG9'
            M = 2; k = 20; l = 4;
            mop = benchmark('wfg9', k + l);
            paretoName = ['PF/WFG.', num2str(mop.od), 'D/', mop.name, '.', num2str(mop.od), 'D.pf'];
            PF = dlmread(paretoName);
            mop.ref = [11; 11];
        case 'UF1'
            mop = benchmark('uf1', 30);
            paretoName = 'PF/UF/UF1.pf';
            PF = dlmread(paretoName);
            mop.ref = [11; 11];
%             mop.ref = max(PF)';
        case 'UF2'
            mop = benchmark('uf2', 30);
            paretoName = 'PF/UF/UF2.pf';
            PF = dlmread(paretoName);
            mop.ref = [11; 11];
%             mop.ref = max(PF)';
        case 'UF3'
            mop = benchmark('uf3', 30);
            paretoName = 'PF/UF/UF3.pf';
            PF = dlmread(paretoName);
            mop.ref = [11; 11];
%             mop.ref = max(PF)';
        case 'UF4'
            mop = benchmark('uf4', 30);
            paretoName = 'PF/UF/UF4.pf';
            PF = dlmread(paretoName);
            mop.ref = [11; 11];
%             mop.ref = max(PF)';
        case 'UF5'
            mop = benchmark('uf5', 30);
            paretoName = 'PF/UF/UF5.pf';
            PF = dlmread(paretoName);
            mop.ref = [11; 11];
%             mop.ref = max(PF)';
        case 'UF6'
            mop = benchmark('uf6', 30);
            paretoName = 'PF/UF/UF6.pf';
            PF = dlmread(paretoName);
            mop.ref = [11; 11];
%             mop.ref = max(PF)';
        case 'UF7'
            mop = benchmark('uf7', 30);
            paretoName = 'PF/UF/UF7.pf';
            PF = dlmread(paretoName);
            mop.ref = [11; 11];
%             mop.ref = max(PF)';
        case 'UF8'
            mop = benchmark('uf8', 30);
            paretoName = 'PF/UF/UF8.pf';
            PF = dlmread(paretoName);
            mop.ref = [11; 11; 11];
%             mop.ref = max(PF)';
        case 'UF9'
            mop = benchmark('uf9', 30);
            paretoName = 'PF/UF/UF9.pf';
            PF = dlmread(paretoName);
            mop.ref = [11; 11; 11];
%             mop.ref = max(PF)';
        case 'UF10'
            mop = benchmark('uf10', 30);
            paretoName = 'PF/UF/UF10.pf';
            PF = dlmread(paretoName);
            mop.ref = [11; 11; 11];
%             mop.ref = max(PF)';    
        case 'zdt3_2D'
            mop = benchmark('zdt3', 2);
            paretoName = ['PF/ZDT/', mop.name,'.pf'];
            PF = dlmread(paretoName);
            mop.ref = [11; 11];
        case 'DTLZ1_3D'
            M = 3; k = 5; l = 1;
            mop = benchmark('DTLZ1', M + k - 1);
            paretoName = ['PF/DTLZ.', num2str(mop.od), 'D/', mop.name, '.', num2str(mop.od), 'D.pf'];
            PF = dlmread(paretoName);
            mop.ref = [11; 11];
        case 'DTLZ2_3D'
            M = 3; k = 10; l = 1;
            mop = benchmark('DTLZ2', M + k - 1);
            paretoName = ['PF/DTLZ.', num2str(mop.od), 'D/', mop.name, '.', num2str(mop.od), 'D.pf'];
            PF = dlmread(paretoName);
            mop.ref = [11; 11];
        case 'DTLZ3_3D'
            M = 3; k = 10; l = 1;
            mop = benchmark('DTLZ3', M + k - 1);
            paretoName = ['PF/DTLZ.', num2str(mop.od), 'D/', mop.name, '.', num2str(mop.od), 'D.pf'];
            PF = dlmread(paretoName);
        case 'DTLZ4_3D'
            M = 3; k = 10; l = 1;
            mop = benchmark('DTLZ4', M + k - 1);
            paretoName = ['PF/DTLZ.', num2str(mop.od), 'D/', mop.name, '.', num2str(mop.od), 'D.pf'];
            PF = dlmread(paretoName);
            mop.ref = [11; 11];
        case 'fonseca2D'
            mop = benchmark('fonseca', 2);
            paretoName = 'PF/Fonseca.pf';
            PF = dlmread(paretoName);
            mop.ref = max(PF)';
        otherwise
            fprintf('Uknown MOP');
    end 
end