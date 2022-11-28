function plot_sobol_indices_dt(class_scenario,undg_vax,fac_grad_vax,...
                                      param_set,print_contact_info,...
                                      print_ylabel,SubSampN,font_size)
% function: (dorm_living_scenario,undg_vax,fac_grad_vax,fig_width,fig_height)
% class cap  
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
% fig_width  = fig_dim(1);
% fig_height = fig_dim(2);

% INFECTION PARAMETERS
% Beta:(9) transmission rate
% sig:(11) 1/E(time in exposed state)
% phi:(12) Pr(becoming asymptomatic)
% alp:(13) Pr(quarantine|symptomatic)

% SOCIAL CONTACT PARAMETERS
% M:(2) Fraction infected people in Merced
% m:(10) mask usage
% c:(1) Community contact multiplier
% p:(3) party size
% w:(4) Weekend contact rate multiplier

% INITIAL CONDITIONS
% I^s_u:(5)  symptomatic live-off campus undergraduate
% I^a_u:(6)  asymptomatic live off-campus undergraduate
% I^s_d:(7)  symptomatic live on-campus undergraduate
% I^a_d:(8)  asymptomatic live on-campus undergraduate
% I^s_g:(9)  symptomatic graduate students
% I^a_g:(10) asymptomatic graduate students
% I^s_f:(11) symptomatic faculty 
% I^a_f:(12) asymptomatic faculty

% Order of 19 parameters
SI_order = [5 7:8 18 19    2 1 4 3 6 9    10:17];% ALL 19 parameters
IP = 1:5;   % Infection parameters
SC = 6:11;  % Social contact parameters
IC = 12:19; % Initial conditions 

% Load parameters
[prms_info,paper_prms_str] = fun_model_parameter_ranges;
paper_prms_str = paper_prms_str(prms_info(:,1)==1);
paper_prms_str = paper_prms_str(SI_order);
% SubSampN  =  1500*(sum(prms_info(:,1))+2);

C = linspecer(2);
y_sol_idx = 1; % Cumulative Infections     

% Load baseline case to determine expected cumulative infections
load('./sobol_indices_dt/Sobol_dt_sol1sub31500rep2000CC_50000Vu0Vd0Vg0Vf0.mat','Stats');

Baseline_CI = mean(Stats.mu.CI); % Save mean cumulative infectiosn

vG = fac_grad_vax; vF = fac_grad_vax;
vU = undg_vax; vD = undg_vax;

load(['sobol_indices_dt/Sobol_dt_sol',...
    num2str(y_sol_idx),...
        'sub',num2str(SubSampN),...
            'rep2000','CC_',num2str(class_scenario),...
            'Vu',num2str(vU),'Vd',num2str(vD),'Vg',num2str(vG),'Vf',num2str(vF),'.mat'])

% SOBOL: infection and contact parameters
if contains(param_set,'inf') && contains(param_set,'con')
    FOI.error   = FirstOrderIdx.error(SI_order); 
    FOI.error   = FOI.error([IP SC]);
    FOI.indices = FirstOrderIdx.indices(SI_order);
    FOI.indices = FOI.indices([IP SC]);
    TOI.error   = TotalOrderIdx.error(SI_order); 
    TOI.error   = TOI.error([IP SC]);
    TOI.indices = TotalOrderIdx.indices(SI_order);
    TOI.indices = TOI.indices([IP SC]);
    leg_bool = false;
    axlab_bool=0;
    barplot_sobol_indices(FOI,TOI,paper_prms_str([IP SC]),leg_bool,axlab_bool);
    if print_contact_info
        text(1,0.78,['$\Delta T$ = ',sprintf('%0.1f',mean(Stats.mu.DT)/7),...
            '(',sprintf('%0.2f',sqrt(mean(Stats.var.DT)/7)/mean(Stats.mu.DT/7)),') weeks'],...
            'Interpreter','latex','FontSize',font_size,'HorizontalAlignment','left',...
            'VerticalAlignment','middle');
        text(1,0.90,['$E\left[C(t_\textit{end})\right]$ = ',sprintf('%0.1f',mean(Stats.mu.CI)),...
                            '(',sprintf('%0.2f',sqrt(mean(Stats.var.CI))/mean(Stats.mu.CI)),')'],...
            'Interpreter','latex','FontSize',font_size,'HorizontalAlignment','left',...
            'VerticalAlignment','middle');
        if vF == 0 && vU == 0 && class_scenario == 50000
            RECB = 0;
        else
            RECB = 1-mean(Stats.mu.CI)/Baseline_CI;
        end
        text(1,0.68,['RECI $= ',sprintf('%0.1f',RECB*100),'\%$'],...
                    'HorizontalAlignment','left',...
                        'VerticalAlignment','middle',...
                            'Interpreter','latex','FontSize',font_size)
    end
    set(gca,'TickLabelInterpreter','latex','FontSize',font_size+1)
end

% SOBOL: initial conditions
if contains(param_set,'init')
    FOI.error   = FirstOrderIdx.error(SI_order); 
    FOI.error   = FOI.error(IC);
    FOI.indices = FirstOrderIdx.indices(SI_order);
    FOI.indices = FOI.indices(IC);
    TOI.error   = TotalOrderIdx.error(SI_order); 
    TOI.error   = TOI.error(IC);
    TOI.indices = TotalOrderIdx.indices(SI_order);
    TOI.indices = TOI.indices(IC);
    leg_bool = false;
    axlab_bool=0;
    barplot_sobol_indices(FOI,TOI,paper_prms_str([IC]),leg_bool,axlab_bool);
    if print_contact_info
        text(1,0.78,['$\Delta T$ = ',sprintf('%0.1f',mean(Stats.mu.DT)/7),...
            '(',sprintf('%0.2f',sqrt(mean(Stats.var.DT)/7)/mean(Stats.mu.DT/7)),') weeks'],...
            'Interpreter','latex','FontSize',font_size,'HorizontalAlignment','left');
        text(1,0.90,['$E\left[C(t_\textit{end})\right]$ = ',sprintf('%0.1f',mean(Stats.mu.CI)),...
                            '(',sprintf('%0.2f',sqrt(mean(Stats.var.CI))/mean(Stats.mu.CI)),')'],...
            'Interpreter','latex','FontSize',font_size,'HorizontalAlignment','left');
        if vF == 0 && vU == 0 && class_scenario == 50000
            RECB = 0;
        else
            RECB = 1-mean(Stats.mu.CI)/Baseline_CI;
        end
        text(1,0.68,['RECI $= ',sprintf('%0.1f',RECB*100),'\%$'],...
                    'HorizontalAlignment','left',...
                        'VerticalAlignment','middle',...
                            'Interpreter','latex','FontSize',font_size)
    end
    set(gca,'TickLabelInterpreter','latex','FontSize',font_size+1)
end

% Add descriptive y-label to figure
if print_ylabel
    ylabel(['\boldmath$v_{f,0} =  v_{g,0} = ',num2str(fac_grad_vax),'\%$',...
        char(10),'\boldmath$v_{u,0} = v_{d,0} = ',num2str(undg_vax),'\%$',...
            char(10),'Fraction of Variance'],...
                'Interpreter','latex','FontSize',font_size+2)
end 