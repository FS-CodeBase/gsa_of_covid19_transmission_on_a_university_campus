function plot_sobol_indices_time(class_scenario,undg_vax,fac_grad_vax,...
                                	param_set,print_contact_info,...
                                    print_ylabel,Sobol_order,SubSampN,...
                                    font_size)
% Script: plot_sobol_indices_in_time
% Author: Fabian Santiago
% Date: 8/11/21
% Description: This script plots the Sobol Indices for the ``cumulative
%              number'' of infections during the course of the semester.
%              First 30 days are refined then, every 3.5 days are shown.
% Inputs: 
% clear
% close all
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
load('./colorblind_colormap/colorblind_colormap.mat','colorblind')
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

% Order of 19 parameters from param_vax_model_parameter_ranges_cmt.
SI_order = [5 7:8 18 19    2 1 4 3 6 9    10:17];% ALL 19 parameters
IP = 1:5;   % Infection parameters
SC = 6:11;   % Social contact parameters
IC = 12:19; % Initial conditions 

[prms_info,paper_prms_str] = fun_model_parameter_ranges;
% [prms_info,paper_prms_str] = paper_vax_model_parameter_ranges_cmt;
paper_prms_str = paper_prms_str(prms_info(:,1)==1);
paper_prms_str = paper_prms_str(SI_order);

% Model solution length in days
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DaysRefined = 30; % Number of days model is solved refined
MaxDays = 7*15; % End of semester (15 week semester)

% Refine solutions during the first 30 days, then every 3 days
tspan = 24*[0:DaysRefined (DaysRefined+3):3:MaxDays];
tidx = [1:3:DaysRefined 31:56];
% C = linspecer(2);
% SubSampN  = 1500*(sum(prms_info(:,1))+2);
y_sol_idx = 1; % Cumulative Infections
LW = linspace(1,1,100); % For different line widths (not very useful)

% Labels for weeks since semester started
xlabs = repmat({''},1,16);
xlabs(1:3:16) = compose('%0.0f',0:3:15); % Weeks

% Load baseline case to determine expected cumulative infections
% load('./sobol_indices_dt/Sobol_dt_sol1sub6rep2CC_50000Vu0Vd0Vg0Vf0.mat','Stats');
load('./sobol_indices_dt/Sobol_dt_sol1sub31500rep2000CC_50000Vu0Vd0Vg0Vf0.mat','Stats');
Baseline_CI = mean(Stats.mu.CI); % Save mean cumulative infectiosn

% for dorm_living_scenario = CON_MATS(IDX,:)
vG = fac_grad_vax; vF = fac_grad_vax;
vU = undg_vax; vD = undg_vax;

load(['sobol_indices_in_time/SobolT_'...
    'sol',num2str(y_sol_idx),...
        'sub',num2str(SubSampN),...
            'rep2000','CC_',num2str(class_scenario),...
            'Vu',num2str(vU),'Vd',num2str(vD),'Vg',num2str(vG),'Vf',num2str(vF),'.mat'])

if contains(lower(Sobol_order),'fo')
    SOBOL.error   = FirstOrderIdx.error(:,SI_order); 
    SOBOL.indices = FirstOrderIdx.indices(:,SI_order);
elseif contains(lower(Sobol_order),'te')
    SOBOL.error   = TotalOrderIdx.error(:,SI_order); 
    SOBOL.indices = TotalOrderIdx.indices(:,SI_order);
else
    error(['Sobol order must be ''FO'' for first order',...
            'or ''TE'' for total effect']);
end


class_cap_str = num2str(class_cap);
if 50000 == class_cap
    class_cap_str = 'None';
end

% figure
%%%%%%%%%%%%%%%%%%%%%%% PLOT INFECTION PARAMETERS %%%%%%%%%%%%%%%%%%%%%
if contains(param_set,'inf') % Check if these should be plotted
	hIP = zeros(1,6);
    for i = IP
        hIP(i) = errorbar(tspan(tidx)/(24*7),SOBOL.indices(tidx,i),SOBOL.error(tidx,i),'LineWidth',LW(10-i),...
            'Color',colorblind(13-i,:),'MarkerSize',5); hold on
    end
    % CASE DETAILS in TEXT ON PLOT
    load(['sobol_indices_dt/Sobol_dt_sol',...
    num2str(y_sol_idx),...
        'sub',num2str(SubSampN),...
            'rep2000','CC_',num2str(class_scenario),...
                'Vu',num2str(vU),'Vd',num2str(vD),...
                    'Vg',num2str(vG),'Vf',num2str(vF),'.mat'],'Stats')
    if print_contact_info
    text(0.5,0.78,['$\Delta T$ = ',sprintf('%0.1f',mean(Stats.mu.DT)/7),...
                    '(',sprintf('%0.2f',sqrt(mean(Stats.var.DT)/7)/mean(Stats.mu.DT/7)),') weeks'],...
    'Interpreter','latex','FontSize',font_size,'HorizontalAlignment','left',...
    'VerticalAlignment','middle');
    text(0.5,0.90,['$E\left[C(t_\textit{end})\right]$ = ',sprintf('%0.1f',mean(Stats.mu.CI)),...
                        '(',sprintf('%0.2f',sqrt(mean(Stats.var.CI))/mean(Stats.mu.CI)),')'],...
        'Interpreter','latex','FontSize',font_size,'HorizontalAlignment','left',...
    'VerticalAlignment','middle');
    if vF == 0 && vU == 0 && class_scenario == 50000
        RECB = 0;
    else
        RECB = 1-mean(Stats.mu.CI)/Baseline_CI;
    end
    text(.5,0.68,['RECI $= ',sprintf('%0.1f',RECB*100),'\%$'],...
                'HorizontalAlignment','left',...
                    'VerticalAlignment','middle',...
                        'Interpreter','latex','FontSize',font_size)
    end % END If print contact info
                
    hold off
    grid on
    ytickformat('%0.02f');
    xticks(0:15); xlim([0,15]);
    xticklabels(xlabs)
	yticks(0:.1:1)
    ylim([-0.01 1])
    set(gca,'TickLabelInterpreter','latex','FontSize',font_size+1)
end % END plot of sensitivity of infection parameters
    
%%%%%%%%%%%%%%%%%%%% PLOT SOCIAL CONTACT PARAMETERS %%%%%%%%%%%%%%%%%%%
if contains(param_set,'con') || contains(param_set,'soc')% Check if these should be plotted
    hSC = zeros(1,5);
    for i = 1:6
        sc = SC(i);
        hSC(i) = errorbar(tspan(tidx)/(24*7),SOBOL.indices(tidx,sc),SOBOL.error(tidx,sc),'LineWidth',LW(8-i),...
            'Color',colorblind(10-i,:),'MarkerSize',5); hold on
    end
    % CASE DETAILS in TEXT ON PLOT
    load(['sobol_indices_dt/Sobol_dt_sol',...
    num2str(y_sol_idx),...
        'sub',num2str(SubSampN),...
            'rep2000','CC_',num2str(class_scenario),...
                'Vu',num2str(vU),'Vd',num2str(vD),...
                    'Vg',num2str(vG),'Vf',num2str(vF),'.mat'],'Stats')
    if print_contact_info
    text(0.5,0.78,['$\Delta T$ = ',sprintf('%0.1f',mean(Stats.mu.DT)/7),...
                    '(',sprintf('%0.2f',sqrt(mean(Stats.var.DT)/7)/mean(Stats.mu.DT/7)),') weeks'],...
    'Interpreter','latex','FontSize',font_size,'HorizontalAlignment','left',...
    'VerticalAlignment','middle');
    text(0.5,0.90,['$E\left[C(t_\textit{end})\right]$ = ',sprintf('%0.1f',mean(Stats.mu.CI)),...
                        '(',sprintf('%0.2f',sqrt(mean(Stats.var.CI))/mean(Stats.mu.CI)),')'],...
        'Interpreter','latex','FontSize',font_size,'HorizontalAlignment','left',...
    'VerticalAlignment','middle');
    if vF == 0 && vU == 0 && class_scenario == 50000
        RECB = 0;
    else
        RECB = 1-mean(Stats.mu.CI)/Baseline_CI;
    end
    text(.5,0.68,['RECI $= ',sprintf('%0.1f',RECB*100),'\%$'],...
                'HorizontalAlignment','left',...
                    'VerticalAlignment','middle',...
                        'Interpreter','latex','FontSize',font_size)
    end % END If print contact info
    hold off
    grid on
    ytickformat('%0.02f');
    xticks(0:15); xlim([0,15]);
    xticklabels(xlabs)
	yticks(0:.1:1)
    ylim([-0.01 1])
    set(gca,'TickLabelInterpreter','latex','FontSize',font_size+1)
end % END plot of sensitivity of social contact parameters
    
% figure
%%%%%%%%%%%%%%%%%%%% PLOT INITIAL CONDITIONS %%%%%%%%%%%%%%%%%%%
if contains(param_set,'init')
    hICs = zeros(1,8);
    for i = 1:8
        ic = IC(i);
        hICs(i) = errorbar(tspan(tidx)/(24*7),SOBOL.indices(tidx,ic),SOBOL.error(tidx,ic),'LineWidth',LW(i),...
            'Color',colorblind(11-i,:),'MarkerSize',5); hold on
    end
    % CASE DETAILS in TEXT ON PLOT
    load(['sobol_indices_dt/Sobol_dt_sol',...
    num2str(y_sol_idx),...
        'sub',num2str(SubSampN),...
            'rep2000','CC_',num2str(class_scenario),...
                'Vu',num2str(vU),'Vd',num2str(vD),...
                    'Vg',num2str(vG),'Vf',num2str(vF),'.mat'],'Stats')
    if print_contact_info
    text(0.5,0.78,['$\Delta T$ = ',sprintf('%0.1f',mean(Stats.mu.DT)/7),...
                    '(',sprintf('%0.2f',sqrt(mean(Stats.var.DT)/7)/mean(Stats.mu.DT/7)),') weeks'],...
    'Interpreter','latex','FontSize',font_size,'HorizontalAlignment','left',...
    'VerticalAlignment','middle');
    text(0.5,0.90,['$E\left[C(t_\textit{end})\right]$ = ',sprintf('%0.1f',mean(Stats.mu.CI)),...
                        '(',sprintf('%0.2f',sqrt(mean(Stats.var.CI))/mean(Stats.mu.CI)),')'],...
        'Interpreter','latex','FontSize',font_size,'HorizontalAlignment','left',...
    'VerticalAlignment','middle');
    if vF == 0 && vU == 0 && class_scenario == 50000
        RECB = 0;
    else
        RECB = 1-mean(Stats.mu.CI)/Baseline_CI;
    end
    text(.5,0.68,['RECI $= ',sprintf('%0.1f',RECB*100),'\%$'],...
                'HorizontalAlignment','left',...
                    'VerticalAlignment','middle',...
                        'Interpreter','latex','FontSize',font_size)
    end % END If print contact info
    hold off
    grid on
    ytickformat('%0.02f');
    xticks(0:15); xlim([0,15]);
    xticklabels(xlabs)
    yticks(0:.1:1)
    ylim([-0.01 1])
    set(gca,'TickLabelInterpreter','latex','FontSize',font_size+1)
end % END plot of sensitivity of initial conditions

% Add descriptive y-label to figure
if print_ylabel
    ylabel(['\boldmath$v_{f,0} =  v_{g,0} = ',num2str(fac_grad_vax),'\%$',...
        char(10),'\boldmath$v_{u,0} = v_{d,0} = ',num2str(undg_vax),'\%$',...
            char(10),'Fraction of Variance'],...
                'Interpreter','latex','FontSize',font_size+2)
end