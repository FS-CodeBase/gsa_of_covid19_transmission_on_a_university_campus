function [BPLOT,LEG] = barplot_sobol_indices(FirstOrderIdx,TotalOrderIdx,...
                    prms_str,leg_bool,axis_lab_bool)
LEG = 0;
C = linspecer(2);

    SobolMat = [FirstOrderIdx.indices;TotalOrderIdx.indices];
    BPLOT = bar(SobolMat','grouped','FaceColor','flat'); hold on
    for c_idx = 1:size(SobolMat',2)
        BPLOT(c_idx).CData = repmat(C(c_idx,:),[numel(FirstOrderIdx.indices),1]);
    end
    ylim([0 1])
    yticks(0:.1:1)
    xticks(1:numel(prms_str))
    xticklabels(prms_str)
    
    % Plot error bars
    ERRS = [FirstOrderIdx.error;TotalOrderIdx.error];
    errorbar([BPLOT(1).XEndPoints;BPLOT(2).XEndPoints],...
            SobolMat,ERRS,'k','linestyle','none','LineWidth',1);  
    hold off
    if leg_bool
        LEG = legend({'FO','TE'},'Location',...
                        'northwest',...
                        'Interpreter','latex','FontSize',12,...
                        'NumColumns',2);
    end
    if axis_lab_bool
        ylabel('Fraction of Variance','Interpreter','latex','FontSize',13)
        xlabel('Model Factor','Interpreter','latex','FontSize',13)
    end
    grid on
    set(gcf,'Position',[10 10 600 320])
    set(gca,'TickLabelInterpreter','latex','FontSize',13)
    ytickformat('%0.02f');
