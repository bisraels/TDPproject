% function cP_plotter(P, time)


%-------------------------------------------------------------------------
% Plot the Conditional Probabilities
%-------------------------------------------------------------------------
if plotCondProbMode == 1
    if plotMode == 1
        
        figure(3);
        clf;
        set(gcf,'Color','w');
        hold on;
        logx;
        set(gca,'xscale','log');
        
        [sample_description, ~] = sample_descriptionGetter();
        title_str = ['Conditional Proababilities: P_{i\rightarrowj}(t)' ...
            10 sample_description];
        title(title_str,'fontsize',14);
        %         title('Conditional Proababilities','FontSize',18)
        xlabel('Time (sec)','FontSize',14);
        ylabel('Probability','FontSize',14);
        
        % FOR THE THREE STATE CYCLICAL MODEL
        % 1 := C    Coiled (collapsed?)
        % 2 := PC   Partially Coiled (collapsed?)
        % 3 := E    Extended
        
%         state_1_color = 'k';%[0, 0.4470, 0.7410];        % "C"
%         state_2_color = 'b';%[0.8500, 0.3250, 0.0980];   % "I"
%         state_3_color = 'r';%[0.9290, 0.6940, 0.1250];   % "FE"
        state_color = ['k','b','r'];  % Plot based on i
        state_style  = ['-',':','--']; % Plot based on j
        timeUB = 1e-2;
        Pshape = size(P);
        N = Pshape(1);
        
        for j = 1:N
        for  i = 1:N
            plot(time, reshape(P(j,i,:),[1,length(time)]),'color',state_color(i),'LineStyle',state_style(j));
        end
        end
    end
end
        
        
        
        
        
        