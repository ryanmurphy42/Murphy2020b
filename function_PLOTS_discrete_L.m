function function_discrete_PLOTS_L(t_hist_sim_run, x_hist_sim_run, filepath_save_figs,num_sims, plot_times,mixing_pop_id_hist_sim_run)

%% Plot the mean and standard deviation of the interface position


%Create a matrix to store the N for each sim at the times
intpos_plot = zeros(num_sims,size(plot_times,2));

%For each time detmerine the interface position
for intpos_plot_timeloop = 1:size(plot_times,2)
    for sim_run1=1:num_sims
        %find the closest time
        timing_diff = t_hist_sim_run{sim_run1}(1:end) - plot_times(intpos_plot_timeloop);
        timing_diff(timing_diff <-0.0001) = inf;
        if sum(timing_diff==inf) == size(timing_diff,2)
            time_index_1 = size(timing_diff,2);
        else
            [~, time_index_1] = min(timing_diff);
        end
        
        %use the population id to determine where the interface position is
        nocellpop1 = length(cell2mat(mixing_pop_id_hist_sim_run{sim_run1}(time_index_1))) - sum(cell2mat(mixing_pop_id_hist_sim_run{sim_run1}(time_index_1)));
        if nocellpop1 > 0
            x_int_time =  cell2mat(x_hist_sim_run{sim_run1}(time_index_1));
            intpos_plot(sim_run1,intpos_plot_timeloop) = x_int_time(nocellpop1+1);
        else
            intpos_plot(sim_run1,intpos_plot_timeloop) = 0;
        end
        
    end
end

%Determine the mean for each time point
intpos_plot_mean = mean(intpos_plot,1);

%Determine the standard deviation for each time point
intpos_plot_std = std(intpos_plot,0,1);

%plot figures
figure
hold on
errorbar(plot_times,intpos_plot_mean,intpos_plot_std,'-s','MarkerSize',5,'MarkerEdgeColor','blue','MarkerFaceColor','blue','Color','blue','LineStyle','-')
xlabel('Time')
ylabel('L(t)')
title('Evolution of L(t)')

box on

%save figures
print(gcf,'-depsc2',[filepath_save_figs  'Evolution_L_dis.eps']);
saveas(gcf, [filepath_save_figs  'Evolution_L_dis.fig'])
saveas(gcf, [filepath_save_figs  'Evolution_L_dis.jpg'])

end
