function function_PLOTS_ctm_total_cell_number(t_hist, q_hist,L_hist,L,dz, filepath_load_figs,filepath_save_figs, plot_times, plot_ctm_over_dis)

%% Plot the total cell number from the continuum model

mass_hist = [];
for time_vector_plot_counter = 1:length(plot_times)
    tcomparing = plot_times(time_vector_plot_counter);
    [~, first_compare_index]   = min(abs(t_hist(1:end) -tcomparing));
    mass_hist= [mass_hist; trapz((0:dz:1)*L_hist(:,first_compare_index),q_hist(:,first_compare_index)')];
end


%% Plots


if plot_ctm_over_dis==1 %plot ctm model over the discrete model
    
    %open the discrete figure
    openfig([filepath_load_figs '\' 'N_evolution_dis.fig'])
    hold on
    
    %plot the total cell number for population 1 from the continuum model
    plot(plot_times',mass_hist,'g','linewidth',3)
    
    legend('Dis','Ctm')
    
    title('Evolution of N(t) - dis & ctm')
    xlabel('t')
    ylabel('N(t)')
    
    
    %save figure
    print(gcf,'-depsc2',[filepath_save_figs '\' 'N_evolution_dis_ctm' '.eps'])
    saveas(gcf,[filepath_save_figs '\' 'N_evolution_dis_ctm' '.fig'])
    saveas(gcf,[filepath_save_figs '\' 'N_evolution_dis_ctm' '.jpg'])
    
    
else
    
    %plot the total cell number
    figure
    plot(plot_times',mass_hist,'m','linewidth',3)
    title('Evolution of N(t) - ctm')
    xlabel('t')
    ylabel('N(t)')
    
    
    %Save figures
    print(gcf,'-depsc2',[filepath_save_figs '\' 'N_evolution_ctm' '.eps'])
    saveas(gcf,[filepath_save_figs '\' 'N_evolution_ctm' '.fig'])
    saveas(gcf,[filepath_save_figs '\' 'N_evolution_ctm' '.jpg'])
    
end

end