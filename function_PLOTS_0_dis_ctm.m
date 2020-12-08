function function_PLOTS_dis_ctm(simulation_id, num_sims_realiz,chem_1_inc,plot_dis,plot_ctm,plot_ctm_over_dis)

%% PLOTS
% 1.1) Discrete characteristics
% 1.2) Discrete - evolution of average of many discrete realisations for total cell number, N(t), and boundary position, L(t)
% 2) Continuum model - evolution of total cell number, N(t), and boundary position, L(t)

if plot_dis == 1
    %% 1.1) PLOT - DISCRETE - characteristics
    
    filepath_save_figs = [pwd '\' simulation_id '_Discrete\'];
    load([filepath_save_figs 'Results_Sim' simulation_id '_Discrete_Characteristics.mat']);
    
    %check that the all the data is available to plot characteristics without error. If the below returns a result decreate dt.
    function_PLOTS_characteristic_plot_check(num_sims,t_hist_sim_run,x_hist_sim_run);
    
    %define figure properties
    
    %plot the characteristics
    for sim_run = 1:num_sims
        
        cell_with_proliferation_event_hist=cell_with_proliferation_event_hist_sim_run{sim_run};
        cell_with_EMT_event_hist=cell_with_EMT_event_hist_sim_run{sim_run};
        cell_event_hist=cell_event_hist_sim_run{sim_run};
        prolif_EMT_indicator_hist=prolif_EMT_indicator_hist_sim_run{sim_run};
        t_hist=t_hist_sim_run{sim_run};
        x_hist=x_hist_sim_run{sim_run};
        k_hist =k_hist_sim_run{sim_run};
        a_hist =a_hist_sim_run{sim_run};
        G_hist =G_hist_sim_run{sim_run};
        mixing_pop_id_hist =mixing_pop_id_hist_sim_run{sim_run};
        chem_1_hist = chem_1_hist_sim_run{sim_run};
        
        %density colouring
        colouring =1;
        function_PLOTS_characteristics(t_hist, x_hist, k_hist, mixing_pop_id_hist, a_hist,...
            filepath_save_figs,  colouring, cell_with_proliferation_event_hist,cell_with_EMT_event_hist,...
            cell_event_hist,prolif_EMT_indicator_hist,sim_run,chem_1_hist,chem_1_inc,timestop)
        
        if chem_1_inc == 1
            %chemical 1 colouring
            colouring =2;
             function_PLOTS_characteristics(t_hist, x_hist, k_hist, mixing_pop_id_hist, a_hist,...
            filepath_save_figs,  colouring, cell_with_proliferation_event_hist,cell_with_EMT_event_hist,...
            cell_event_hist,prolif_EMT_indicator_hist,sim_run,chem_1_hist,chem_1_inc,timestop)
       
        end
        
        close all
    end
    
    
    %% 1.2) PLOT - DISCRETE - evolution of average of many discrete realisations for total cell number, N(t), and boundary position, L(t)
    
    filepath_save_figs = [pwd '\' simulation_id '_Discrete\'];
    load([filepath_save_figs 'Results_Sim' simulation_id '_Discrete_R' num2str(num_sims_realiz) '.mat']);
    
    if exist('time_record_vector','var') == 0
        plot_times = 0:5:timestop;
    else
        plot_times = time_record_vector;
    end
    
    % 2.5.3) Evolution of N(t)
    function_PLOTS_discrete_total_cell_number(t_hist_sim_run, x_hist_sim_run, filepath_save_figs,num_sims_realiz,...
        plot_times,mixing_pop_id_hist_sim_run)
    
    % 2.5.4) Evolution of L(t)
    function_PLOTS_discrete_L(t_hist_sim_run, x_hist_sim_run, filepath_save_figs,num_sims_realiz, plot_times,mixing_pop_id_hist_sim_run)
    
end


%% 2) Overlay continuum density snapshots, total cell number discrete, and interface position
if plot_ctm == 1
    filepath_save_figs = [pwd '\' simulation_id '_Continuum\'];
    load([filepath_save_figs 'Results_Sim' simulation_id '_Continuum.mat']);
    
    file_save_name = ['Results_Sim' simulation_id '_Continuum'];
    
    folder_name = [simulation_id '_Continuum'];
    if ~exist([folder_name'], 'dir')
        mkdir([folder_name]);
    end
    
    filepath_load_figs = [pwd '\' simulation_id '_Discrete\'];
    
    
    if exist('time_record_vector_ctm','var') == 0
        plot_times = 0:5:timestop;
    else
        plot_times = time_record_vector_ctm;
    end
    
    % 2.6.3) Total cell number
    function_PLOTS_ctm_total_cell_number(t_hist, q_hist,L_hist,L,dz, filepath_load_figs,filepath_save_figs, time_record_vector_ctm, plot_ctm_over_dis)
    
    % 2.6.4) Interface position
    function_PLOTS_ctm_L(t_hist,L_hist, filepath_load_figs,filepath_save_figs, time_record_vector_ctm, plot_ctm_over_dis)
    
end

end