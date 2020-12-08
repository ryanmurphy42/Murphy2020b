function function_DISCRETE_realisations_generate_many(file_save_name, num_sims, k,a, prolif_law, ...
    beta_prolif, characteristics_include, timestop, time_record_vector, N_init,...
    mixing_pop_id ,eta,dt_discrete,q_init_condition,folder_name, EMT_law, EMT_param1, EMT_param2, EMT_param3,...
    init_chem_1_particle_distribution,D_chem_1,chem_1_inc,simulation_seed_start)

%% Set of many realisations of the discrete model and save the results

%% Initialise arrays for storage of many discrete realisations

t_hist_sim_run={};
x_hist_sim_run={};
k_hist_sim_run={};
a_hist_sim_run={};
G_hist_sim_run={};
cell_with_proliferation_event_hist_sim_run={};
cell_event_hist_sim_run = {};
cell_event_times_sim_run = {};
prolif_EMT_indicator_hist_sim_run = {};
cell_death_lr_hist_sim_run= {};
mixing_pop_id_hist_sim_run={};


%% Initialise variables
N=N_init;

%Initial density condition
if q_init_condition == 1
    % Mechanical equilibrium initial condition - cells all the same size
    x_current = [0;cumsum(a)]';
elseif q_init_condition ==2
    % All cells initially compressed to half of resting cell length
    x_current = [0;0.5*cumsum(a)]';
elseif q_init_condition == 3
    % All cells initially stretched to double resting cell length
    x_current = [0;2*cumsum(a)]';
end


for sim_run=1:num_sims
    
    %% Discrete numerical solver for one realisation
    
    rng(simulation_seed_start + sim_run,'twister') % Set the seed for realisation dependent on the simulation number
    
    tic
    [t_hist, x_hist_m, k_hist_m, a_hist_m, G_hist_m, cell_with_proliferation_event_hist,cell_with_EMT_event_hist,cell_event_hist, cell_event_times,...
        prolif_EMT_indicator_hist,mixing_pop_id_m_hist, chem_1_hist_m] = ...
        function_DISCRETE_generate_one(N, k, a, eta, beta_prolif, x_current, prolif_law, ...
        timestop, time_record_vector, mixing_pop_id,characteristics_include,dt_discrete, EMT_law, EMT_param1, EMT_param2, EMT_param3,...
        init_chem_1_particle_distribution, D_chem_1,chem_1_inc);
    toc
    
    
    %% Store variables
    
    t_hist_sim_run{sim_run} = t_hist;
    x_hist_sim_run{sim_run} = x_hist_m;
    k_hist_sim_run{sim_run} = k_hist_m;
    a_hist_sim_run{sim_run} = a_hist_m;
    G_hist_sim_run{sim_run} = G_hist_m;
    cell_with_proliferation_event_hist_sim_run{sim_run} = cell_with_proliferation_event_hist;
    cell_with_EMT_event_hist_sim_run{sim_run} = cell_with_EMT_event_hist;
    cell_event_hist_sim_run{sim_run} = cell_event_hist;
    cell_event_times_sim_run{sim_run} = cell_event_times;
    prolif_EMT_indicator_hist_sim_run{sim_run} = prolif_EMT_indicator_hist;
    mixing_pop_id_hist_sim_run{sim_run} = mixing_pop_id_m_hist;
    chem_1_hist_sim_run{sim_run} = chem_1_hist_m;
    
    
end


%% Save .mat file with outputs

save([pwd '/' folder_name '/' file_save_name],'-v7.3',...
    't_hist_sim_run',...
    'x_hist_sim_run',...
    'k_hist_sim_run',...
    'a_hist_sim_run',...
    'eta',...
    'num_sims',...
    'mixing_pop_id_hist_sim_run',...
    'cell_with_proliferation_event_hist_sim_run',...
    'cell_with_EMT_event_hist_sim_run',...
    'cell_event_hist_sim_run',...
    'cell_event_times_sim_run',...
    'prolif_EMT_indicator_hist_sim_run',...
    'cell_death_lr_hist_sim_run',...
    'G_hist_sim_run',...
    'timestop',...
    'chem_1_hist_sim_run',...
    'time_record_vector')
