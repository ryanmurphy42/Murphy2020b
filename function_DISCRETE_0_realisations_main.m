function function_DISCRETE_0_realisations_main(simulation_id,N_cells_init,  beta_1,k1,a1,eta1,omega,phi,...
    num_sims_char,num_sims_realiz,prolif_law,D_chem_1,timestop,q_init_condition,EMT_law,dt_discrete,...
    time_record_vector_char,time_record_vector_realiz,chem_1_inc,simulation_seed_start)

%% Main script to generate discrete realisations

%% Create folder to save results

folder_name = [simulation_id '_Discrete'];
if ~exist([folder_name'], 'dir')
    mkdir([folder_name]);
end

%% 3.1) Set up parameters

N_init=N_cells_init+1; %number of cells + 1 (the plus 1 is for the cell which represents space)
mixing_pop_id = [0*ones(N_init-1,1); 1]; %number of cells (the plus 1 is for the cell which represents space)
k=[k1*ones(N_cells_init,1);0]; %0 is for cell which represents space
a=[a1*ones(N_cells_init,1);0]; %0 is for cell which represents space

beta_prolif=[beta_1*ones(N_cells_init,1);0]; % proliferation rate for each cell

init_chem_1_particle_distribution =  zeros(N_init,1); % zero chemical in any of the cells initially

%%  EMT chemical

% Initialse parameters for EMT
EMT_param1 = omega; %EMT constant rate
EMT_param2 = -1; % EMT concentration threshold
EMT_param3 = -1; % EMT source rate

particles_EMTcrit = 50;
time_EMT1 = (phi)/omega;
time_EMT2 = (1-phi)/omega;

chem_EMTcrit = particles_EMTcrit/a(1); % this is C in the manuscript.
chem_part_source = particles_EMTcrit/time_EMT1; % this is S in the paper
chem_EMTcrateabovecrit = 1/time_EMT2;

% Reset the EMT variables

% EMT rate above threshold
if chem_1_inc == 1 % include the chemical
    if EMT_law == 2 % chemically dependent EMT
        EMT_param1 = chem_EMTcrateabovecrit;
    end
    EMT_param2 = chem_EMTcrit; % EMT concentration threshold
    EMT_param3 = chem_part_source; % EMT source rate 
end

%% 3.2) Generate realisations to plot characteristic diagrams

characteristics_include=1; % include timepoints when a cell event (proliferation or EMT) occurs

file_save_name = ['Results_Sim' simulation_id '_Discrete_Characteristics'];

tic
function_DISCRETE_realisations_generate_many(file_save_name, num_sims_char, k,a, prolif_law, ...
    beta_prolif, characteristics_include, timestop, time_record_vector_char, N_init,...
    mixing_pop_id ,eta1,dt_discrete,q_init_condition,folder_name, EMT_law, EMT_param1, EMT_param2, EMT_param3,...
    init_chem_1_particle_distribution,D_chem_1,chem_1_inc,simulation_seed_start);
toc

%% 3.3) Generate discrete realizations to plot averaged data

characteristics_include=0;

file_save_name = ['Results_Sim' simulation_id '_Discrete_' num2str(simulation_seed_start) '_R' num2str(num_sims_realiz)];

tic
function_DISCRETE_realisations_generate_many(file_save_name, num_sims_realiz, k,a, prolif_law, ...
    beta_prolif, characteristics_include, timestop, time_record_vector_realiz, N_init,...
    mixing_pop_id ,eta1,dt_discrete,q_init_condition,folder_name, EMT_law, EMT_param1, EMT_param2, EMT_param3,...
    init_chem_1_particle_distribution,D_chem_1,chem_1_inc,simulation_seed_start);
toc


end