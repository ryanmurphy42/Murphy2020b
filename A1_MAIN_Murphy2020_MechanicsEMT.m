%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 2020--Murphy et al--The role of mechanical interactions in EMT

% Main Script

% Key algorithms used to generate paper figures in Figure 3a,c and Figure 6b,e,h.

% Please run individual sections separately due to the time required to generate results.

% Discrete realisations in code set to 50 for speed 
% (manuscript uses 2000 discrete realisations and the corresponding .mat files are available on GitHub)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 0) Reset Matlab

clear
clc

%% 1) Loop through Figure 3 and Figure 6

% FigureID=3 for Figures 3a,c
% FigureID=6 for Figures 6b,e,h

for FigureID = [3,6]


%% 2) Obtain the parameterse

% Parameters                 -- Description

% simulation_id              -- ID for simulation to save files.

% k1                         -- cell stiffness
% eta1                       -- mobility coefficient
% a1                         -- resting cell length

% prolif_law                 -- Proliferation mechanism: 1 - cell-length-independent ; 2 - cell-length-dependent
% beta_1                     -- cell-length-independent proliferation rate

% N_cells_init               -- initial number of cells
% q_init_condition           -- 1 - mechanical equilbrium, 2 - uniform compression to half of resting cell length, 3 - uniform stretch to double resting cell length,

% EMT_law                    -- EMT mechanism: 1 - chemically-independent, 2 - chemically-dependent
% omega (MT1)                -- chemically-independent EMT rate
% phi   (0 <= phi <= 1)      -- phi/omega - time to reach chemical threshold when cell at resting cell length, (1-phi)/omega - average time for boundary cell to detach when have invasive phenotype

% D_chem_1                   -- Diffusivity of chemical
% chem_1_inc                 -- Is chemical included in simulation
% chem_1_init_condition      -- Initial condition for chemical concentration

% timestop                   -- Time stop for numerical simulations
% dt_discrete                -- Time step for realisations of discrete model, usually set to 1e-4, for D=1 set to 1e-5.
% dt_ctm                     -- Biggest time step for continuum model - time step is adaptive in the continuum model
% dz_ctm                     -- spatial step in continuum model discretisation
% err_tol_ctm                -- error tolerance for Newton solver for continuum model

% num_sims_char              -- Number of simulations of discrete model to produce characteristic plots
% num_sims_realiz            -- Number of simulations of discrete model for averaging to compare with continuum model
% simulation_seed_start      -- Simulation seed - used if want to repeat a simulation

% time_record_vector_char    -- Time points when to record for characteristic plots
% time_record_vector_realiz  -- Time points when to record for discrete plots
% time_record_vector_ctm     -- Time points when to record for continuum plots


% plot_dis                   -- 1 - plot the results from the discrete model, 0 - do not plot
% plot_ctm                   -- 1 - plot the results from the continuum model, 0 - do not plot
% plot_ctm_over_dis          -- 1 - plot the continuum results over the discrete model results, 0 - plot continuum results without discrete model


if FigureID == 3
    %% For Figure 3a,c
    
    simulation_id='Fig3';
    
    k1=1;
    a1 =0.1;
    eta1=1;
    
    prolif_law=1;
    beta_1=0.0025;
    
    N_cells_init=20;
    q_init_condition=1;
    
    EMT_law=1;
    omega=0.1;
    phi=0.9;
    
    chem_1_inc=0;
    D_chem_1=1; % Not used in this simulation since chem_1_inc=0
    chem_1_init_condition=1;  % Not used in this simulation chem_1_inc=0
    
    timestop=200;
    dt_discrete=1e-4;
    dt_ctm=1e-3;
    dz_ctm=1e-5;
    err_tol_ctm = 1e-8;
    
    num_sims_char=1;
    num_sims_realiz=50;
    simulation_seed_start=1;
    
    time_record_vector_char=0:0.5:timestop;
    time_record_vector_realiz=0:1:timestop;
    time_record_vector_ctm=0:1:timestop;
    
    plot_dis=1;
    plot_ctm=1;
    plot_ctm_over_dis=1;
    
elseif FigureID == 6
    %% For Figure 6b,e,h
    
    simulation_id='Fig6';
    
    k1=1;
    a1 =0.1;
    eta1=1;
    
    prolif_law=1;
    beta_1=0.0025;
    
    N_cells_init=20;
    q_init_condition=1;
    
    EMT_law=2;
    omega=0.1;
    phi=0.9;
    
    chem_1_inc=1;
    D_chem_1=1e-2;
    chem_1_init_condition=1;
    
    timestop=400;
    dt_discrete=1e-4;
    dt_ctm=1e-3;
    dz_ctm=1e-5;
    err_tol_ctm = 1e-8;
    
    num_sims_char=1;
    num_sims_realiz=50;
    simulation_seed_start=1;
    
    time_record_vector_char=0:0.5:timestop;
    time_record_vector_realiz=0:1:timestop;
    time_record_vector_ctm=0:1:timestop;
    
    plot_dis=1;
    plot_ctm=1;
    plot_ctm_over_dis=1;
    
end

%% 3) Run the discrete model (produces characteristics and num_sims_realiz realisations of the discrete model)

tic
function_DISCRETE_0_realisations_main(simulation_id,N_cells_init,  beta_1,k1,a1,eta1,omega,phi,...
    num_sims_char,num_sims_realiz,prolif_law,D_chem_1,timestop,q_init_condition,EMT_law,dt_discrete,...
    time_record_vector_char,time_record_vector_realiz,chem_1_inc,simulation_seed_start )
toc

% For Figure 3 - approximately 20 seconds for one realisation (manuscript uses 2000 simulations).
% For Figure 6 (including the chemical) approximately 1-1.5 minutes per realisation (manuscript uses 2000 simulations).
% For D=1 run with dt_discrete=1e-5, which is minutes slower per realisation.


%% 4) Run the continuum model

tic
function_CTM_0_main(simulation_id,N_cells_init,  beta_1,k1,a1,eta1,omega,...
    prolif_law,phi,D_chem_1,timestop,q_init_condition,EMT_law,dt_ctm,...
    time_record_vector_ctm,chem_1_inc, chem_1_init_condition,dz_ctm,err_tol_ctm)
toc

% For Figure 3 - minutes to run
% For Figure 6 - order of hours to run


%% 5) Plot the Figures

tic
function_PLOTS_0_dis_ctm(simulation_id, num_sims_realiz,chem_1_inc,plot_dis,plot_ctm,plot_ctm_over_dis)
toc

end
