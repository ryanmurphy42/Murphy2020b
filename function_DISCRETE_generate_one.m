function    [t_hist, x_hist_m, k_hist_m, a_hist_m, G_hist_m, cell_with_proliferation_event_hist,cell_with_EMT_event_hist,cell_event_hist, cell_event_times,...
    prolif_EMT_indicator_hist,mixing_pop_id_m_hist, chem_1_hist_m] = ...
    function_DISCRETE_generate_one(N, k, a, eta, beta_prolif, x_current, prolif_law, ...
    timestop, time_record_vector, mixing_pop_id,characteristics_include,dt_discrete, EMT_law, EMT_param1, EMT_param2, EMT_param3,...
    init_chem_1_particle_distribution, D_chem_1,chem_1_inc)


%% Generate one discrete realisation



% convert the times to record at timesteps
time_record_vector_dt = round(time_record_vector/dt_discrete,0);


%% Initialise arrays and vectors for storage

t_hist = [];
x_hist_m = {};
k_hist_m = {};
a_hist_m = {};
G_hist_m = {};
prolif_times=[];
EMT_times=[];
cell_event_times=[];
cell_with_proliferation_event_hist = [];
cell_with_EMT_event_hist = [];
cell_event_hist = [];
prolif_EMT_indicator_hist =[];

chem_1_hist_m = [];
chem_1_m = init_chem_1_particle_distribution./diff(x_current)';
chem_1_m(end)=0;

%% Initialise variables
t=0;
t_record_loop=1;
t_steprec=0;

% cell properties
N_cells=N; % N_cells is number of actual cells + 1.
m_springs_per_cell=1;
N_springs = N_cells*m_springs_per_cell;

%convert the positions to springs per cell
x_current_m=zeros(N_cells*m_springs_per_cell + 1,1);
k_m=zeros(N_cells*m_springs_per_cell,1);
a_m=zeros(N_cells*m_springs_per_cell,1);
beta_prolif_m=zeros(N_cells*m_springs_per_cell,1);
mixing_pop_id_m=zeros(N_cells*m_springs_per_cell,1);


%convert to spring properties (m_springs_per_cell = 1) so cells=springs.
cell_spring_index = 0;
for ii=1:N_cells
    for jj=1:m_springs_per_cell
        cell_spring_index = cell_spring_index +1;
        x_current_m(cell_spring_index) = x_current(ii) + (x_current(ii+1) -x_current(ii))*((jj-1)/m_springs_per_cell);
        k_m(cell_spring_index) = k(ii)*m_springs_per_cell;
        a_m(cell_spring_index) = a(ii)/m_springs_per_cell;
        beta_prolif_m(cell_spring_index) = beta_prolif(ii)/m_springs_per_cell;
        mixing_pop_id_m(cell_spring_index) = mixing_pop_id(ii);
    end
end
x_current_m(N_cells*m_springs_per_cell + 1) = sum(a_m)*100; %set maximum value that wont be reached for plotting.
eta_m = eta/m_springs_per_cell;


%% save the spring properties
t_hist=[t_hist,t];
x_hist_m{t_record_loop} = x_current_m;
k_hist_m{t_record_loop} = k_m;
a_hist_m{t_record_loop} = a_m;
mixing_pop_id_m_hist{t_record_loop} = mixing_pop_id_m;
chem_1_hist_m{t_record_loop} = chem_1_m;
compartments_per_cell_overall = ones(N_cells,1);

% proliferation and EMT propensities
[G_m,~,~] = function_DISCRETE_proliferation_EMT_rates(N_cells,m_springs_per_cell, prolif_law,beta_prolif_m,x_current_m,EMT_law, EMT_param1, EMT_param2, chem_1_m(N_cells-1),a_m);
G_hist_m{t_record_loop} = G_m;



while (t<timestop)
    
    
    %Run the numerical simulation
    
    
    
    if N_cells == 1 %
        
        t_record_loop=t_record_loop+1;
        t_hist=[t_hist,t];
        x_hist_m{t_record_loop} = x_current_m;
        k_hist_m{t_record_loop} = k_m/m_springs_per_cell;
        a_hist_m{t_record_loop} = a_m*m_springs_per_cell;
        G_hist_m{t_record_loop} = G_m*m_springs_per_cell;
        chem_1_hist_m{t_record_loop} = chem_1_m*m_springs_per_cell;
        mixing_pop_id_m_hist{t_record_loop} = mixing_pop_id_m;
        t=timestop+1;
        
    else
        
        %% Update position
        x_prev_m = x_current_m;
        for i=2:N_springs
            x_current_m(i) = x_prev_m(i) + dt_discrete*(1/eta_m)*(    k_m(i-1)*(a_m(i-1) - (abs(x_prev_m(i) - x_prev_m(i-1))))*((x_prev_m(i) - x_prev_m(i-1))/(abs(x_prev_m(i) - x_prev_m(i-1)))) ...
                + (k_m(i)*(a_m(i) - (abs(x_prev_m(i) - x_prev_m(i+1)) )))*((x_prev_m(i) - x_prev_m(i+1))/(abs(x_prev_m(i) - x_prev_m(i+1)))) );
        end
        
        
        %% Diffusion
        
        if chem_1_inc == 1  % chemical included
            
            cell_source_rates_chem_1 = zeros(N_cells,1);
            if N_cells > 1
                cell_source_rates_chem_1(N_cells-1) = (1/(x_current_m(N_cells)-x_current_m(N_cells-1)))*EMT_param3;
            end
            
            [ chem_1_m] =  function_DISCRETE_chemical_concentration(D_chem_1,chem_1_m,x_current_m,x_prev_m,cell_source_rates_chem_1,dt_discrete);
            
        end
        
        
        %% Update time
        t=t+dt_discrete;
        
        t_steprec = t_steprec + 1;
        
        %% Save data
        
        if sum(time_record_vector_dt == t_steprec) > 0
            
            t_record_loop=t_record_loop+1;
            t_hist=[t_hist,t];
            x_hist_m{t_record_loop} = x_current_m;
            k_hist_m{t_record_loop} = k_m/m_springs_per_cell;
            a_hist_m{t_record_loop} = a_m*m_springs_per_cell;
            G_hist_m{t_record_loop} = G_m*m_springs_per_cell;
            mixing_pop_id_m_hist{t_record_loop} = mixing_pop_id_m;
            
            chem_1_hist_m{t_record_loop} = chem_1_m*m_springs_per_cell;
            
            fprintf('\n t= %f; Ncells = %f; length = %f; max chem conc = %f\n', t,N_cells-1,x_current_m(end-1), max(abs(chem_1_m)));
            
            
        end
        
        
    end
    
    if N_cells == 1 % If N_cells=1 then no real cells left - this cell corresponds to a fake cell for simulation purposes.
        cell_event_happening =0;
    else
        x_new_spring_positions = 0:(1/(2*m_springs_per_cell)):1;
        x_new_spring_positions_trunc = x_new_spring_positions(2:end-1)';
        
        %% update the propensities
        
        [G_m,~,~] = function_DISCRETE_proliferation_EMT_rates(N_cells,m_springs_per_cell, prolif_law,beta_prolif_m,x_current_m,EMT_law, EMT_param1, EMT_param2, chem_1_m(N_cells-1),a_m);
        
        
        %% calculate if there is a cell event (proliferation or EMT)
        
        G_m_cumsum = cumsum(G_m)/sum(G_m);
        
        cell_event_happening =0;
        r1 = rand(1);
        if  r1 < sum(G_m*dt_discrete)
            cell_event_happening =1;
        end
        
    end
    
    
    %% save data
    if cell_event_happening == 1
        if characteristics_include==1
            if sum(time_record_vector_dt == t_steprec) == 0 %if havent already recorded this timestep then record here as there is an event
                time_record_vector = sort([time_record_vector,t]);
                time_record_vector_dt = round(time_record_vector/dt_discrete,0);
                
                t_record_loop=t_record_loop+1;
                t_hist=[t_hist,t];
                x_hist_m{t_record_loop} = x_current_m;
                k_hist_m{t_record_loop} = k_m/m_springs_per_cell;
                a_hist_m{t_record_loop} = a_m*m_springs_per_cell;
                G_hist_m{t_record_loop} = G_m*m_springs_per_cell;
                mixing_pop_id_m_hist{t_record_loop} = mixing_pop_id_m;
                
                chem_1_hist_m{t_record_loop} = chem_1_m*m_springs_per_cell;
            end
            
        end
    end
    
    if  cell_event_happening ==1
        
        r2 = rand(1);
        %% Determine whether have cell proliferation OR cell death OR EMT
        
        G_m_cum_sum_inf = G_m_cumsum(1:end) -r2;
        
        G_m_cum_sum_inf(G_m_cum_sum_inf < 0) = inf;
        
        [~ , spring_event]   = min(G_m_cum_sum_inf);
        
        %determine cell event
        
        cell_event = floor((spring_event-1)/m_springs_per_cell)+ 1;
        
        if spring_event <=  N_springs %cell proliferation event
            
            if cell_event == 1
                
                %Update the positions and enforce the cell proliferation event by dividing the cell equally
                %Relabel the x positions
                
                x_current_m = [x_current_m(1);
                    x_current_m(1) + x_new_spring_positions_trunc*(x_current_m(cell_event*m_springs_per_cell+1)-x_current_m(1));
                    x_current_m(cell_event*m_springs_per_cell+1:end)];
                
                
            elseif cell_event == N_cells
                
                %Update the positions and enforce the cell proliferation event by dividing the cell equally
                %Relabel the x positions
                
                x_current_m = [x_current_m(1:(cell_event-1)*m_springs_per_cell +1);
                    x_current_m((cell_event-1)*m_springs_per_cell +1) + x_new_spring_positions_trunc*(x_current_m(end)-x_current_m((cell_event-1)*m_springs_per_cell +1)); ...
                    x_current_m(end)];
                
            else
                
                %Update the positions and enforce the cell proliferation event by dividing the cell equally
                %Relabel the x positions
                
                x_current_m = [x_current_m(1:(cell_event-1)*m_springs_per_cell +1);
                    x_current_m((cell_event-1)*m_springs_per_cell + 1) + x_new_spring_positions_trunc*(x_current_m(cell_event*m_springs_per_cell+1)-x_current_m((cell_event-1)*m_springs_per_cell+ 1));
                    x_current_m(cell_event*m_springs_per_cell+1:end)];
                
            end
            
            %Update k
            k_m = [k_m(1:cell_event*m_springs_per_cell); k_m(cell_event*m_springs_per_cell)*ones(m_springs_per_cell,1) ;k_m(cell_event*m_springs_per_cell+1:end)];
            
            %Update a
            a_m = [a_m(1:cell_event*m_springs_per_cell); a_m(cell_event*m_springs_per_cell)*ones(m_springs_per_cell,1) ;a_m(cell_event*m_springs_per_cell+1:end)];
            
            %Update beta
            beta_prolif_m = [beta_prolif_m(1:cell_event*m_springs_per_cell); beta_prolif_m(cell_event*m_springs_per_cell)*ones(m_springs_per_cell,1) ...
                ;beta_prolif_m(cell_event*m_springs_per_cell+1:end)];
            
            
            if cell_event  == 1
                % Update the chemical concentrations - chemical 1 - cell level %divide the chemical equally across the cell independent of how it is distributed within the cell if there are smaller subcompartments
                chem_1_m = [chem_1_m(cell_event*m_springs_per_cell)*ones(m_springs_per_cell*2,1) ...
                    ;chem_1_m((cell_event+1)*m_springs_per_cell:end)];
            elseif cell_event == N_cells
                % Update the chemical concentrations - chemical 1 - cell level %divide the chemical equally across the cell independent of how it is distributed within the cell if there are smaller subcompartments
                chem_1_m = [chem_1_m(1:(cell_event-1)*m_springs_per_cell); chem_1_m(cell_event*m_springs_per_cell)*ones(m_springs_per_cell*2,1) ];
            else
                % Update the chemical concentrations - chemical 1 - cell level %divide the chemical equally across the cell independent of how it is distributed within the cell if there are smaller subcompartments
                chem_1_m = [chem_1_m(1:(cell_event-1)*m_springs_per_cell); chem_1_m(cell_event*m_springs_per_cell)*ones(m_springs_per_cell*2,1) ...
                    ;chem_1_m((cell_event+1)*m_springs_per_cell:end)];
            end
            
            %Update mixing_pop_id
            mixing_pop_id_m = [mixing_pop_id_m(1:cell_event*m_springs_per_cell); mixing_pop_id_m(cell_event*m_springs_per_cell)*ones(m_springs_per_cell,1) ;mixing_pop_id_m((cell_event)*m_springs_per_cell + 1:end)];
            
            
            %Update the number of cells
            N_cells = N_cells + 1;
            
            N_springs = N_springs + m_springs_per_cell;
            
            cell_with_proliferation_event_hist = [cell_with_proliferation_event_hist,cell_event];
            
            cell_event_hist = [cell_event_hist,cell_event];
            
            prolif_EMT_indicator_hist = [prolif_EMT_indicator_hist, 1]; %1 correspond cell proliferation
            
            prolif_times = [prolif_times, t];
            
            cell_event_times = [cell_event_times,t];
            
            
        else %EMT
            % final real cell is lost
            
            cell_event_EMT = N_cells-1; %last real cell
            
            if m_springs_per_cell == 1
                x_current_m = [x_current_m(1:cell_event_EMT);
                    x_current_m(end)];
            else
                x_current_m = [x_current_m(1:((cell_event_EMT-1)*m_springs_per_cell+1));
                    x_current_m(((cell_event_EMT-1)*m_springs_per_cell+1)) + x_new_spring_positions_death_coal_trunc*(x_current_m(end)-x_current_m(((cell_event_EMT-1)*m_springs_per_cell+1)));
                    x_current_m(end)];
            end
            
            %Update k
            k_m = [k_m(1:(cell_event_EMT-1)*m_springs_per_cell) ;k_m((cell_event_EMT)*m_springs_per_cell + 1:end)];
            
            %Update a
            a_m = [a_m(1:(cell_event_EMT-1)*m_springs_per_cell) ;a_m((cell_event_EMT)*m_springs_per_cell +1:end)];
            
            %Update beta
            beta_prolif_m = [beta_prolif_m(1:(cell_event_EMT-1)*m_springs_per_cell) ;beta_prolif_m((cell_event_EMT)*m_springs_per_cell + 1:end)];
            
            
            %Update the chemical concentrations - chemical 1 - cell level
            chem_1_m = [chem_1_m(1:(cell_event_EMT-1)*m_springs_per_cell) ;chem_1_m((cell_event_EMT)*m_springs_per_cell + 1:end)];
            
            %Update mixing_pop_id
            mixing_pop_id_m = [mixing_pop_id_m(1:(cell_event_EMT-1)*m_springs_per_cell) ;mixing_pop_id_m((cell_event_EMT)*m_springs_per_cell +1:end)];
            
            %Update the number of cells
            N_cells = N_cells - 1;
            
            N_springs = N_springs - m_springs_per_cell;
            
            cell_with_EMT_event_hist = [cell_with_EMT_event_hist,cell_event_EMT];
            
            prolif_EMT_indicator_hist = [prolif_EMT_indicator_hist, 3]; %2 correspond to EMT
            
            cell_event_hist = [cell_event_hist,cell_event_EMT];
            
            cell_event_times = [cell_event_times, t];
            
            EMT_times = [EMT_times, t];
            
            
        end
    end
    
    
    
end