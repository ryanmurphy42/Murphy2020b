function function_PLOTS_characteristics(t_hist, x_hist, k_hist, mixing_pop_id_hist, a_hist,...
    filepath_save_figs,  colouring, cell_with_proliferation_event_hist,cell_with_EMT_event_hist,...
    cell_event_hist,prolif_EMT_indicator_hist,    sim_run,chem_1_hist,chem_1_inc,timestop)

%function_characteristcs_plot - plot characteristics for each realisation sequentially

%t_hist - times recorded
%x_hist - positions
%k_hist - cell stiffness
%a_hist - resting spring length
% mixing_pop_id_hist - 1 if a cell, 0 if not used
%filepath_save_figs
%colouring %1- density, 2 - chemical

m_springs_per_cell=1;

%% check if x_hist changes by more than one in any timestep

x_hist_double_change =[];
for t_step = 2:size( t_hist,2)
    if t_step == 2
        size_prev_step = size(x_hist{t_step},1);
    else
        size_prev_step = size_now;
    end
    
    size_now = size(x_hist{t_step},1);
    
    if size_prev_step == size_now -m_springs_per_cell
    elseif size_prev_step == size_now + m_springs_per_cell
    elseif size_prev_step == size_now
    else
        x_hist_double_change = [x_hist_double_change,t_hist(t_step) ];
    end
end


%% Calculate the proliferation times

prolif_times =[];
for t_step = 2:size( t_hist,2)
    if t_step == 2
        size_prev_step = size(x_hist{t_step},1);
    else
        size_prev_step = size_now;
    end
    
    size_now = size(x_hist{t_step},1);
    
    if size_prev_step < size_now
        prolif_times = [prolif_times,t_hist(t_step) ];
    end
end


%% Calculate the cell event times

event_times =[];
t_hist_cell_change_index = [1]; %for characteristic surf plot
for t_step = 2:size(t_hist,2)
    if t_step == 2
        size_prev_step = size(x_hist{t_step},1);
    else
        size_prev_step = size_now;
    end
    
    size_now = size(x_hist{t_step},1);
    
    if size_prev_step > size_now %death event
        event_times = [event_times,t_hist(t_step) ];
        t_hist_cell_change_index = [t_hist_cell_change_index, t_step];
    elseif size_prev_step < size_now %prolif event
        event_times = [event_times,t_hist(t_step) ];
        t_hist_cell_change_index = [t_hist_cell_change_index, t_step];
    end
end
if t_hist_cell_change_index(end) ~= t_step
    t_hist_cell_change_index = [t_hist_cell_change_index,t_step];
end

%% Calculate traceback history of particles

ch_plot_ord = {};

t_step_ch=1;
cell_prolif_count=0;
cell_EMT_count = 0;
cell_death_index_count=0;

for t_step_ch = 1:(size(event_times,2)+1)
    
    if t_step_ch == 1 %if initial time label the cells
        ch_plot_ord{t_step_ch} = 1:size(x_hist{1},1);
    else %if not the initial time update the cell ordering
        
        %the cell ordering from the previous time step
        int_ch_plot_ord = ch_plot_ord{t_step_ch-1};
        
        %cell event
        cell_ev = cell_event_hist(t_step_ch-1);
        
        %determine whether prolif or death event
        if prolif_EMT_indicator_hist(t_step_ch-1) ==1 %if a proliferation event
            
            cell_prolif_count = cell_prolif_count + 1; %increase loop counter for cell prolif
            
            cell_prolif_ev = cell_with_proliferation_event_hist(cell_prolif_count);
            
            cell_prolif_ev_m = cell_prolif_ev*m_springs_per_cell;
            
            for int_ch_plot_ord_loop = 1:size(int_ch_plot_ord,2)
                
                if int_ch_plot_ord(int_ch_plot_ord_loop) <= cell_prolif_ev_m %if to the left index the same
                    int_ch_plot_ord(int_ch_plot_ord_loop) = int_ch_plot_ord(int_ch_plot_ord_loop);
                    
                elseif int_ch_plot_ord(int_ch_plot_ord_loop) > cell_prolif_ev_m  %if to the right index increase
                    int_ch_plot_ord(int_ch_plot_ord_loop) = int_ch_plot_ord(int_ch_plot_ord_loop) + 1;
                end
                
            end
            %add on the cell which proliferated
            int_ch_plot_ord = [int_ch_plot_ord,cell_prolif_ev+1];
            ch_plot_ord{t_step_ch} = int_ch_plot_ord;
            
            
        elseif prolif_EMT_indicator_hist(t_step_ch-1) == 3  %if a EMT event
            
            cell_EMT_count = cell_EMT_count + 1; %increase loop counter to cell death
            
            cell_EMT_ev = cell_with_EMT_event_hist(cell_EMT_count);
            
            for int_ch_plot_ord_loop = 1:size(int_ch_plot_ord,2)
                
                if int_ch_plot_ord(int_ch_plot_ord_loop) > 0
                    if int_ch_plot_ord(int_ch_plot_ord_loop) <= cell_EMT_ev %if to the left index the same
                        int_ch_plot_ord(int_ch_plot_ord_loop) = int_ch_plot_ord(int_ch_plot_ord_loop);
                        
                    elseif int_ch_plot_ord(int_ch_plot_ord_loop) > cell_EMT_ev+1  %if to the right index decrease
                        int_ch_plot_ord(int_ch_plot_ord_loop) = int_ch_plot_ord(int_ch_plot_ord_loop) -1 ;
                        
                    elseif int_ch_plot_ord(int_ch_plot_ord_loop) == cell_EMT_ev+1  %if to the right index decrease
                        cell_death_index_count = cell_death_index_count -1;
                        int_ch_plot_ord(int_ch_plot_ord_loop) = cell_death_index_count ;
                    end
                end
            end
            ch_plot_ord{t_step_ch} = int_ch_plot_ord;
        end
        
    end
end

%%
% traceback from the final state to determine the cells at each timestep
% determine the starting point for each characteristic

x_order_hist = {};
traceback_Nfin_loop=1;
prolif_counter=0;
for traceback_Nfin_loop = 1:size(ch_plot_ord{end},2)
    inter_hist = [];
    if traceback_Nfin_loop <= size(x_hist{1},1)
        for k_loop = 1:size(ch_plot_ord,2)
            inter_hist = [inter_hist,ch_plot_ord{k_loop}(traceback_Nfin_loop)];
        end
        x_order_hist{traceback_Nfin_loop} =   inter_hist;
    else
        %determine number of death events before cell boundary was generated then know the start time of this event
        
        %a new characteristic for each proliferation event
        %at what position did the proliferation event occur
        
        prolif_counter=prolif_counter+1;
        
        cum_sum_prolif_EMT_indicator_hist  =cumsum(prolif_EMT_indicator_hist==1);
        
        [~,time_begin_for_traj]= min(abs(cum_sum_prolif_EMT_indicator_hist-prolif_counter));
        
        for k_loop=1:(size(ch_plot_ord,2) -time_begin_for_traj)
            inter_hist = [inter_hist,ch_plot_ord{size(ch_plot_ord,2) - k_loop + 1}(traceback_Nfin_loop)];
        end
        x_order_hist{traceback_Nfin_loop} =   fliplr(inter_hist);
    end
end

%% determine the start times for each of the characteristic (based on the final positions)
%% and end times of each characteristics

%for each cell determine the positions
x_hist_traceback = {};
t_hist_traceback = {};
q_hist_traceback = {};
k_hist_traceback = {};
a_hist_traceback = {};
G_net_hist_traceback = {};
mixing_pop_id_hist_traceback = {};
chem_1_hist_traceback = {};

for traceback_Nfin_loop = 1:size(ch_plot_ord{end},2) %for each of the cells
    x_hist_follow_boundary =[];
    t_hist_follow_boundary =[];
    q_hist_follow_boundary =[];
    k_hist_follow_boundary =[];
    a_hist_follow_boundary =[];
    G_net_hist_follow_boundary = [];
    mixing_pop_id_hist_follow_boundary = [];
    chem_1_hist_follow_boundary =[];

    for t_step = 1:size(t_hist,2) %for each timestep
        
        t_val = t_hist(t_step);
        
        %determine where t_val sits within event_times
        event_times_inf = event_times - t_val;
        event_times_inf(event_times_inf <= 0) = inf;
        
        if sum(event_times_inf== inf) == size(event_times_inf,2)
            event_time_pos = size(event_times,2)+1;
        else
            [~, event_time_pos] = min(event_times_inf);
        end
        
        %does the cell exist at this time
        if event_time_pos >= size(event_times,2) - size(x_order_hist{traceback_Nfin_loop},2) + 2 %cell exists
            
            event_time_pos_with_prolif = event_time_pos +  (size(x_order_hist{traceback_Nfin_loop},2)-size(x_order_hist{1},2));
            
            ch_plot_order_x = x_order_hist{traceback_Nfin_loop}(event_time_pos_with_prolif);
            %ch_plot_ord{prolif_time_pos}(traceback_Nfin_loop);
            
            if ch_plot_order_x > 0
                
                x_val = x_hist{t_step}(ch_plot_order_x);
                x_hist_follow_boundary = [x_hist_follow_boundary,x_val];
                t_hist_follow_boundary = [t_hist_follow_boundary,t_val];
                
                %density
                if  x_hist{t_step}(ch_plot_order_x) ==  x_hist{t_step}(1) %if at the first position
                    q_val = (1/(x_hist{t_step}(ch_plot_order_x+1) - x_hist{t_step}(ch_plot_order_x)));
                    q_hist_follow_boundary = [q_hist_follow_boundary,q_val];
                    
                elseif x_hist{t_step}(ch_plot_order_x) ==  x_hist{t_step}(end) %at the last position
                    q_val = (1/(x_hist{t_step}(ch_plot_order_x) - x_hist{t_step}(ch_plot_order_x-1)));
                    q_hist_follow_boundary = [q_hist_follow_boundary,q_val];
                    
                else %interior cell
                    q_val = 0.5*( (1/(x_hist{t_step}(ch_plot_order_x+1) - x_hist{t_step}(ch_plot_order_x))) + (1/(x_hist{t_step}(ch_plot_order_x) - x_hist{t_step}(ch_plot_order_x-1))));
                    q_hist_follow_boundary = [q_hist_follow_boundary,q_val];
                end
                
                %cell stiffness
                if  x_hist{t_step}(ch_plot_order_x) ==  x_hist{t_step}(1) %if at the first position
                    k_val = k_hist{t_step}(ch_plot_order_x);
                    k_hist_follow_boundary = [k_hist_follow_boundary,k_val];
                    
                elseif x_hist{t_step}(ch_plot_order_x) ==  x_hist{t_step}(end) %at the last position
                    k_val = k_hist{t_step}(ch_plot_order_x-1);
                    k_hist_follow_boundary = [k_hist_follow_boundary,k_val];
                    
                else %interior cell
                    k_val = 0.5*( k_hist{t_step}(ch_plot_order_x-1) + k_hist{t_step}(ch_plot_order_x));
                    k_hist_follow_boundary = [k_hist_follow_boundary,k_val];
                end
                
                %resting cell length
                if  x_hist{t_step}(ch_plot_order_x) ==  x_hist{t_step}(1) %if at the first position
                    a_val = a_hist{t_step}(ch_plot_order_x);
                    a_hist_follow_boundary = [a_hist_follow_boundary,a_val];
                    
                elseif x_hist{t_step}(ch_plot_order_x) ==  x_hist{t_step}(end) %at the last position
                    a_val = a_hist{t_step}(ch_plot_order_x-1);
                    a_hist_follow_boundary = [a_hist_follow_boundary,a_val];
                    
                else %interior cell
                    a_val = 0.5*( a_hist{t_step}(ch_plot_order_x-1) + a_hist{t_step}(ch_plot_order_x));
                    a_hist_follow_boundary = [a_hist_follow_boundary,a_val];
                end
                
                %mixing populations id
                
                if iscell(mixing_pop_id_hist) == 0
                    mixing_pop_id_hist_follow_boundary=-1;
                else
                    if m_springs_per_cell == 1
                        if  x_hist{t_step}(ch_plot_order_x) ==  x_hist{t_step}(1) %if at the first position
                            mixing_pop_id_val = mixing_pop_id_hist{t_step}(ch_plot_order_x);
                            mixing_pop_id_hist_follow_boundary = [mixing_pop_id_hist_follow_boundary,mixing_pop_id_val];
                            
                        elseif x_hist{t_step}(ch_plot_order_x) ==  x_hist{t_step}(end) %at the last position
                            mixing_pop_id_val = mixing_pop_id_hist{t_step}(ch_plot_order_x-1);
                            mixing_pop_id_hist_follow_boundary = [mixing_pop_id_hist_follow_boundary,mixing_pop_id_val];
                            
                        else
                            mixing_pop_id_val = 0.5*( mixing_pop_id_hist{t_step}(ch_plot_order_x-1) + mixing_pop_id_hist{t_step}(ch_plot_order_x) );
                            mixing_pop_id_hist_follow_boundary = [mixing_pop_id_hist_follow_boundary,mixing_pop_id_val];
                        end
                    else
                        if  x_hist{t_step}(ch_plot_order_x) ==  x_hist{t_step}(1) %if at the first position
                            mixing_pop_id_val = mixing_pop_id_hist{t_step}(ch_plot_order_x);
                            mixing_pop_id_hist_follow_boundary = [mixing_pop_id_hist_follow_boundary,mixing_pop_id_val];
                            
                        elseif x_hist{t_step}(ch_plot_order_x) ==  x_hist{t_step}(end) %at the last position
                            mixing_pop_id_val = mixing_pop_id_hist{t_step}(ch_plot_order_x-1);
                            mixing_pop_id_hist_follow_boundary = [mixing_pop_id_hist_follow_boundary,mixing_pop_id_val];
                            
                        else
                            mixing_pop_id_val = 0.5*( mixing_pop_id_hist{t_step}(ch_plot_order_x-1) + mixing_pop_id_hist{t_step}(ch_plot_order_x) );
                            mixing_pop_id_hist_follow_boundary = [mixing_pop_id_hist_follow_boundary,mixing_pop_id_val];
                        end
                    end
                end
                
                if chem_1_inc == 1
                    % chemical 1 colouring
                    
                    if  x_hist{t_step}(ch_plot_order_x) ==  x_hist{t_step}(1) %if at the first position
                        chem_1_val = chem_1_hist{t_step}(ch_plot_order_x);
                        chem_1_hist_follow_boundary = [chem_1_hist_follow_boundary,chem_1_val];
                        
                    elseif x_hist{t_step}(ch_plot_order_x) ==  x_hist{t_step}(end) %at the last position
                        chem_1_val = chem_1_hist{t_step}(ch_plot_order_x-1);
                        chem_1_hist_follow_boundary = [chem_1_hist_follow_boundary,chem_1_val];
                        
                    else %interior cell
                        chem_1_val = 0.5*( chem_1_hist{t_step}(ch_plot_order_x-1) + chem_1_hist{t_step}(ch_plot_order_x));
                        chem_1_hist_follow_boundary = [chem_1_hist_follow_boundary,chem_1_val];
                    end
                    
                end
                
                
            end
        end
    end
    
    x_hist_traceback{traceback_Nfin_loop} = x_hist_follow_boundary;
    t_hist_traceback{traceback_Nfin_loop} = t_hist_follow_boundary;
    q_hist_traceback{traceback_Nfin_loop} = q_hist_follow_boundary;
    k_hist_traceback{traceback_Nfin_loop} = k_hist_follow_boundary;
    a_hist_traceback{traceback_Nfin_loop} = a_hist_follow_boundary;
    mixing_pop_id_hist_traceback{traceback_Nfin_loop} = mixing_pop_id_hist_follow_boundary;
    
    if chem_1_inc == 1
        chem_1_hist_traceback{traceback_Nfin_loop} = chem_1_hist_follow_boundary;
    end
end


%% Calculate the value inside the cell to fill the inside of the characteristics

if colouring == 1
    
    
    %calculate time points when cell events occur
    traceback_Nfin_loop_vec = 1:size(ch_plot_ord{end},2);
    traceback_Nfin_loop_vec = [traceback_Nfin_loop_vec(1:length(x_hist{1})-1),traceback_Nfin_loop_vec((length(x_hist{1})+1):end)];
    
    figure
    for traceback_Nfin_loop = traceback_Nfin_loop_vec
        x_plot =  x_hist_traceback{traceback_Nfin_loop};
        y_plot =  t_hist_traceback{traceback_Nfin_loop};
        plot3(x_plot, y_plot, 1000*ones(size(x_plot)), ...
            'k', 'LineWidth', 0.1)
        hold on
    end
    
    timesteps_with_same_cell=[];
    if t_hist_cell_change_index(end) ~= length(t_hist)
        timesteps_with_same_cell = [ diff(t_hist_cell_change_index),length(t_hist)];
    else
        timesteps_with_same_cell =  diff(t_hist_cell_change_index);
    end
    
    for ii=1:length(t_hist_cell_change_index)-1
        
        x_hist_ii_more = zeros(length(x_hist{t_hist_cell_change_index(ii)}),timesteps_with_same_cell(ii));
        
        for jj1=1:timesteps_with_same_cell(ii)
            x_hist_jj1 = x_hist{t_hist_cell_change_index(ii)-1+jj1};
            for jj2=1:length(x_hist{t_hist_cell_change_index(ii)})
                x_hist_ii_more(jj2,jj1) = x_hist_jj1(jj2);
            end
        end
        jj1=timesteps_with_same_cell(ii)+1;
        x_hist_jj1 = x_hist{t_hist_cell_change_index(ii)-1+jj1-1};
        for jj2=1:length(x_hist{t_hist_cell_change_index(ii)})
            x_hist_ii_more(jj2,jj1) = x_hist_jj1(jj2);
        end
        
        
        den_hist_ii_more =zeros(length(x_hist{t_hist_cell_change_index(ii)})-1,timesteps_with_same_cell(ii)+1);
        for jj1=1:timesteps_with_same_cell(ii)+1
            if jj1 < timesteps_with_same_cell(ii)+1
                
                x_hist_jj1 = x_hist{t_hist_cell_change_index(ii)-1+jj1};
            end
            for jj2=1:(length(x_hist{t_hist_cell_change_index(ii)})-1)
                den_hist_ii_more(jj2,jj1) = 1/(x_hist_jj1(jj2+1)-x_hist_jj1(jj2));
            end
            
        end
        
        
        for jj = 1:(length(x_hist{t_hist_cell_change_index(ii)})-2)
            surf(x_hist_ii_more(jj:jj+1,:), ...
                [t_hist(t_hist_cell_change_index(ii):(t_hist_cell_change_index(ii+1)));t_hist(t_hist_cell_change_index(ii):(t_hist_cell_change_index(ii+1)))]...
                ,[den_hist_ii_more(jj:jj,:);den_hist_ii_more(jj:jj,:)], ...
                'EdgeColor', 'none')
            hold on
        end
        
    end
    view(2)
    colorbar
    shg
    
    title(['Characteristics - Density Colouring - Sim Run ' num2str(sim_run)])
    
elseif colouring == 2
    
    traceback_Nfin_loop_vec = 1:size(ch_plot_ord{end},2);
    traceback_Nfin_loop_vec = [traceback_Nfin_loop_vec(1:length(x_hist{1})-1),traceback_Nfin_loop_vec((length(x_hist{1})+1):end)];
    
    figure
    for traceback_Nfin_loop = traceback_Nfin_loop_vec
        x_plot =  x_hist_traceback{traceback_Nfin_loop};
        y_plot =  t_hist_traceback{traceback_Nfin_loop};
        plot3(x_plot, y_plot, 10000000*ones(size(x_plot)), ...
            'k', 'LineWidth', 1)
        hold on
    end
    
    % t_hist_cell_change_index
    timesteps_with_same_cell=[];
    if t_hist_cell_change_index(end) ~= length(t_hist)
        timesteps_with_same_cell = [ diff(t_hist_cell_change_index),length(t_hist)];
    else
        timesteps_with_same_cell =  diff(t_hist_cell_change_index);
    end
    
    for ii=1:length(t_hist_cell_change_index)-1
        
        x_hist_ii_more = zeros(length(x_hist{t_hist_cell_change_index(ii)}),timesteps_with_same_cell(ii));
        
        for jj1=1:timesteps_with_same_cell(ii)
            x_hist_jj1 = x_hist{t_hist_cell_change_index(ii)-1+jj1};
            for jj2=1:length(x_hist{t_hist_cell_change_index(ii)})
                x_hist_ii_more(jj2,jj1) = x_hist_jj1(jj2);
            end
        end
        jj1=timesteps_with_same_cell(ii)+1;
        x_hist_jj1 = x_hist{t_hist_cell_change_index(ii)-1+jj1-1};
        for jj2=1:length(x_hist{t_hist_cell_change_index(ii)})
            x_hist_ii_more(jj2,jj1) = x_hist_jj1(jj2);
        end
        
        
        chem_1_hist_ii_more =zeros(length(x_hist{t_hist_cell_change_index(ii)})-1,timesteps_with_same_cell(ii)+1);
        for jj1=1:timesteps_with_same_cell(ii)+1
            if jj1 < timesteps_with_same_cell(ii)+1
                
                chem_1_hist_jj1 = chem_1_hist{t_hist_cell_change_index(ii)-1+jj1};
            end
            for jj2=1:(length(x_hist{t_hist_cell_change_index(ii)})-1)
                chem_1_hist_ii_more(jj2,jj1) = chem_1_hist_jj1(jj2);
            end
            
        end
        
        
        for jj = 1:(length(x_hist{t_hist_cell_change_index(ii)})-2)
            surf(x_hist_ii_more(jj:jj+1,:), ...
                [t_hist(t_hist_cell_change_index(ii):(t_hist_cell_change_index(ii+1)));t_hist(t_hist_cell_change_index(ii):(t_hist_cell_change_index(ii+1)))]...
                ,[chem_1_hist_ii_more(jj:jj,:);chem_1_hist_ii_more(jj:jj,:)], ...
                'EdgeColor', 'none')
            hold on
        end
        
    end
    view(2)
    colorbar
    colormap(flipud(hot))
    shg
    
    title(['Characteristics - Chemical Colouring - Sim Run ' num2str(sim_run)])
    
end

%flip the y axis axis
ax = gca;
ax.YDir = 'reverse';
ylim([0,timestop])

box on

%% Save the plot

print(gcf,'-depsc2',  [filepath_save_figs 'Characteristicsfill_with_colourbar_type_' num2str(colouring) '_sim_run_' num2str(sim_run) '2.eps']);
saveas(gcf, [filepath_save_figs 'Characteristicsfill_with_colourbar_type_' num2str(colouring) '_sim_run_' num2str(sim_run) '.fig'])
saveas(gcf, [filepath_save_figs 'Characteristicsfill_with_colourbar_type_' num2str(colouring) '_sim_run_' num2str(sim_run) '.jpg'])


end
