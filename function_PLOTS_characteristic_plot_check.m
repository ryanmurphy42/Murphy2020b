function function_characteristic_plot_check(num_sims,t_hist_sim_run,x_hist_sim_run)

%Check to see if dt was chosen to small enough to reproduce characteristics
%If there is output from this then choose dt smaller.

for sim_run = 1:num_sims
    
	t_hist=t_hist_sim_run{sim_run};
    x_hist=x_hist_sim_run{sim_run};
    
    
    x_hist_double_change =[];
    for t_step = 2:size( t_hist,2)
        if t_step == 2
            size_prev_step = size(x_hist{t_step},1);
        else
            size_prev_step = size_now;
        end
        
        size_now = size(x_hist{t_step},1);
        
        if size_prev_step == size_now -1
        elseif size_prev_step == size_now + 1
        elseif size_prev_step == size_now
        else
            x_hist_double_change = [x_hist_double_change,t_hist(t_step) ];
            t_step
        end
    end
    
    if isempty(x_hist_double_change) == 0
        disp(num2str(sim_run))
    end
    
end


end
