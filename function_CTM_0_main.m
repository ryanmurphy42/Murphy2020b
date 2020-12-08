function function_CTM_0_main(simulation_id,N_cells_init,  beta_1,k1,a1,eta,omega,...
    prolif_law,phi,D_chem_1,timestop,q_init_condition,EMT_law,dt_ctm,...
    time_record_vector_ctm,chem_1_inc, chem_1_init_condition,dz_ctm,err_tol_ctm)
	
	%% Main script to run continuum model

tic

%% Numerical simulation parameters

dz=dz_ctm;
dt_ori = dt_ctm;
err_tol = err_tol_ctm;
max_iters=10;
nodesz=round(1/dz) +1;

%% cell properties initial conditions

k=k1*ones(nodesz,1);
a=a1*ones(nodesz,1);
beta_prolif=beta_1*ones(nodesz,1);

a_cells=[a1*ones(N_cells_init,1)];

%% q initial condition

if q_init_condition == 1
    %Uniform initial condition
    q0 = zeros(nodesz,1);
    q1=1/a1;
    L=sum(a_cells);
    for i=1:nodesz
        q0(i) = q1;
    end
    
elseif q_init_condition == 2
    %Compressed to half of initial resting cell length
    q0 = zeros(nodesz,1);
    q1=1/(0.5*a1);
    L=0.5*sum(a_cells);
    for i=1:nodesz
        q0(i) = q1;
    end
    
elseif q_init_condition == 3
    %Stretched to double resting cell length
    q0 = zeros(nodesz,1);
    q1=1/(2*a1);
    L=2*sum(a_cells);
    for i=1:nodesz
        q0(i) = q1;
    end
end



if chem_1_init_condition == 1
    % no chemical anywhere
    chem_1=zeros(nodesz,1);
end

%% create a folder to store data

file_save_name = ['Results_Sim' simulation_id '_Continuum'];

folder_name = [simulation_id '_Continuum'];
if ~exist([folder_name'], 'dir')
    mkdir([folder_name]);
end


%% convert time_record_vector_ctm into timesteps

time_record_vector_dt = round(time_record_vector_ctm/dt_ctm,0);

len_time_record_vector_ctm = length(time_record_vector_ctm);

%% storage and save initial values

loop_count_stored = 1;
max_recorded_data_points = length(time_record_vector_dt); %value to initialise for below

q_hist = zeros(nodesz,max_recorded_data_points);
q_hist(:,loop_count_stored) = q0;
q=q0;

k_hist = zeros(nodesz, max_recorded_data_points);
k_hist(:,loop_count_stored) = k;

a_hist = zeros(nodesz, max_recorded_data_points);
a_hist(:,loop_count_stored) = a;

beta_prolif_hist = zeros(nodesz, max_recorded_data_points);
beta_prolif_hist(:,loop_count_stored) = beta_prolif;

node_velocity_hist=zeros(nodesz,max_recorded_data_points);

t_hist  = zeros(1,max_recorded_data_points);
t_hist(1)=0;

L_hist   = zeros(1,max_recorded_data_points);
L_hist(1)=L;

chem_1_hist = zeros(nodesz,max_recorded_data_points);
chem_1_hist(:,loop_count_stored) = chem_1;

%%

t=0; %initial time

q_update = q;
L_update = L;
k_update = k;
a_update = a;
beta_prolif_update = beta_prolif ;
chem_1_update = chem_1 ;

time_record_vector_ctm_loop = 1;


while t < timestop
    
    
    
    dt = dt_ori; % reset time step
    
	% update the values for the previous time step
    q_old = q; 
    L_old = L;
    k_old = k;
    a_old = a;
    beta_prolif_old =beta_prolif;
    chem_1_old = chem_1;
    
	
    res = err_tol + 1; % set the residual error bigger than the error tolerance so simulation starts
     
    w=1; % set the initial Newton-Raphson iterate.
    
    dt_level=1;
    
    
    if L_old > 0.025 % Run until close to extinction
        %Newton iteration
        while w <= max_iters && res > err_tol
            
            if w > (max_iters - 2)
                w=1;
                dt=dt/10;
                dt_level = dt_level + 1;
            end
            
            q = q_update;
            L = L_update;
            k = k_update;
            a = a_update;
            beta_prolif = beta_prolif_update;
            chem_1 = chem_1_update;
            
            %% Solve for l
            if EMT_law== 2
                % use chemical to update EMT1
                
                if chem_1_update(end) >= 500
                    EMT1 = omega/(1-phi);
                else
                    EMT1 = 0;
                end
                
                [L_update] = function_CTM_boundary_position(L_old,q_update,k_update,a_update,dt,eta,nodesz, EMT1);
                
                dL=L_update-L_old;
                
            elseif EMT_law== 1
                
                EMT1= omega;
                [L_update] = function_CTM_boundary_position(L_old,q_update,k_update,a_update,dt,eta,nodesz, EMT1);
            end
            
            
            
            %% Solve for density, q
            [L_diag,D_diag,U_diag] = function_CTM_jacobian_q(q,k,nodesz,dt,dz,eta,a,L_old,L);
            rhs = -function_CTM_discretised_func_q(q,q_old,k,a,nodesz,dt,dz,eta,L_old,L,prolif_law,beta_prolif);
            dq = function_CTM_tridia(nodesz,L_diag,D_diag,U_diag,rhs); % solve solution using thomas algorthim
            q_update = q + dq;
            
            %% Solve for chem_1
            
            if chem_1_inc == 1
                
                % calculate node velocity
                [node_velocity] = function_CTM_node_velocity(q,k,a,eta,nodesz,dz,L_old,L,dt);
                
                % update the chem_1
                [L_diag,D_diag,U_diag] = function_CTM_jacobian_chem_1(q,k,nodesz,dt,dz,eta,a,L_old,L,chem_1, chem_1_old,D_chem_1, node_velocity);
                rhs = -function_CTM_discretised_func_chem_1(q,k,a,nodesz,dt,dz,eta,L_old,L,prolif_law,beta_prolif,chem_1, chem_1_old,D_chem_1,node_velocity,omega, phi);
                dchem_1 = function_CTM_tridia(nodesz,L_diag,D_diag,U_diag,rhs); % solve solution using thomas algorthim
                chem_1_update = chem_1 + dchem_1;
                
            else
                dchem_1=0;
            end
            
            
            %% Calculate the residual
            
            res = norm([dq;dchem_1],inf); % update residual
            
            w=w+1;
            
            
            
        end
        
        
        total_cells_update = trapz((0:dz:1)*L_update,q_update');
        
        if w >= max_iters
            fprintf('\n MAX ITERS REACHED! :t = %f; iters = %f; res = %f; res divide errtol = %f\n', t,w,res, res/err_tol*100);
        end
        
    end
    
    t = t + dt;
    
    fprintf('\n t = %f; dt_level = %f; iters = %f ; mass = %f; length = %f\n', t,dt_level, w,total_cells_update,L_update);
    
    
    %% Save if near timestep
    
    
    if time_record_vector_ctm_loop <= len_time_record_vector_ctm
        
        if abs(time_record_vector_ctm(time_record_vector_ctm_loop)-t) < 2*dt_ori
            time_record_vector_ctm_loop = time_record_vector_ctm_loop + 1;
            
            loop_count_stored = loop_count_stored + 1;
            q_hist(:,loop_count_stored) = q_update;
            k_hist(:,loop_count_stored) = k_update;
            a_hist(:,loop_count_stored) = a_update;
            beta_prolif_hist(:,loop_count_stored) = beta_prolif_update;
            t_hist(loop_count_stored) = t;
            L_hist(loop_count_stored) = L_update;
            chem_1_hist(:,loop_count_stored) = chem_1_update;
            
            fprintf('\n t = %f; iters = %f ; mass = %f; length = %f\n', t,w,total_cells_update,L_update);
            
        end
        
    end
    if sum(isnan(q_update)) > 1
        fprintf('\n chem_1_update is NAN t = %f; iters = %f ; mass = %f; length = %f\n', t,w,total_cells_update,L_update);
        break
    end
    
    
    if sum(isnan(chem_1_update)) > 1
        fprintf('\n chem_1_update is NAN t = %f; iters = %f ; mass = %f; length = %f\n', t,w,total_cells_update,L_update);
        break
    end
    
end

toc

%% Save the full .mat file.

save([pwd '/' folder_name '/' file_save_name],'-v7.3');

end