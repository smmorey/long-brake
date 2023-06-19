'''
Charlie Shobe & Susannah Morey
June 19, 2023

longBRaKE: Longitudinal Profile Blocky River And Knickpoint Evolution

Starting from the 1-D reach model published in GRL 2016

Expanded the 1-D reach model into a full 1-D long profile,
added depositional block event capabilities.

'''
#Packages:
from __future__ import division
import numpy as np
import sys
import os
np.set_printoptions(threshold=np.inf)

class Brake(object):

    def __init__(self):
        pass
    
    def import_data(self, suffix):
        self.time = np.load('time_record_' + suffix + '.npy')
        self.elev = np.load('elev_record_' + suffix + '.npy')
        self.slope = np.load('slope_record_' + suffix + '.npy')
        #self.block = np.load('block_count_record_' + suffix + '.npy')
        #self.cover = np.load('cover_frac_record_' + suffix + '.npy')
        
    def save_to_txt(self, suffix):
        np.savetxt('time_record_' + suffix + '.txt', self.time)
        np.savetxt('elev_record_' + suffix + '.txt', self.elev)
        np.savetxt('slope_record_' + suffix + '.txt', self.slope)
        np.savetxt('block_count_record_' + suffix + '.txt', self.block)
        np.savetxt('cover_frac_record_' + suffix + '.txt', self.cover)
            
#    def model_response_data(self, domain_length):
#        self.integrated_elev_lost = sum(self.elev[0, :] - self.elev[-1, :])
#        self.domain_length = domain_length
        
    def save_model_params(self, gamma, bl_drop, delay_timescale, d, timestep, run_time, suffix, total_volume_lost, n_cells, dx, bed_k, hillslope_block_k, fluvial_block_k):
        params_list = ['gamma', 'bl_drop', 'delay_timescale', 'block_side_length',
                       'timestep', 'run_time', 'bed_k', 'hillslope_block_k', 'fluvial_block_k', 'total_volume_lost', 'domain_length', 'total_blocks']
        params = {'gamma': gamma, 'bl_drop': bl_drop, 'delay_timescale': delay_timescale,
                  'block_side_length': d, 'timestep': timestep, 'run_time': run_time, 'bed_k': bed_k,
                  'hillslope_block_k': hillslope_block_k, 'fluvial_block_k': fluvial_block_k, 'total_volume_lost': total_volume_lost,
                  'domain_length': n_cells * dx, 'total_blocks': self.blocks}
        os.system("touch model_params" + suffix + ".txt")
        for i in range(len(params)):
            line = str(params_list[i]) + ' = ' + str(params[params_list[i]]) + '\n'
            with open("model_params" + suffix + ".txt", "a") as f: 
                f.write(line)  
    
#    def import_discharges(self, file_name, num_iterations):
#        q_distribution_initial = np.genfromtxt(file_name, delimiter = ',')
#        q_distribution = q_distribution_initial
#        while len(q_distribution) <= num_iterations:
#            q_distribution = np.append(q_distribution, q_distribution_initial)
#        return q_distribution
    
    def calc_shear_stress_with_roughness(self, pre_adjusted_tau, tracking_mat, drag_cube, flow_depth, dx, roughness_height, slicing_index, cell):
        #blocks_above_flow = (self.is_block_in_cell) #& (tracking_mat[0:slicing_index, 2] >= flow_depth)        
        a1 = 6.5 #ferguson 2007
        a2 = 2.5 #ferguson 2007        
        if np.count_nonzero(self.is_block_in_cell) == 0:
            sigma_d_blocks = 0
        else:
            #blocks_above_flow = (self.is_block_in_cell) #& (tracking_mat[0:slicing_index, 2] >= flow_depth)
            beta = (a1 * (flow_depth / roughness_height)) / np.power(np.power(flow_depth / roughness_height, 5 / 3) + np.power(a1 / a2, 2), 1/2)
            avg_diam_blocks = np.average(tracking_mat[0:slicing_index, 1][self.is_block_in_cell]) #!
            submerged_block = (self.is_block_in_cell) & (tracking_mat[0:slicing_index, 1] < flow_depth)#!
            emergent_block = (self.is_block_in_cell) & (tracking_mat[0:slicing_index, 1] >= flow_depth)#!
            tracking_mat[0:slicing_index, 3][submerged_block] = tracking_mat[0:slicing_index, 1][submerged_block]#!
            tracking_mat[0:slicing_index, 3][emergent_block] = flow_depth           #! 
            avg_submerged_height_blocks = np.average(tracking_mat[0:slicing_index, 3][self.is_block_in_cell]) #!
            #avg_spacing_blocks =  dx / np.count_nonzero(blocks_above_flow)
            spacing_squared = (self.dx * self.dy[cell]) / np.count_nonzero(self.is_block_in_cell)
            sigma_d_blocks = (1 / 2) * drag_cube * np.power(beta, 2) *(avg_submerged_height_blocks * avg_diam_blocks / spacing_squared)
        self.sigma_d_array[cell] = sigma_d_blocks
        adjusted_shear_stress = pre_adjusted_tau / (1 + sigma_d_blocks)
        return adjusted_shear_stress
        
    def roughness_bisection(self, q, s, d, g):
        if s <= 0:
            os.system("rm results.out")
            os.system("touch results.out")
            with open("results.out", "a") as f: 
                f.write('FAIL\n') 
            sys.exit("NEGATIVE SLOPE-- KILLING MODEL")
        else:
            pass
        a1 = 6.5
        a2 = 2.5
        coef_1 = (g * s * np.power(a1, 2)) / (np.power(q, 2) * np.power(d, 2))
        coef_2 = 0.0
        coef_3 = -1 / np.power(d, 5 / 3)
        coef_4 = -np.power(a1 / a2, 2)
        coef_array = np.array([coef_1, coef_2, coef_3, coef_4])
        #print coef_array
        roots = np.roots(coef_array)
        is_real = np.isreal(roots)
        root = np.real(roots[is_real])[0]
        #print root        
        h = np.power(root, 3 / 5)
        if np.isnan(np.sum(h)):
            print(h)
            sys.exit("NAN FOUND IN ROOT'S H-- KILLING MODEL")
        #rough_iter = 1
        #error = 1
        #while error > tol:
        #    rough_iter += 1
        #    if rough_iter > 100000:
        #        os.system("rm results.out")
         #       os.system("touch results.out")
         #       with open("results.out", "a") as f: 
         #           f.write('FAIL\n')
         #       print 'real_q = ', q
         #       print 'slope = ', s
         #       print 'd = ', d
         #       print 'h_up = ', h_up
         #       sys.exit("ROUGHNESS CALCULATION STUCK-- KILLING MODEL")
         #   else:
         #       pass
         #   h_mid = (h_low + h_up)/2
         #   q_mid = h_mid * np.sqrt(g * h_mid * s) * ((a1 * h_mid / d)) / np.sqrt(np.power(h_mid / d, 5/3) + np.power(a1 / a2, 2))
         #   error = abs(q_mid - q)
         #   if q_mid > q and error > tol:
         #       h_up = h_mid
         #   elif q_mid < q and error > tol:
         #       h_low = h_mid
         #   else:
         #       pass
        return h
    
    def convert_years_to_seconds(self, value):
        value_years = value * 365 * 24 * 3600
        return value_years
        
    def calculate_timing_and_recording_vars(self, record_time_int, time_to_run_yrs, timestep_yrs, imposed_rock_uplift_rate, delay_timescale):
        self.record_time_seconds = self.convert_years_to_seconds(record_time_int)#self.record_time_int * 365 * 24 * 3600
        #spinup_time = self.convert_years_to_seconds(spinup_time_yrs)#spinup_time_yrs * 365 * 24 * 3600
        #perturb_time = self.convert_years_to_seconds(time_to_run_yrs)#time_to_run_yrs * 365 * 24 * 3600
        #spindown_time = self.convert_years_to_seconds(spindown_time_yrs)#spindown_time_yrs * 365 * 24 * 3600
        self.timestep_yrs = timestep_yrs
        #self.time_to_run = spinup_time + perturb_time + spindown_time#945900000000 #units of seconds
        self.time_to_run = self.convert_years_to_seconds(time_to_run_yrs)
        self.timestep = self.convert_years_to_seconds(timestep_yrs)#timestep_yrs * 365 * 24 * 3600#30000000 #units of seconds.        
        #number_iterations_spinup = spinup_time / self.timestep
        #number_iterations_perturb = perturb_time / self.timestep
        #number_iterations_spindown = spindown_time / self.timestep
        #self.number_records = (spinup_time_yrs + time_to_run_yrs + spindown_time_yrs) / self.record_time_int
        self.number_records = time_to_run_yrs / record_time_int
        print('Saving ', str(int(self.number_records) + 1), ' records (including one for the initial condition)')
        #self.number_iterations = number_iterations_spinup + number_iterations_perturb + number_iterations_spindown #number of model runs
        self.number_iterations = self.time_to_run / self.timestep
        #self.baselevel_drop_annual = np.concatenate((np.repeat(bl_drop_spinup, number_iterations_spinup), np.repeat(bl_drop_perturb, number_iterations_perturb), np.repeat(bl_drop_spindown, number_iterations_spindown))) / 365 / 24 / 3600 #now in m/s
        self.rock_uplift_rate = imposed_rock_uplift_rate / 365 / 24 / 3600 #now in m/s
        #mem_time = delay_timescale #half_width / bl_drop_perturb #becomes units of years
        mem_timesteps = int(np.round(delay_timescale / (self.timestep / 3600 / 24 / 365))) #whole number of years
        self.incision_memory = np.zeros((mem_timesteps, self.n_cells))
        
    def instantiate_common_constants(self):  
        np.random.seed(50) #TURN THIS LINE ON TO ELIMINATE INHERENT VARIABILITY
        self.dens_water = 1000 #kg/m^3
        self.dens_sediment = 2650 #kg/m^3
        self.g = 9.81 #acceleration due to gravity (m/s^2)
        self.drag_cube = .8 #from Carling 1998
        self.coeff_static_fric = .6 #from Carling 1998
        
    def set_shear_stress_params(self, bed_k, tau_c):
        self.ke_br_bed = bed_k #bedrock is not easy to erode obviously
        #self.ke_br_block = block_k
        self.a_br = 1 #bedrock
        self.tau_c_br = tau_c #pascals. bedrock detachment critical shear stress
        
    def set_roughness_depth_calc_params(self, z_0):
        self.tolerance = .001 #1 mm
        self.z0 = z_0 #roughness height NEED TO MAKE THIS DEPENDENT ON SED SIZE DISTRIBUTION
        self.upper_test_lim = 10000 #upper q limit to test e.g. maximum possible flow height
        self.lower_test_lim = 0 #lowest possible flow height
        
    def define_cell_centers_and_edges(self, n_cells, dx, initial_slope, starting_height_adjustment):        
        cell_indices = np.arange(n_cells, dtype = np.float64)
        self.x_array = np.arange(0, n_cells, 1)
        self.cell_centers = (cell_indices * dx) + (dx / 2)
        self.upstream_edges_of_cell = self.cell_centers - (dx / 2)
        self.downstream_edges_of_cell = self.cell_centers + (dx / 2)
        self.surface_elev = -initial_slope * cell_indices * dx + starting_height_adjustment #starting surface elevation (m)
        #self.surface_elev = np.load('long_prof_3000m_equil.npy') # after run to equilibrium, can upload a longtiduinal profile here

    def instantiate_tracking_matrix(self, number_of_blocks, gamma, starting_dist, fluvial_block_k):
        number_of_blocks = int(number_of_blocks)
        #boulders = int(boulders)
        #cubes = int(cubes)
        
        #explanation of tracking matrix columns:
        #0: distance along model domain where block is located [m]
        #1: block side length (blocks are assumed to be cubes) [m]
        #2: block volume, or side length cubed [m^3]
        #3: used for tracking whether blocks are submerged or not in flow calcs
        #4: block erodibility (k)
        
        if number_of_blocks > 0:
            if gamma == 0:
                self.tracking_mat = np.zeros((number_of_blocks, 5), dtype = np.float64)#!
            else:
                self.tracking_mat = np.zeros((number_of_blocks * 1000, 5), dtype = np.float64) #able to track 6 attributes #!
        else:
        	#print('nothing in tracking mat')
        	self.tracking_mat = np.zeros((100000, 5), dtype=np.float64) #!
        self.tracking_mat[:, :] = np.nan#-9999. #THIS IS SO THAT I CAN DISTINGUISH UN-ENTERED VALUES FROM ACTUAL ZEROS
        #self.tracking_mat[boulders:boulders + cubes + 1, 0] = 0 #this column is worthless, why is it here? #!
        self.tracking_mat[0:number_of_blocks, 0] = starting_dist #starting distance #!
        self.tracking_mat[0:number_of_blocks + 1, 1] = self.cube_side_length #!
        self.tracking_mat[0:number_of_blocks + 1, 2] = np.power(self.tracking_mat[0:number_of_blocks + 1, 1], 3) #!
        self.tracking_mat[0:number_of_blocks + 1, 3] = 0#np.nan #TO BE USED FOR SUBMERGENCE/EMERGENCE #!
        self.tracking_mat[0:number_of_blocks + 1, 4] = fluvial_block_k
        self.cover_frac_array = np.zeros(self.n_cells, dtype = np.float64)
        #self.tracking_mat[0:number_of_pieces, 0] = np.random.random_integers(0, 99 * self.dx, sum(~np.isnan(self.tracking_mat[:, 0]))) #!
        self.tracking_mat[0:number_of_blocks, 0] = np.random.random_integers(0, 99 * self.dx, number_of_blocks) #! SMM: replace this line with whatever I want; I could import an array that is scaled to the model domain
        #print(self.tracking_mat)
        
    def instantiate_single_dt_arrays(self, n_cells, initial_slope, surface_elevation):
        self.slope = np.zeros(n_cells, dtype = np.float64)
        self.flow_depth = np.zeros(n_cells, dtype = np.float64)
        self.flow_velocity = np.zeros(n_cells, dtype = np.float64)
        self.surface_elev_array = np.zeros(n_cells, dtype = np.float64)
        self.old_surf_elev_array = np.zeros(n_cells, dtype = np.float64)   
        self.sigma_d_array = np.zeros(n_cells, dtype = np.float64)   
        self.uncorrected_tau_array = np.zeros(n_cells, dtype = np.float64)
        self.corrected_tau_array = np.zeros(n_cells, dtype = np.float64)        
        self.time_avg_inc_rate_array = np.zeros(n_cells, dtype = np.float64)
        #self.blocks_in_cells = np.zeros(n_cells, dtype=np.float64)
        #self.zero_to_one_temp = np.zeros(n_cells, dtype=np.float64)
        #self.one_to_two_temp = np.zeros(n_cells, dtype=np.float64)
        #self.two_to_three_temp = np.zeros(n_cells, dtype=np.float64)
        #self.three_to_four_temp = np.zeros(n_cells, dtype=np.float64)
        #self.four_to_five_temp = np.zeros(n_cells, dtype=np.float64)
        #self.five_to_six_temp = np.zeros(n_cells, dtype=np.float64)
        self.record_count = 0
        self.slope[:] = initial_slope
        self.old_surf_elev_array[:] = surface_elevation
        self.surface_elev_array[:] = surface_elevation
        
    def instantiate_record_keeping_arrays(self, time_to_run, record_time_seconds, number_records, n_cells, surface_elev_array, slope):
        number_records = int(number_records)
        self.time_record = np.arange(0, time_to_run + record_time_seconds, record_time_seconds) #records time
        self.surface_elev_record = np.zeros((number_records + 1, n_cells), dtype = np.float64) #records surface elevation
        self.slope_record = np.zeros((number_records + 1, n_cells), dtype = np.float64) #records slope
        self.block_count_record = np.zeros((number_records + 1, n_cells), dtype = np.float64) #records number of blocks in each cell
        self.mean_block_size_record = np.zeros((number_records + 1, n_cells), dtype = np.float64) #records mean size of blocks in each cell
        self.cover_frac_record = np.zeros((number_records + 1, n_cells), dtype = np.float64) #records cover fraction
        self.time_avg_inc_rate_record = np.zeros((number_records + 1, n_cells), dtype = np.float64)        
        self.uncorrected_tau_record = np.zeros((number_records + 1, n_cells), dtype = np.float64)
        self.corrected_tau_record = np.zeros((number_records + 1, n_cells), dtype = np.float64)
        #arrays for grain size
        #self.zero_to_one = np.zeros((number_records + 1, n_cells), dtype = np.float64)       
        #self.one_to_two = np.zeros((number_records + 1, n_cells), dtype = np.float64)
        #self.two_to_three = np.zeros((number_records + 1, n_cells), dtype = np.float64)
        #self.three_to_four = np.zeros((number_records + 1, n_cells), dtype = np.float64)
        #self.four_to_five = np.zeros((number_records + 1, n_cells), dtype = np.float64)
        #self.five_to_six = np.zeros((number_records + 1, n_cells), dtype = np.float64)
        #initial values
        self.time_record[0] = 0
        self.surface_elev_record[0, :] = surface_elev_array
        self.slope_record[0, :] = slope
        
    def instantiate_counting_parameters(self):
        self.time = 0
        self.run_count = 0
        self.blocks = 0
        
    def calculate_initial_array_values(self, number_of_blocks, x_array, dx, dy, side_length):
        self.for_slicing = int(number_of_blocks + 1)
        for x in x_array:
            self.is_block_in_cell = (self.tracking_mat[0:self.for_slicing, 0] >= self.upstream_edges_of_cell[x]) & (self.tracking_mat[0:self.for_slicing, 0] < self.downstream_edges_of_cell[x]) #!
            #arr_padding = np.repeat(False, len(self.tracking_mat[:, 0]) - self.for_slicing)  
            #self.is_block_in_cell = np.append(self.is_block_in_cell, arr_padding)            
            #self.is_block_in_cell = is_in_cell #& (self.tracking_mat[0:self.for_slicing, 0] == 0)  #!  
            block_cover = sum(np.power(self.tracking_mat[0:self.for_slicing, 1][self.is_block_in_cell], 2)) #area covered by blocks
            num_starting_blocks = sum(self.tracking_mat[0:self.for_slicing, 1][self.is_block_in_cell]) / side_length #!
            cover_frac = block_cover / (dx * dy[x]) #how much of cell is COVERED by big bits
            self.cover_frac_array[x] = 1 - np.exp(-cover_frac)
            self.block_count_record[0, x] = num_starting_blocks
            if np.count_nonzero(self.tracking_mat[0:self.for_slicing, 1][self.is_block_in_cell]) == 0:
            	self.mean_block_size_record[0, x] = 0
            else:
            	self.mean_block_size_record[0, x] = np.average(self.tracking_mat[0:self.for_slicing, 1][self.is_block_in_cell]) #!
            self.cover_frac_record[0, x] = self.cover_frac_array[x]
            #self.abbrev_track_mat = self.tracking_mat[0:self.for_slicing, :]

        #get rid of nans in mean block size (these occur because of division by 0 in np.average when there are no blocks)
        #self.mean_block_size_record[0, :][np.isnan(self.mean_block_size_record[self.record_count, :])] = 0
            
    def calc_num_new_blocks(self, x):
        mean_incision_rate = np.mean(self.incision_memory[:, x])
        self.time_avg_inc_rate_array[x] = mean_incision_rate
        lam_rockfall = mean_incision_rate * self.gamma #GAMMA CONTROLS RATE OF BLOCK INPUT
        num_new_blocks = int(np.random.poisson(lam_rockfall * self.timestep, 1)) #number of new pieces
        return num_new_blocks
        
    def track_new_blocks(self, num_new_blocks, x, hillslope_block_k):
        for new in range(0, num_new_blocks):
            self.blocks += 1
            try:
                next_entry = max(max(np.where(np.isfinite(self.tracking_mat[:, 0])))) + 1 #adding one for next open entry #!
            except ValueError:
                next_entry = 0
            self.for_slicing = next_entry + 1 #b/c when you slice it takes the one before the end                
            if self.for_slicing >= self.tracking_mat.shape[0]:
                addition = np.zeros((1000, 5), dtype = np.float64) #!
                addition[:, :] = -9999
                self.tracking_mat = np.concatenate((self.tracking_mat, addition))
            else:
                pass
            #self.tracking_mat[next_entry, 0] = 0 #piece is a block #!
            self.tracking_mat[next_entry, 0] = self.cell_centers[x] #piece is dropped in at center of cell #!
            self.tracking_mat[next_entry, 1] = self.cube_side_length #!
            self.tracking_mat[next_entry, 2] = np.power(self.tracking_mat[next_entry, 1], 3) #calculate volume #!
            
            self.tracking_mat[next_entry, 4] = hillslope_block_k

        #self.abbrev_track_mat = self.tracking_mat[0:self.for_slicing, :] # should look into tis
        
    def add_block_deposition_event(self, event_number, depo_blocks_info, fluvial_block_k):
    	#depo_blocks_info is 2n (where n is # of depo events)-column array with col 0 as distance and col 1 as size
    	event_entry_distance = 2 * event_number
    	event_entry_size = 2 * event_number + 1
    	#print(depo_blocks_info)
    	#print(event_entry_distance)
    	try:
    		finite_entries = max(max(np.where(np.isfinite(depo_blocks_info[:, event_entry_distance])))) + 1
    	except ValueError:
            finite_entries = 0
    	#print(finite_entries)
    	num_new_blocks = len(depo_blocks_info[:finite_entries, event_entry_distance])
    	for new in range(0, num_new_blocks):
            self.blocks += 1
            try:
                next_entry = max(max(np.where(np.isfinite(self.tracking_mat[:, 0])))) + 1 #adding one for next open entry #!
            except ValueError:
                next_entry = 0
            self.for_slicing = next_entry + 1 #b/c when you slice it takes the one before the end                
            if self.for_slicing >= self.tracking_mat.shape[0]:
                addition = np.zeros((1000, 5), dtype = np.float64) #!
                addition[:, :] = -9999
                self.tracking_mat = np.concatenate((self.tracking_mat, addition))
            else:
                pass
            #self.tracking_mat[next_entry, 0] = 0 #piece is a block #!
            self.tracking_mat[next_entry, 0] = depo_blocks_info[new, event_entry_distance] 
            self.tracking_mat[next_entry, 1] = depo_blocks_info[new, event_entry_size]
            self.tracking_mat[next_entry, 2] = np.power(self.tracking_mat[next_entry, 1], 3) #calculate volume #!
        	#column 3 is reserved for submergence/emergence calculations
            self.tracking_mat[next_entry, 4] = fluvial_block_k
    
    def calc_flow_depth_and_velocity(self, x):
        h = self.roughness_bisection(self.q[x], self.slope[x], self.z0, self.g)  
        if np.isnan(np.sum(h)):
            print(h)
            sys.exit("NAN FOUND IN FLOW DEPTH-- KILLING MODEL")
        v = self.q[x] / h
        tau_initial = self.dens_water * self.g * h * self.slope[x] #shear stress at each node(Pa)
        self.uncorrected_tau_array[x] = tau_initial  
        if np.isnan(np.sum(self.uncorrected_tau_array)):
            print(self.uncorrected_tau_array)
            sys.exit("NAN FOUND IN UNCORRECTED TAU-- KILLING MODEL")                     
        tau = self.calc_shear_stress_with_roughness(tau_initial, self.tracking_mat, self.drag_cube, h, self.dx, self.z0, self.for_slicing, x)
        self.corrected_tau_array[x] = tau
        if np.isnan(np.sum(self.corrected_tau_array)):
            print(self.corrected_tau_array)
            sys.exit("NAN FOUND IN CORRECTED TAU-- KILLING MODEL")

        self.flow_depth[x] = h
        self.flow_velocity[x] = v
        return (h, v, tau)
        
    def erode_bed(self, rock_uplift_rate, n_cells):
        #self.surface_elev_array[-1] -= (baselevel_drop * self.timestep) #adjust baselevel node
        excess_tau = self.corrected_tau_array[0:-1] - self.tau_c_br
        #excess_tau = excess_tau.clip(min=0)
        #f_open = 1 - self.cover_frac_array[0:-1]
        #f_open = f_open.clip(min=0)
        #self.surface_elev_array[0:-1] -= self.ke_br_bed * (excess_tau) * (f_open) * self.timestep
        #print 'excess tau = ', excess_tau
        #print 'ideal erosion = ', (self.ke_br_bed * (excess_tau) * (f_open) * self.timestep)
        #print baselevel_drop
        #print self.surface_elev_array
        self.surface_elev_array[:-1] += (rock_uplift_rate * self.timestep) #do rock uplift at all but baselevel node
        #self.surface_elev_array[:20] += (baselevel_drop *self.timestep) # SMM either use this one (and the one below)or 360
        #self.surface_elev_array[20:-1] += (baselevel_drop*self.timestep) # SMM goes with top one
        #f_open = 1 - self.cover_frac_array
        #f_open = f_open.clip(min=0)
        #self.surface_elev_array[:-1] -= self.ke_br_bed * f_open[:-1] * self.corrected_tau_array[:-1] * self.timestep
        for cell in range(n_cells - 2, -1, -1):
            if excess_tau[cell] <= 0:
                pass #new elev is the same as old b/c no erosion
            else:
                f_open = 1 - self.cover_frac_array[cell]  
                f_open = f_open.clip(min=0)
                self.surface_elev_array[cell] = (self.surface_elev_array[cell] + \
                    (self.surface_elev_array[cell + 1] * f_open * self.ke_br_bed * \
                    self.dens_water * self.g * self.flow_depth[cell] * self.timestep / \
                    (self.dx * (1 + self.sigma_d_array[cell]))) + (self.timestep * f_open * \
                    self.ke_br_bed * self.tau_c_br)) / (1 + (f_open * self.ke_br_bed * \
                    self.dens_water * self.g * self.flow_depth[cell] * self.timestep / \
                    (self.dx * (1 + self.sigma_d_array[cell]))))
        if np.isnan(np.sum(self.surface_elev_array)):
            print(self.surface_elev_array)
            sys.exit("NAN FOUND IN ELEV ARRAY-- KILLING MODEL")
                
        #print self.surface_elev_array
    def calc_force_balance(self, h, v, x):
        lam_hop = (1. / 365 / 24 / 3600) #average rate of motion = 1 grain diameter per year   
        drag_force = np.zeros(len(self.tracking_mat[0:self.for_slicing, 0]), dtype = np.float64) #!
        weight_force = np.zeros(len(self.tracking_mat[0:self.for_slicing, 0]), dtype = np.float64) #!
        shear_stress_force = np.zeros(len(self.tracking_mat[0:self.for_slicing, 0]), dtype = np.float64) #!
        hop_length = np.zeros((len(self.tracking_mat[0:self.for_slicing, 0]))) #!
        
        #motion of blocks
        submerged_block = (self.is_block_in_cell) & (self.tracking_mat[0:self.for_slicing, 1] < h)#!
        emergent_block = (self.is_block_in_cell) & (self.tracking_mat[0:self.for_slicing, 1] >= h)#!
        self.tracking_mat[0:self.for_slicing, 3][submerged_block] = self.tracking_mat[0:self.for_slicing, 1][submerged_block]#!
        self.tracking_mat[0:self.for_slicing, 3][emergent_block] = h#!
        
        #FORCE BALANCE
        friction_angle = np.radians(20)
        drag_force[self.is_block_in_cell] = (1. / 2) * self.drag_cube * \
            self.tracking_mat[0:self.for_slicing, 3][self.is_block_in_cell] * \
            self.tracking_mat[0:self.for_slicing, 1][self.is_block_in_cell] * \
            self.dens_water * np.power(v, 2) #!
        weight_force[self.is_block_in_cell] = self.g * ((self.dens_sediment - \
            self.dens_water) * (self.tracking_mat[0:self.for_slicing, 3][self.is_block_in_cell] * \
            np.power(self.tracking_mat[0:self.for_slicing, 1][self.is_block_in_cell], 2)) + \
            self.dens_sediment * (self.tracking_mat[0:self.for_slicing, 1][self.is_block_in_cell] - \
            self.tracking_mat[0:self.for_slicing, 3][self.is_block_in_cell])\
            * np.power(self.tracking_mat[0:self.for_slicing, 1][self.is_block_in_cell], 2))#!
        shear_stress_force[emergent_block] = 0         
        #not_emergent_over_h = (submerged_block)# & (h > self.tracking_mat[0:self.for_slicing, 3])#!
        #emergent_over_h = (submerged_block) & (h <= self.tracking_mat[0:self.for_slicing, 3])    #!          
        shear_stress_force[submerged_block] = (self.dens_water * self.g * \
            self.slope[x] * (h - self.tracking_mat[0:self.for_slicing, 3][submerged_block])) * \
            np.power(self.tracking_mat[0:self.for_slicing, 1][submerged_block], 2)  #!        
        #shear_stress_force[emergent_over_h] = 0
        lift_force = 0.85 * (drag_force + shear_stress_force)
        slope_rad = np.arctan(self.slope[x]) 
        total_motion_force = weight_force * np.sin(slope_rad) + shear_stress_force + drag_force
        total_resist_force = (weight_force * np.cos(slope_rad) - lift_force) * np.tan(friction_angle)

        is_moving_block = (total_motion_force > total_resist_force) & (self.is_block_in_cell)
        not_moving_block = (total_motion_force <= total_resist_force) & (self.is_block_in_cell)
        hop_length[is_moving_block] = np.random.poisson(lam_hop * self.timestep, \
            sum(1 for x in is_moving_block if x)) * self.tracking_mat[0:self.for_slicing, 1][is_moving_block]#!
        hop_length[not_moving_block] = 0
        self.tracking_mat[0:self.for_slicing, 0] = self.tracking_mat[0:self.for_slicing, 0] + hop_length #adds hop length to downstream distance #!
        hop_length[:] = 0
    
    def erode_blocks(self, x):
        excess_block_tau = (self.uncorrected_tau_array[x] - self.corrected_tau_array[x] - self.tau_c_br)
        if excess_block_tau < 0:
            excess_block_tau = 0
        else:
            pass
        self.tracking_mat[0:self.for_slicing, 1][self.is_block_in_cell] -= self.tracking_mat[0:self.for_slicing, 4][self.is_block_in_cell] * excess_block_tau * self.timestep#!
        self.tracking_mat[0:self.for_slicing, 1][self.is_block_in_cell] = self.tracking_mat[0:self.for_slicing, 1][self.is_block_in_cell].clip(min = 0) #no negative blocks #!
        self.tracking_mat[0:self.for_slicing, 2][self.is_block_in_cell] = np.power(self.tracking_mat[0:self.for_slicing, 1][self.is_block_in_cell], 3)  #!
        #self.tracking_mat[0:self.for_slicing, 1:3] = self.abbrev_track_mat[:, 1:3]
            
    def calc_cover_frac(self, x):
        block_cover = sum(np.power(self.tracking_mat[0:self.for_slicing, 1][self.is_block_in_cell], 2)) #area covered by blocks    #! 
        cover_frac = block_cover / (self.dx * self.dy[x]) #how much of cell is COVERED by big bits
        self.cover_frac_array[x] = 1 - np.exp(-cover_frac)
        #get number of blocks in cell for plotting
        #self.blocks_in_cells[x] = sum(self.is_block_in_cell)
        
    #def get_block_size_classes(self, is_block_in_cell, x):
    #    sizes = self.tracking_mat[0:self.for_slicing, 1][is_block_in_cell] #!
    #    if np.isnan(np.sum(sizes)):
    #        print(sizes)
    #        print(self.tracking_mat[:, 1]) #!
    #        sys.exit("NAN FOUND IN SIZES ARRAY-- KILLING MODEL")
    #    self.zero_to_one_temp[x] = sum(sizes[sizes <= 1])
    #    self.one_to_two_temp[x] = sum(sizes[(sizes > 1) & (sizes <= 2)])
    #    self.two_to_three_temp[x] = sum(sizes[(sizes > 2) & (sizes <= 3)])
    #    self.three_to_four_temp[x] = sum(sizes[(sizes > 3) & (sizes <= 4)])
    #    self.four_to_five_temp[x] = sum(sizes[(sizes > 4) & (sizes <= 5)])
    #    self.five_to_six_temp[x] = sum(sizes[(sizes > 5) & (sizes <= 6)])
        
    def track_vertical_erosion(self, rock_uplift_rate):
        #COMMENTED OUT 11/25/2016 to see what results look like
        #if np.any(self.surface_elev_array > self.old_surf_elev_array):
        #    print self.surface_elev_array
        #    print self.old_surf_elev_array
        #    print self.cover_frac_array
        #    print self.ke_br_bed
        #    print self.flow_depth
        #    print self.sigma_d_array
        #    print self.tau_c_br
        #    print self.timestep
        #    print self.dx
        #    print self.g
        #    print self.dens_water
        #    os.system("rm results.out")
        #    os.system("touch results.out")
        #    with open("results.out", "a") as f: 
        #        f.write('FAIL\n')
        #    sys.exit("NUMERICAL AGGRADATION-- KILLING MODEL")
        #self.incision_rate_array = abs(self.surface_elev_array - self.old_surf_elev_array)/ self.timestep
        self.incision_rate_array = abs(((self.old_surf_elev_array + np.append(rock_uplift_rate, 0) * self.timestep) - self.surface_elev_array))/ self.timestep
        #print(self.incision_rate_array)
        self.incision_memory = np.roll(self.incision_memory, 1, 0)
        self.incision_memory[0, :] = self.incision_rate_array
        self.old_surf_elev_array[:] = self.surface_elev_array #update old, for incision rate calculations
        
    def display_model_progress(self):
        if self.run_count == 1 or np.remainder(self.run_count, 1) == 0: # display run count
            print('Run ' + str(self.run_count) + ' of ' + str(self.number_iterations))
            #print(self.blocks)
        else: 
            pass
        
    def update_slope(self):
        self.slope[0:-1] = -(self.surface_elev_array[1:] - self.surface_elev_array[0:-1]) / (self.dx)
        self.slope[-1] = self.initial_slope 
        
    def update_record_keeping_arrays(self):
    	self.record_count += 1            
    	self.surface_elev_record[self.record_count, :] = self.surface_elev_array
    	self.slope_record[self.record_count, :] = self.slope
    	self.uncorrected_tau_record[self.record_count, :] = self.uncorrected_tau_array
    	self.corrected_tau_record[self.record_count, :] = self.corrected_tau_array
    	self.time_avg_inc_rate_record[self.record_count, :] = self.time_avg_inc_rate_array
    	#self.block_count_record[self.record_count, :] = sum(self.is_block_in_cell)
    	#self.mean_block_size_record[self.record_count, :] = np.average(self.tracking_mat[0:self.for_slicing, 1][self.is_block_in_cell])
    	self.cover_frac_record[self.record_count, :] = self.cover_frac_array
    	
    	#spatial loop, needed to do block count and mean size at every cell ONLY DURING RECORDED TIMES
    	for x in self.x_array:
    		self.is_block_in_cell = (self.tracking_mat[0:self.for_slicing, 0] >= self.upstream_edges_of_cell[x]) & (self.tracking_mat[0:self.for_slicing, 0] < self.downstream_edges_of_cell[x]) #!
    		self.block_count_record[self.record_count, x] = sum(self.is_block_in_cell)
    		if np.count_nonzero(self.tracking_mat[0:self.for_slicing, 1][self.is_block_in_cell]) == 0:
    			self.mean_block_size_record[self.record_count, x] = 0
    		else:
    			self.mean_block_size_record[self.record_count, x] = np.average(self.tracking_mat[0:self.for_slicing, 1][self.is_block_in_cell]) #!
    		#self.mean_block_size_record[self.record_count, x] = np.average(self.tracking_mat[0:self.for_slicing, 1][self.is_block_in_cell])
    		
    	#get rid of nans in mean block size (these occur because of division by 0 in np.average when there are no blocks)
    	#self.mean_block_size_record[self.record_count, :][np.isnan(self.mean_block_size_record[self.record_count, :])] = 0
        
        #self.one_to_two[self.record_count, :] = self.one_to_two_temp
        #self.two_to_three[self.record_count, :] = self.two_to_three_temp
        #self.three_to_four[self.record_count, :] = self.three_to_four_temp
        #self.four_to_five[self.record_count, :] = self.four_to_five_temp
        #self.five_to_six[self.record_count, :] = self.five_to_six_temp
        
    def delete_eroded_or_gone_blocks(self):
        #now find and delete any tracking rows where blocks have degraded to 0.
        for_abrasion_deletion = np.where(self.tracking_mat[:, 1] == 0)[0] #for_deletion is an array of vertical indices #!
        self.tracking_mat = np.delete(self.tracking_mat, for_abrasion_deletion, axis = 0) #deletes relevant rows
        
        #now find and delete any tracking rows where blocks have left domain
        #if self.run_count > 260:
        #    print(self.tracking_mat[:, 0])
        #    print(self.downstream_edges_of_cell[-1])
        for_transport_deletion = np.where(self.tracking_mat[:, 0] > self.downstream_edges_of_cell[-1])[0] #!
        self.tracking_mat = np.delete(self.tracking_mat, for_transport_deletion, axis = 0)
        
    def calc_and_return_model_metrics(self, dx, dy, domain_length):
        #total volume lost
        total_volume_lost = sum(self.elev[0, :] - self.elev[-1, :]) * dx * dy[:]
        
        #final max elevation
        final_max_elevation = np.max(self.elev[-1, :])      
        
        #final mean elevation
        final_mean_elevation = np.average(self.elev[-1, :])        
        
        #standard deviation of elevation
        final_stdev_elevation = np.std(self.elev[-1, :])
        
        #final maximum slope
        final_max_slope = np.max(self.slope[-1, :])
        
        #final mean slope
        final_mean_slope = np.average(self.slope[-1, :])        
        
        #final standard deviation of slope
        final_stdev_slope = np.std(self.slope[-1, :])
        
        #total number of blocks
        total_blocks = self.blocks
        
        #final maximum block size
        
        #final average block size
        
        #final stdev of block size
        
        #something about cover
        
        #something about shear stress reduction
        outputs_list = [total_volume_lost, final_max_elevation, final_mean_elevation,
            final_stdev_elevation, final_max_slope, final_mean_slope, final_stdev_slope, total_blocks]
        os.system("touch outputs_for_analysis.txt")
        for i in range(len(outputs_list)):
            line = str(outputs_list[i]) + '\n'
            with open("outputs_for_analysis.txt", "a") as f: 
                f.write(line)          
        
        #with open('outputs_for_analysis.txt', 'a') as fp:
        #    fp.write(total_volume_lost+'\n'+final_max_elevation+'\n'+final_mean_elevation+'\n'+final_stdev_elevation+'\n'+final_max_slope+'\n'+final_mean_slope+'\n'+final_stdev_slope+'\n')
    
        return (total_volume_lost, final_max_elevation, final_mean_elevation,
            final_stdev_elevation, final_max_slope, final_mean_slope, final_stdev_slope, total_blocks)    
    
    def initialize(self, time_to_run_yrs, timestep_yrs, record_time_int, imposed_rock_uplift_rate, delay_timescale, gamma, n_cells, dx, side_length, bed_k, hillslope_block_k, fluvial_block_k, z0, tau_c, depo_blocks_info, depo_events_timing, chan_width, k_w, hack_c, hack_exp):
        #inputs
        self.n_cells = n_cells #number of cells
        self.dx = dx #cell side length (m)
        
        #channel width is an array now
        self.x_increase_from_left = np.arange(0, n_cells * dx, dx) #cell width (m)
        self.dy = chan_width + k_w * np.power(self.x_increase_from_left / hack_c, 0.5 / hack_exp)
        
        self.record_time_int = record_time_int #years
        self.calculate_timing_and_recording_vars(self.record_time_int, time_to_run_yrs, timestep_yrs, imposed_rock_uplift_rate, delay_timescale)
        self.instantiate_common_constants()        
        self.initial_slope = .00001 #initial bed slope at each point (m/m)
        self.gamma = gamma#0 #degree of channel-hillslope connectivity. Higher means stronger block effects. set to 0 for no block effects.
        starting_dist = 0
        starting_height_adjustment = 10 #starting elev of base of profile [m]
        self.cube_side_length = side_length #m #for cubes
        #boulders = 0. #how many boulders? tracked as 1's
        number_of_blocks = 0. #how many cubes? tracked as 0's
        #number_of_pieces = boulders + cubes #This is how many pieces are initially distributed

        self.instantiate_tracking_matrix(number_of_blocks, self.gamma, starting_dist, fluvial_block_k)

        #import discharge distribution        
        #self.discharges = self.import_discharges('discharges_for_model.csv', self.number_iterations)        
        
        #Shear stress formulation parameters:
        self.set_shear_stress_params(bed_k, tau_c)
        
        #user-defined variables relating to bed layering
        self.define_cell_centers_and_edges(self.n_cells, self.dx, self.initial_slope, starting_height_adjustment)
        
        #instantiate arrays that change with each timestep
        self.instantiate_single_dt_arrays(self.n_cells, self.initial_slope, self.surface_elev)
        
        #instantiate record-keeping arrays
        self.instantiate_record_keeping_arrays(self.time_to_run, self.record_time_seconds, self.number_records, self.n_cells, self.surface_elev_array, self.slope)
        
        #set roughness calculation parameters (for roughness generated by bed, not boulders)
        self.set_roughness_depth_calc_params(z0)        
        
        #instantiate counters for loops
        #from cfuncs import instantiate_counting_parameters
        self.instantiate_counting_parameters()

        #calculate initial values of arrays such as block_count_record and cover_frac_record
        self.calculate_initial_array_values(number_of_blocks, self.x_array, self.dx, self.dy, self.cube_side_length)        
        self.extra_q_iterator = 0 #instantiate here
    def update(self, weibull_k, weibull_scale, hillslope_block_k, fluvial_block_k, depo_blocks_info, depo_events_timing, area_0, hack_c, hack_exp):
        #from cfuncs import calc_num_new_blocks 
        self.time += self.timestep
        self.run_count += 1
        #baselevel_drop = self.baselevel_drop_annual[self.run_count - 1]
        #rock_uplift_rate = self.rock_uplift_rate
        #pull discharge from q_distribution
        #q_volume = self.discharges[self.run_count - 1]
        #weibull_mean = 10.0#15.0384 #m3/s discharge at orodell gauge
        #weibull_k = 0.5 #unitless, c in Rossi et al (2016)
        #weibull_scale = 5.0#16.94664 #calc'ed outside with sp.special.gamma()
        q_volume = weibull_scale * np.random.weibull(weibull_k) #Q at left hand side
        r = q_volume / area_0 #get runoff
        self.area_array = area_0 + np.power(self.x_increase_from_left / hack_c, 1 / hack_exp)
        self.q = (r * self.area_array) / self.dy
        
        for x in self.x_array: #spatial loop        
            #allow rockfall to put in new blocks
            #print self.q
            num_new_blocks = self.calc_num_new_blocks(x)
            self.track_new_blocks(num_new_blocks, x, hillslope_block_k)
        
            #now, need to redo is_block, is_boulder, etc
            self.is_block_in_cell = (self.tracking_mat[0:self.for_slicing, 0] >= self.upstream_edges_of_cell[x]) & (self.tracking_mat[0:self.for_slicing, 0] < self.downstream_edges_of_cell[x]) #!
            #arr_padding = np.repeat(False, len(self.tracking_mat[:, 0]) - self.for_slicing)  
            #is_block_in_cell = np.append(is_block_in_cell, arr_padding)
            #is_block_in_cell = is_in_cell #& (self.tracking_mat[0:self.for_slicing, 0] == 0)   #! 
            
            #calculate flow depth and flow velocity using VPE
            self.calc_flow_depth_and_velocity(x)   
            #calc cover frac so it is using same information as roughness 
            self.calc_cover_frac(x)
        #erode the bed (Fastscape-style algorithm)
        self.erode_bed(self.rock_uplift_rate, self.n_cells)
            
        for x in self.x_array: #spatial loop  
            self.calc_force_balance(self.flow_depth[x], self.flow_velocity[x], x)
            
            #calculate change in cover fraction based on blocks that moved
            self.is_block_in_cell = (self.tracking_mat[0:self.for_slicing, 0] >= self.upstream_edges_of_cell[x]) & (self.tracking_mat[0:self.for_slicing, 0] < self.downstream_edges_of_cell[x]) #!
            #is_block_in_cell = is_in_cell #& (self.tracking_mat[0:self.for_slicing, 0] == 0)  #!
            #arr_padding = np.repeat(False, len(self.tracking_mat[:, 0]) - self.for_slicing)  
            #is_block_in_cell = np.append(is_block_in_cell, arr_padding)
            if np.isnan(np.sum(self.is_block_in_cell)):
                print(self.is_block_in_cell)
                print(self.tracking_mat[:, 1]) #!
                sys.exit("NAN FOUND IN IS_BLOCK ARRAY-- KILLING MODEL")
            #abrade blocks with excess shear stress
            self.erode_blocks(x)   

            #calculate change in cover fraction based on blocks that were abraded          
            #self.calc_cover_frac(is_block_in_cell, x)
            #11/25/16 commenting this out b/c it isnt the cover frac that is affecting erosion
            #get size classes of blocks in cells
            #self.get_block_size_classes(is_block_in_cell, x)
            
        
        #track vertical erosion to keep incision rate matrix up-to-date
        self.track_vertical_erosion(self.rock_uplift_rate)
        
        #if time equals one of the times in the megaflood occurrence array:
        if self.time / (3600*24*365) in depo_events_timing:
        	print('depo event triggered')
        	event_number = np.where(depo_events_timing == (self.time / (3600*24*365)))[0][0]
        	self.add_block_deposition_event(event_number, depo_blocks_info, fluvial_block_k)
        
        #print model progress to screen at given intervals
        self.display_model_progress()
        
        #update slope for next iteration (downwind scheme)
        self.update_slope()

        #update record keeping arrays   
        if np.remainder(self.run_count, (self.record_time_int / self.timestep_yrs)) == 0: #
            self.update_record_keeping_arrays()
        else:
            pass
        
        self.delete_eroded_or_gone_blocks()
        
    def finalize(self, suffix):
        #once everything has run, spit out arrays as .npy binaries
        #print self.blocks
        np.save('time_record_' + suffix, self.time_record)
        np.save('elev_record_' + suffix, self.surface_elev_record)
        np.save('slope_record_' + suffix, self.slope_record)
        np.save('block_count_record_' + suffix, self.block_count_record)
        np.save('mean_block_size_record_' + suffix, self.mean_block_size_record)
        np.save('cover_frac_record_' + suffix, self.cover_frac_record)
        np.save('time_avg_inc_rate_record_' + suffix, self.time_avg_inc_rate_record)
        np.save('uncorrected_tau_record_' + suffix, self.uncorrected_tau_record)
        np.save('corrected_tau_record_' + suffix, self.corrected_tau_record)
        #np.save('zero_to_one_' + suffix, self.zero_to_one)
        #np.save('one_to_two_' + suffix, self.one_to_two)
        #np.save('two_to_three_' + suffix, self.two_to_three)
        #np.save('three_to_four_' + suffix, self.three_to_four)
        #np.save('four_to_five_' + suffix, self.four_to_five)
        #np.save('five_to_six_' + suffix, self.five_to_six)
        
    ##############################################################################
    def run_model(self, time_to_run, timestep, record_time_int, imposed_rock_uplift_rate, delay_timescale, gamma, side_length, suffix, n_cells, dx, bed_k, hillslope_block_k, fluvial_block_k, z_0, tau_c, weibull_k, weibull_scale, depo_blocks_info, depo_events_timing, chan_width, k_w, hack_c, hack_exp, area_0):
        self.initialize(time_to_run, timestep, record_time_int, imposed_rock_uplift_rate, delay_timescale, gamma, n_cells, dx, side_length, bed_k, hillslope_block_k, fluvial_block_k, z_0, tau_c, depo_blocks_info, depo_events_timing, chan_width, k_w, hack_c, hack_exp)
        while self.time < self.time_to_run:
            self.update(weibull_k, weibull_scale, hillslope_block_k, fluvial_block_k, depo_blocks_info, depo_events_timing, area_0, hack_c, hack_exp)
        self.finalize(suffix)
