"""
Goal/Purpose: To simulate kaon decays to pi+ and pi- particles and reconstruct the z positions
of particles that successfully enter the trigger space of the calorimeter. These values are then
placed in a histogram for comparison to exiting histogram data sets.
Then, to find a minimum chi-squared value between two sets of bin arrays and the 
normalization and tau values for which this minimum occurs. This gives an estimate for lifetime.
    
Assignment: Unit 10 Problem 2: Fitting for the kaon mean lifetime
    
Author: Jamie Howell

Collaborators: Professor George Gollin

File: Kaon Decay Project.py

Date: December 7, 2017

References: Physics 298owl course notes / in-class exercise solutions

"""

import time
import unit08_generatedecayObject as GDO
import unit08_histogramObject as hb
import numpy as np

# instantiate the kaon decay generator
KDG = GDO.Kdecay()

# initialize variables and z positions
number_of_decays_so_far = 0
number_of_decays_to_generate = 100000
number_of_triggers = 0
z_target = 0
z_chamber1 = z_target + 61.
z_chamber2 = z_chamber1 + 10
z_calorimeter = z_chamber1 + 15
z_hole = z_calorimeter

# half-widths
halfwidth_target = 0
halfwidth_chamber1 = 0.6
halfwidth_chamber2 = 0.7
halfwidth_calorimeter = 0.75
halfwidth_hole = 0.25

# total kaons
number_of_kaons_total = 100000

# tracking chamber resolution
chamber_resolution = 0.2e-3
        
# pi plus and pi minus position arrays
piplus_x = np.array([np.nan, np.nan, np.nan, np.nan, np.nan])
piplus_y = np.array([np.nan, np.nan, np.nan, np.nan, np.nan])
piminus_x = np.array([np.nan, np.nan, np.nan, np.nan, np.nan])
piminus_y = np.array([np.nan, np.nan, np.nan, np.nan, np.nan])

# create arrays for z positions and half-widths
z_array = np.array([z_target, z_chamber1, z_chamber2, z_calorimeter, z_hole])
half_width_array = np.array([halfwidth_target, halfwidth_chamber1, halfwidth_chamber2, \
halfwidth_calorimeter, halfwidth_hole]) 

# create our histogram
histo1 = hb.histo("z reconstructed, funny lifetime", 300, 0., 30.00)
histo1.hsetlabels("reconstructed decay z, meters", "N")

#starting time
print("start generating events at ", time.ctime())

# loop over each decay
while number_of_decays_so_far < number_of_decays_to_generate:
    
    number_of_decays_so_far = number_of_decays_so_far + 1     
          
    # get the vertex and four momentum values for pi plus and minus
    vertex, pmu_plus, pmu_minus = KDG.getdecay()

   
    # loop over each tracker's plane
    for plane in range(1,5):
        
        # z separation between decay point and plane under consideration
        dz = z_array[plane] - vertex[2]
        
        # change in x for the two pions between vertex and this z location
        dxplus = (pmu_plus[1] / pmu_plus[3]) * dz
        dxminus = (pmu_minus[1] / pmu_minus[3]) * dz

        # change in y for the two pions between vertex and this z location
        dyplus = (pmu_plus[2] / pmu_plus[3]) * dz
        dyminus = (pmu_minus[2] / pmu_minus[3]) * dz
        
        # store these in arrays now:
        piplus_x[plane] = vertex[0] + dxplus
        piplus_y[plane] = vertex[1] + dyplus
        piminus_x[plane] = vertex[0] + dxminus
        piminus_y[plane] = vertex[1] + dyminus

    
    #check to see if the trigger was satisfied for pi plus nd minus
    piplus_OK = \
    np.abs(piplus_x[1]) < half_width_array[1] and \
    np.abs(piplus_y[1]) < half_width_array[1] and \
    np.abs(piplus_x[2]) < half_width_array[2] and \
    np.abs(piplus_y[2]) < half_width_array[2] and \
    np.abs(piplus_x[3]) < half_width_array[3] and \
    np.abs(piplus_y[3]) < half_width_array[3] and \
    (np.abs(piplus_x[4]) > half_width_array[4] or \
    np.abs(piplus_y[4]) > half_width_array[4])
    
    piminus_OK = \
    np.abs(piminus_x[1]) < half_width_array[1] and \
    np.abs(piminus_y[1]) < half_width_array[1] and \
    np.abs(piminus_x[2]) < half_width_array[2] and \
    np.abs(piminus_y[2]) < half_width_array[2] and \
    np.abs(piminus_x[3]) < half_width_array[3] and \
    np.abs(piminus_y[3]) < half_width_array[3] and \
    (np.abs(piminus_x[4]) > half_width_array[4] or \
    np.abs(piminus_y[4]) > half_width_array[4])
 
    trigger_OK = piplus_OK and piminus_OK and vertex[2] < z_chamber2
                                                    
    pair_mass = 1000.0 * np.sqrt((pmu_plus[0] + pmu_minus[0])**2 - \
                                (pmu_plus[1] + pmu_minus[1])**2 - \
                                (pmu_plus[2] + pmu_minus[2])**2 - \
                                (pmu_plus[3] + pmu_minus[3])**2 )

    if trigger_OK:

        number_of_triggers = number_of_triggers + 1

        # smear chamber positions
        xA1 = piplus_x[1] + np.random.normal(0.0, chamber_resolution) 
        xA2 = piplus_x[2] + np.random.normal(0.0, chamber_resolution) 
        xB1 = piminus_x[1] + np.random.normal(0.0, chamber_resolution) 
        xB2 = piminus_x[2] + np.random.normal(0.0, chamber_resolution) 
        yA1 = piplus_y[1] + np.random.normal(0.0, chamber_resolution) 
        yA2 = piplus_y[2] + np.random.normal(0.0, chamber_resolution) 
        yB1 = piminus_y[1] + np.random.normal(0.0, chamber_resolution) 
        yB2 = piminus_y[2] + np.random.normal(0.0, chamber_resolution) 

        new_z_chamber1 = z_chamber1
        new_z_chamber2 = z_chamber2

        # calculate z values of the crossings of tracks
        zx = (xA1 - xB1) * (new_z_chamber2 - new_z_chamber1) / (xB2 - xB1 - xA2 + xA1) \
        + new_z_chamber1
        zy = (yA1 - yB1) * (new_z_chamber2 - new_z_chamber1) / (yB2 - yB1 - yA2 + yA1) \
        + new_z_chamber1

        # find x and y angles
        theta_x = ((xA2 - xB2) - (xA1 - xB1)) / (new_z_chamber2 - new_z_chamber1)
        theta_y = ((yA2 - yB2) - (yA1 - yB1)) / (new_z_chamber2 - new_z_chamber1)

        # find the reconstructed value
        z_reconstructed = (zx * theta_x**2 + zy * theta_y**2) / (theta_x**2 + theta_y**2)
        
        # fill our histogram
        histo1.hfill(z_reconstructed)

# bin population for unknown events
in_binpop = np.array([0] * 300)
in_data_file = open("unit10_homework2_GG.dat", "rb")
in_binpop = np.load(in_data_file)

#initialize more variables
chi_squared_min = 1.0e12
gamma = KDG.e_beam / KDG.m_kaon
z_bin_width = 0.1
trial_tau_start = 8.0e-11
trial_tau_stop = 15.0e-11
trial_tau_increment = 0.1e-11
normalization_start = 3.0
normalization_stop = 30.0
normalization_increment = 0.1

trial_tau = trial_tau_start
    
while trial_tau < trial_tau_stop:

    normalization_ratio = normalization_start

    while normalization_ratio < normalization_stop:
        # mystery sample-to-my sample ratio:
        
        # calculate delta values
        given_delta = gamma * KDG.clight * trial_tau
        my_delta = gamma * KDG.clight * 8.954e-11

        chi_squared_sum = 0.0
        
        #loop over bins
        for i in range(0, 300):       
            #z at center of bin
            z_center = (i + 0.5) * z_bin_width
            multiplier = normalization_ratio * np.exp(-z_center * (1.0/given_delta - 1.0/my_delta))
        
            if histo1.binpop[i] > 0:
                chi_squared_sum = chi_squared_sum + \
                (histo1.binpop[i] - \
                 in_binpop[i] / multiplier)**2 / histo1.binpop[i]
                
        if chi_squared_min > chi_squared_sum:
            chi_squared_min = chi_squared_sum
            best_tau = trial_tau
            best_normalization = normalization_ratio
               
        normalization_ratio = normalization_ratio + normalization_increment 
    
    trial_tau = trial_tau + trial_tau_increment
        
    
print("Chi-squared minimum is: ", chi_squared_min)
print("Best guess at lifetime is: ", best_tau)
print("Best guess at relative normalization is: ", best_normalization)
print("Ending time is: ", time.ctime())
