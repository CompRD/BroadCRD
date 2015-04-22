#!/usr/bin/env python

###############################################################################
#                    SOFTWARE COPYRIGHT NOTICE AGREEMENT                      #
#        This software and its documentation are copyright (2011) by the      #
#    Broad Institute.  All rights are reserved.  This software is supplied    #
#    without any warranty or guaranteed support whatsoever. The Broad         #
#    Institute is not responsible for its use, misuse, or functionality.      #
###############################################################################

# A simulation program to help debug the problems with calculating error
# bars on the latent variance (aka "jaffe statistic") of read starts computed
# by BadCoverage.

# Original author: Michael G. Ross <mgross@broadinstitute.org>

import random
from scipy.stats import kurtosis, poisson

def generate_event(region_boundaries):
    which_region = random.randint(0, len(region_boundaries) - 1)
    base = 0
    if which_region > 0:
        base = region_boundaries[which_region - 1]
    return random.randint(base, region_boundaries[which_region] - 1)

def generate_counts(region_boundaries, num_events):
    counts = [0] * region_boundaries[-1]
    for e in range(num_events):
        counts[generate_event(region_boundaries)] += 1
    return counts

def compute_sample_mean(samples):
    sample_mean = 0
    for s in samples:
        sample_mean += s
    sample_mean /= float(len(samples))
    return sample_mean

def compute_sample_variance(samples):
    sample_mean = compute_sample_mean(samples)
    sample_var = 0
    for s in samples:
        sample_var += (s - sample_mean) ** 2
    sample_var /= float(len(samples) - 1)
    return sample_var

def compute_latent_variance(region_boundaries):
    mean_prob = 1 / float(region_boundaries[-1])
    var = 0
    prev_bound = 0
    for r in region_boundaries:
        rsize = r - prev_bound
        rprob = 1 / float(len(region_boundaries))
        var += rsize * (rprob / rsize - mean_prob) ** 2
        prev_bound = r
    var /= region_boundaries[-1]
    return var

def main():
    prob_region_boundaries = [500, 750, 2500, 5000, 10000]
    num_reads = 50000
    genome_size = prob_region_boundaries[-1]
    latent_variance = compute_latent_variance(prob_region_boundaries)
    latent_mean = 1 / float(genome_size)
    print 'true latent mean: {0}, latent variance: {1}'.format(latent_mean,
        latent_variance)
    print 'true jaffe stat: {0}'.format(latent_variance / float(latent_mean**2))
    var_samples = []
    est_var_samples_var = []
    jaffe_samples = []
    jaffe_est_vars = []
    jaffe_est_vars2 = []
    for k in range(100):
        #count_samples = generate_counts(prob_region_boundaries, num_reads)
        count_samples = poisson.rvs(num_reads / float(genome_size),
            size=num_reads)
        
        count_mean = compute_sample_mean(count_samples)
        count_var = compute_sample_variance(count_samples)
        var_samples.append(count_var)

        jaffe_stat = (((count_var / count_mean**2) - (1 / count_mean)) +
            (1 / float(num_reads))) * (float(num_reads) / float(num_reads - 1))
        jaffe_samples.append(jaffe_stat)
   
   
        count_kurtosis = kurtosis(count_samples)
        e_vsv = count_var ** 2 * (2 / float(len(count_samples) - 1) +
            count_kurtosis / float(len(count_samples)))
        est_var_samples_var.append(e_vsv)
        
        jvar = (float(num_reads**2) / float((num_reads - 1)**2) *
            count_mean**(-4) * e_vsv)

        jaffe_est_vars.append(jvar)

    var_samples_mean = compute_sample_mean(var_samples)
    var_samples_var = compute_sample_variance(var_samples)
    
    print ('count var mean: {0}, count var variance: {1}'
        .format(var_samples_mean, var_samples_var))
    mean_est_var_samples_var = compute_sample_mean(est_var_samples_var)
    print 'mean var estimate: {0}'.format(mean_est_var_samples_var)
    
    #for v in range(len(est_var_samples_var)):
    #    print (est_var_samples_var[v] / var_samples_var)
    
    mean_jaffe_stat = compute_sample_mean(jaffe_samples)
    var_jaffe_stat = compute_sample_variance(jaffe_samples)
    print 'mean jaffe stat: {0}, var jaffe stat: {1}'.format(mean_jaffe_stat,
        var_jaffe_stat)
    mean_jaffe_est_vars = compute_sample_mean(jaffe_est_vars)
    var_jaffe_est_vars = compute_sample_variance(jaffe_est_vars)
    print ('mean jaffe var est: {0} var jaffe var est: {1}'
        .format(mean_jaffe_est_vars, var_jaffe_est_vars))
    
    return 0

if __name__ == '__main__':
    exit(main())