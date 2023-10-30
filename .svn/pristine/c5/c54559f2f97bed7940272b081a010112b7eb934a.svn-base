#!/bin/bash
#
# Script to illustrate running batch jobs and passing in arguments.
#
# 
# This script assumes that the following has been run successfully:
# scons co=1 b=GccOpt ts=projects/CryptProliferation2013/test/TestCryptProliferationDistributionLiteratePaperRunner.hpp
#


end_t=2200 #2200
n_wnt_sweeps=10 #10
n_gen_sweeps=10 #10
n_ci_sweeps=10 #10



echo "End Time = $end_t, Num Sweeps = $n_wnt_sweeps"

nice -20 ../build/optimised/TestCryptProliferationDistributionLiteratePaperRunner -CI -end_time $end_t  -CCM 1 -min 1 -max 6 -num_sweeps $n_gen_sweeps  -min_CI 0 -max_CI 1 -num_CI_sweeps $n_ci_sweeps > output/Model1_output.txt 2>&1 &

nice -20 ../build/optimised/TestCryptProliferationDistributionLiteratePaperRunner -CI -end_time $end_t  -CCM 2 -min 0 -max 1 -num_sweeps $n_wnt_sweeps  -min_CI 0 -max_CI 1 -num_CI_sweeps $n_ci_sweeps > output/Model2_output.txt 2>&1 &

nice -20 ../build/optimised/TestCryptProliferationDistributionLiteratePaperRunner -CI -end_time $end_t  -CCM 3 -min 0 -max 1 -num_sweeps $n_wnt_sweeps  -min_CI 0 -max_CI 1 -num_CI_sweeps $n_ci_sweeps > output/Model3_output.txt 2>&1 &

nice -20 ../build/optimised/TestCryptProliferationDistributionLiteratePaperRunner -WDCCD -CI -end_time $end_t  -CCM 1 -min 1 -max 6 -num_sweeps $n_gen_sweeps  -min_CI 0 -max_CI 1 -num_CI_sweeps $n_ci_sweeps > output/Model1a_output.txt 2>&1 &

nice -20 ../build/optimised/TestCryptProliferationDistributionLiteratePaperRunner -WDCCD -CI -end_time $end_t  -CCM 2 -min 0 -max 1 -num_sweeps $n_wnt_sweeps  -min_CI 0 -max_CI 1 -num_CI_sweeps $n_ci_sweeps > output/Model2a_output.txt 2>&1 &

nice -20 ../build/optimised/TestCryptProliferationDistributionLiteratePaperRunner -WDCCD -CI -end_time $end_t  -CCM 3 -min 0 -max 1 -num_sweeps $n_wnt_sweeps  -min_CI 0 -max_CI 1 -num_CI_sweeps $n_ci_sweeps > output/Model3a_output.txt 2>&1 &

echo "Jobs submitted"


