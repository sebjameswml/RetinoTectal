import os

# On corebest, 34 cores was fastest, 36 almost as fast.
for cores in range(32,37):
    print ('\nTesting {0} cores'.format(cores))
    os.environ['OMP_NUM_THREADS'] = '{0}'.format(cores);
    os.system('time ./build/sim/agent/agent1c configs/a1/m_ee_GJ.json  configs/a1/e_wt.json 1000');
