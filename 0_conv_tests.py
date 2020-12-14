"""
automag.0_conv_tests
==================================

Script which runs convergence tests.

.. codeauthor:: Michele Galasso <michele.galasso@skoltech.ru>
"""

from copy import copy

from common.ConvergenceTest import ConvergenceTest

# define the VASP parameters
params = {
    'xc': 'PBE',
    'setups': 'recommended',
    'istart': 0,
    'prec': 'Normal',
    'ncore': 4,
    'ediff': 1e-6,
    'ismear': 1,
    'nbands': 148,
    'lcharg': False,
    'lwave': False,
}

# convergence test w.r.t. encut
params1 = copy(params)
params1['sigma'] = 0.1
params1['kpts'] = 40

encut_values = range(400, 800, 10)

convtest = ConvergenceTest(params1, encut_values=encut_values)
convtest.submit()

# convergence test w.r.t. sigma and kpts
params2 = copy(params)
params2['encut'] = 600

sigma_values = [0.05, 0.1, 0.2]
kpts_values = range(20, 110, 10)

convtest = ConvergenceTest(params2, sigma_values=sigma_values, kpts_values=kpts_values)
convtest.submit()
