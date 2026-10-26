"""Small numerical invariants; full research experiments are not CI tests."""
import unittest
import numpy as np
from physics_model import Column,solve
from fast_column_sensitivity import checks

class NumericalChecks(unittest.TestCase):
    def test_reference_mass_and_fields(self):
        p=Column();r=solve(p,32,times=np.linspace(0,1,21));storage=np.mean(r['c']+p.beta*r['s'],axis=1)
        self.assertLess(np.max(abs(storage-r['integrated_net_flux'])),1e-6)
        self.assertGreaterEqual(float(r['c'].min()),-1e-7)
        self.assertLessEqual(float(r['c'].max()),1+1e-6)
        self.assertTrue(np.all(r['c'][0]==0))
    def test_forward_sensitivity(self):
        result=checks();self.assertLess(result['directional_jacobian_relative_error'],1e-3)

if __name__=='__main__':unittest.main()
