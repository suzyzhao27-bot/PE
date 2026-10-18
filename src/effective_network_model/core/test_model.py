"""Mathematical acceptance checks for the shared effective-field model."""
import math
import unittest
from theory_model import (field,log_resistance_nm,log_integral_exp_field,thickness_log_slope,
                          neutral_profile,finite_bath,normalized_finite_bath)


class PhysicalChecks(unittest.TestCase):
    def test_zero_field_returns_geometric_length(self):
        for h in (65.,350.,2600.):
            self.assertAlmostEqual(math.exp(log_resistance_nm(h,0.,900.,25.)),h+math.pi*25/2,places=8)

    def test_constant_transfer_energy_limit(self):
        # ell -> infinity gives spatially constant w=2*ws.
        h,ws=350.,2.
        logr=log_integral_exp_field(h,ws,1e14)
        self.assertAlmostEqual(math.exp(math.log(h)-logr),math.exp(-2*ws),places=10)

    def test_robin_boundary_and_ode(self):
        # Independent derivative formulas for B=1,C=1/ell,H=2*ws/ell.
        for h in (65.,350.,2600.):
            ell,ws=894.,7.04
            for z in (0.,h/3,h):
                e1,e2=math.exp(-z/ell),math.exp(-(h-z)/ell)
                q=field(z,h,ws,ell)
                qpp=ws*(e1+e2)/(ell*ell)
                self.assertAlmostEqual(qpp*ell*ell,q,places=12)
            qp0=ws*(-1+math.exp(-h/ell))/ell
            self.assertAlmostEqual(-qp0+field(0,h,ws,ell)/ell,2*ws/ell,places=14)

    def test_reflection_symmetry(self):
        for r in (-.5,0.,1.):
            self.assertAlmostEqual(field(23,350,7,894,r),field(327,350,7,894,r),places=12)

    def test_quadrature_refinement(self):
        for h in (65.,350.,2600.,10000.):
            r1=log_integral_exp_field(h,7.04,894.,n=2048)
            r2=log_integral_exp_field(h,7.04,894.,n=4096)
            self.assertLess(abs(math.expm1(r1-r2)),2e-7)

    def test_positive_barrier_adds_resistance(self):
        self.assertGreater(log_resistance_nm(350,3,894),log_resistance_nm(350,0,894))

    def test_thickness_derivative_against_finite_difference(self):
        for h in (65.,350.,700.,2600.):
            eps=1e-5
            left=log_resistance_nm(h*math.exp(-eps),7.04,894.)
            right=log_resistance_nm(h*math.exp(eps),7.04,894.)
            numerical=-(right-left)/(2*eps)
            self.assertLess(abs(numerical-thickness_log_slope(h,7.04,894.)),2e-7)

    def test_interface_partition_not_bath_concentration(self):
        p=neutral_profile(350,3,894,270,1,1000)
        self.assertAlmostEqual(p[0]['concentration']*math.exp(p[0]['w_kBT']),270.,places=10)
        self.assertAlmostEqual(p[-1]['concentration']*math.exp(p[-1]['w_kBT']),1.,places=10)

    def test_detailed_balance_equal_reservoir(self):
        p=neutral_profile(350,3,894,30,30,400)
        v=[row['concentration']*math.exp(row['w_kBT']) for row in p]
        self.assertLess(max(v)-min(v),1e-10)

    def test_finite_bath_salt_conservation(self):
        for t in (0.,60.,7200.,1e6):
            l,r=finite_bath(t,1e-8,1e-4,1e-5,2e-5,270,0.)
            self.assertAlmostEqual(1e-5*l+2e-5*r,1e-5*270,places=14)

    def test_constant_P_shape_is_concave(self):
        for theta in (.001,.2,1.,3.):
            self.assertGreaterEqual(normalized_finite_bath(.5,theta),.5)

    def test_invalid_lengths_rejected(self):
        with self.assertRaises(ValueError):
            field(0,0,7,894)


if __name__=='__main__':
    unittest.main(verbosity=2)
