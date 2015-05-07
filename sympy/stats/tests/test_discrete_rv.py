from sympy.stats.drv_types import (PoissonDistribution, GeometricDistribution,
        Poisson, Geometric)
from sympy.abc import x, z
from sympy import S, Sum
from sympy.stats import E, variance, density
from sympy.utilities.pytest import raises

def test_PoissonDistribution():
    l = 3
    p = PoissonDistribution(l)
    assert abs(p.cdf(10).evalf() - 1) < .001
    assert p.expectation(x, x) == l
    assert p.expectation(x**2, x) - p.expectation(x, x)**2 == l

def test_Poisson():
    l = 3
    x = Poisson('x', l)
    assert E(x) == l
    assert variance(x) == l
    assert density(x) == PoissonDistribution(l)
    assert isinstance(E(x, evaluate=False), Sum)
    assert isinstance(E(2*x, evaluate=False), Sum)

def test_GeometricDistribution():
    p = S.One / 5
    d = GeometricDistribution(p)
    assert d.expectation(x, x) == 1/p
    assert d.expectation(x**2, x) - d.expectation(x, x)**2 == (1-p)/p**2
    assert abs(d.cdf(20000).evalf() - 1) < .001
    raises(ValueError, lambda: d.check(5.0))
    raises(ValueError, lambda: d.check(-1.0))

    # Test _inverse_cdf_expression:
    prob = x
    dist = GeometricDistribution(prob)
    raises(NotImplementedError, lambda: dist._inverse_cdf_expression())

def test_geometric():
    p = S.One / 5
    x = Geometric("x", p)
    assert E(x) == 5
    assert variance(x) == 20
    assert density(x)(z) == (4/5)**(z - 1)/5
