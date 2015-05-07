from sympy.liealgebras.cartan_type import CartanType
from sympy.liealgebras.type_e import TypeE
from sympy.core.compatibility import range
from sympy.matrices import Matrix
from sympy.utilities.pytest import raises

def test_type_E():
    c = CartanType("E6")
    m = Matrix(6, 6, [2, 0, -1, 0, 0, 0, 0, 2, 0, -1, 0, 0,
        -1, 0, 2, -1, 0, 0, 0, -1, -1, 2, -1, 0, 0, 0, 0,
        -1, 2, -1, 0, 0, 0, 0, -1, 2])
    assert c.cartan_matrix() == m
    assert c.dimension() == 8
    assert c.simple_root(6) == [0, 0, 0, -1, 1, 0, 0, 0]
    assert c.roots() == 72
    assert c.basis() == 78
    diag = " "*8 + "2\n" + " "*8 + "0\n" + " "*8 + "|\n" + " "*8 + "|\n"
    diag += "---".join("0" for i in range(1, 6))+"\n"
    diag += "1   " + "   ".join(str(i) for i in range(3, 7))
    assert c.dynkin_diagram() == diag
    posroots = c.positive_roots()
    assert posroots[8] == [1, 0, 0, 0, 1, 0, 0, 0]

    # Test simple root
    assert c.simple_root(1) == [0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, 0.5]
    assert c.simple_root(2) == [1, 1, 0, 0, 0, 0, 0, 0]

    e7 = CartanType("E7")
    positive_roots = e7.positive_roots()
    assert positive_roots[1] == [-1, 1, 0, 0, 0, 0, 0, 0]
    assert positive_roots[10] == [1, 0, 0, 0, 0, 1, 0, 0]
    assert positive_roots[20] == [0, 0, 1, 1, 0, 0, 0, 0]
    assert positive_roots[30] == [0, 0, 0, 0, 1, 1, 0, 0]
    assert positive_roots[40] == [-1/2, -1/2, -1/2, -1/2, -1/2, 1/2, -1/2, 1/2]
    assert positive_roots[50] == [-1/2, -1/2, -1/2, -1/2, -1/2, 1/2, -1/2, 1/2]
    assert positive_roots[60] == [-1/2, -1/2, -1/2, -1/2, -1/2, 1/2, -1/2, 1/2]

    e8 = CartanType("E8")
    positive_roots = e8.positive_roots()
    assert positive_roots[1] == [-1, 1, 0, 0, 0, 0, 0, 0]
    assert positive_roots[10] == [1, 0, 0, 0, 0, 1, 0, 0]
    assert positive_roots[110] == [-1/2, -1/2, -1/2, -1/2, -1/2, 1/2, 1/2, 1/2]

    # Test roots
    assert e7.roots() == 126
    assert e8.roots() == 240

    # Test basis
    assert e7.basis() == 133
    assert e8.basis() == 248

def test_Exceptions():
    raises(ValueError, lambda: TypeE(1))
    raises(ValueError, lambda: TypeE(6).simple_root(7))
    raises(ValueError, lambda: TypeE(7).simple_root(8))

