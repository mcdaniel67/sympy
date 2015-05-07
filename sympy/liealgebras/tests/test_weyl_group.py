from sympy.liealgebras.weyl_group import WeylGroup
from sympy.matrices import Matrix
from mpmath import mpf

def test_weyl_group_a():
    # A WeylGroups
    c = WeylGroup("A3")
    assert c.matrix_form('r1*r2') == Matrix([[0, 0, 1, 0], [1, 0, 0, 0],
        [0, 1, 0, 0], [0, 0, 0, 1]])
    assert c.generators() == ['r1', 'r2', 'r3']
    assert c.group_order() == 24.0
    assert c.group_name() == "S4: the symmetric group acting on 4 elements."
    assert c.coxeter_diagram() == "0---0---0\n1   2   3"
    assert c.element_order('r1*r2*r3') == 4
    assert c.element_order('r1*r3*r2*r3') == 3

def test_weyl_group_b():
    # B WeylGroups
    d = WeylGroup("B5")
    assert d.group_order() == 3840
    assert d.element_order('r1*r2*r4*r5') == 12
    assert d.matrix_form('r2*r3') ==  Matrix([[0, 0, 1, 0, 0], [1, 0, 0, 0, 0],
        [0, 1, 0, 0, 0], [0, 0, 0, 1, 0], [0, 0, 0, 0, 1]])
    assert d.element_order('r1*r2*r1*r3*r5') == 6
    assert d.group_name() == "The hyperoctahedral group acting on 10 elements."
    assert d.coxeter_diagram() == '0---0---0---0===0\n1   2   3   4   5'

def test_weyl_group_d():
    # D WeylGroups
    e = WeylGroup("D5")
    assert e.element_order('r2*r3*r5') == 4
    assert e.matrix_form('r2*r3*r5') == Matrix([[1, 0, 0, 0, 0], [0, 0, 0, 0, -1],
        [0, 1, 0, 0, 0], [0, 0, 1, 0, 0], [0, 0, 0, -1, 0]])
    assert e.group_order() == mpf('1920.0')
    assert e.group_name() == "The symmetry group of the 5-dimensional demihypercube."

def test_weyl_group_g():
    # G WeylGroups
    f = WeylGroup("G2")
    assert f.element_order('r1*r2*r1*r2') == 3
    assert f.element_order('r2*r1*r1*r2') == 1
    assert f.matrix_form('r1*r2*r1*r2') == Matrix([[0, 1, 0], [0, 0, 1], [1, 0, 0]])
    assert f.group_order() == 12
    assert f.group_name() == "D6, the dihedral group of order 12, and symmetry group of the hexagon."
    assert f.coxeter_diagram() == '0\xe2\x89\xa1\xe2\x89\xa1\xe2\x89\xa10\n1   2'
    assert f.element_order('r1*r2*r3*r5*r1') == 2
    assert f.element_order('r1') == 2

def test_weyl_group_f():
    # F WeylGroups
    g = WeylGroup("F4")
    assert g.matrix_form('r2*r3') == Matrix([[1, 0, 0, 0], [0, 1, 0, 0],
        [0, 0, 0, -1], [0, 0, 1, 0]])
    assert g.element_order('r2*r3') == 4
    assert g.group_order() == 1152
    assert g.group_name() == "The symmetry group of the 24-cell, or icositetrachoron."
    assert g.coxeter_diagram() == '0---0===0---0\n1   2   3   4'
    assert g.element_order('r1*r5') == 2

def test_weyl_group_e():
    # E WeylGroups
    we6 = WeylGroup("E6")
    assert we6.group_order() == 51840
    assert we6.group_name() == "The symmetry group of the 6-polytope."

    we7 = WeylGroup("E7")
    assert we7.group_order() == 2903040
    assert we7.group_name() == "The symmetry group of the 7-polytope."
    assert we7.element_order('r56') == 2
    mtrx_form_r1r2 = Matrix([
        [-1/4, -3/4,  1/4,  1/4,  1/4,  1/4,  1/4, -1/4],
        [-3/4, -1/4, -1/4, -1/4, -1/4, -1/4,  1/4, -1/4],
        [1/4, -1/4,  3/4, -1/4, -1/4, -1/4, -1/4,  1/4],
        [1/4, -1/4, -1/4,  3/4, -1/4, -1/4, -1/4,  1/4],
        [1/4, -1/4, -1/4, -1/4,  3/4, -1/4, -1/4,  1/4],
        [1/4, -1/4, -1/4, -1/4, -1/4,  3/4, -1/4,  1/4],
        [1/4, -1/4, -1/4, -1/4, -1/4, -1/4, -3/4,  1/4],
        [1/4, -1/4, -1/4, -1/4, -1/4, -1/4, -1/4,  3/4]])
    assert we7.matrix_form('r1*r2') == mtrx_form_r1r2

    we8 = WeylGroup("E8")
    assert we8.group_order() == 696729600
    assert we8.group_name() == "The symmetry group of the 8-polytope."