from sympy.core.numbers import ilcm, Integer
from sympy.core.expr import Expr
from sympy.core.function import expand
from sympy.ntheory.primetest import isprime
from sympy.polys.polytools import LC, Poly, poly

from exprtools import reduce_coefficients

import random

class FF:
    """
        represents a splitted finite field of some polynomials
        over a finite field of which modulus is 'mod'.
    """

    def __init__(self, rel, mod):
        """
        Instance variables:
            * rel_list: list of relational equations
                        of which element is a dict {'var': variable, 'rep': equation}
            * mod: characteristic number
            * var_list: list of variables of relational equations
            * exdeg: extension degree
            * num: number of elements
            * gens: list of generators as vector space

        Example:
            a_1 = symbols('a_1') 
            rel = {'var': a_1, 'rep': a_1 ** 2 - 3, 'deg': 2, 'uni': True}

            self.rel_list = [rel_1, rel_2, ..., rel_n]
            self.mod = p
            self.var_list = [a_1, a_2, ..., a_m]
            self.exdeg = e
            self.num = p ** e
            self.gens = [g_1, g_2, ..., g_(p ** e)]
        """
        self.rel_list = []
        _dom = 'FF(' + str(mod) + ')'
        if rel == 0:
            self.is_prime = True
        elif isinstance(rel, Expr):
            self.is_prime = False
            if not LC(rel.as_poly()) == 1:
                rel = poly(rel.as_expr() * pow(LC(rel.as_poly()), mod - 2, mod), domain=_dom).as_expr()
            if len(poly(rel).gens) == 1:
                self.rel_list.append({'var': poly(rel).gens[0],
                                      'rep': rel.as_expr(), 
                                      'deg': poly(rel).degree(), 
                                      'is_uni': True})
            else:
                self.rel_list.append({'var': poly(rel).gens[0], 
                                      'rep': rel.as_expr(), 
                                      'deg': poly(rel).degree(), 
                                      'is_uni': False})
        elif isinstance(rel, list):
            self.is_prime = False
            if rel == [] or rel == [0] * len(rel):
                self.is_prime = True
            else:
                for _p in rel:
                    if _p == 0:
                        continue
                    if not LC(_p.as_poly()) == 1:
                        _p = poly(_p * pow(LC(_p.as_poly()), mod - 2, mod), domain=_dom).as_expr()
                    if len(poly(_p).gens) == 1:
                        self.rel_list.append({'var': poly(_p).gens[0], 
                                              'rep': _p.as_expr(), 
                                              'deg': poly(_p).degree(), 
                                              'is_uni': True})
                    else:
                        self.rel_list.append({'var': poly(_p).gens[0], 
                                              'rep': _p.as_expr(), 
                                              'deg': poly(_p).degree(), 
                                              'is_uni': False})
        else: 
            raise ValueError("first argument needs to be 0, an Expr instance or list.")

        if isprime(mod):
            self.mod = mod
        else:
            raise ValueError("modulus needs to be a prime number")

        if self.is_prime:
            self.var_list = []
            self.exdeg = 1
            self.num = self.mod
            self.gens = [1]
        else:
            _var_list, _degs, _gens = [], [], []
            for _rel in self.rel_list:
                _degs.append(_rel['deg'])
                if not _rel['var'] in _var_list and _rel['is_uni']:
                    _var_list.append(_rel['var'])
            self.var_list = _var_list
            self.exdeg = ilcm(1, *_degs)
            self.num = mod ** self.exdeg
            """
                Note:
                Below codes include errors.
                When degree of two extension variables are each 2 and 4 then returns 8 generators,
                however acutually there exist only 4 generators.
            """
            _expr = 1
            for _var in self.var_list:
                _deg = next(item for item in self.rel_list if item['var'] == _var and item['is_uni'])['deg']
                _expr_1 = 0
                for i in range(_deg):
                    _expr_1 += _var ** i
                _expr *= _expr_1
            _tuple = expand(_expr).as_coeff_add()
            self.gens = tuple([_tuple[0]] + [item for item in _tuple[1]])

    def __str__(self):
        return self.as_SFF()

    def __repr__(self):
        return self.as_SFF()

    def as_FF(self):
        return "FF(%s**%s)" % (self.mod, self.exdeg)

    def as_sympy_FF(self):
        return "FF(%s)" % self.mod

    def as_SFF(self):
        if self.is_prime:
            return self.as_sympy_FF()
        else:
            return "FiniteField(%s**%s) splitting %s" % (self.mod, self.exdeg, [rel['rep'] for rel in self.rel_list if rel['is_uni']])

    def rel_deg(self, **args):
        if not args:
            return poly(self.rel_list[0]['rep']).degree()
        elif 'var' in args:
            for rel in self.rel_list:
                if args['var'] == rel['var']:
                    return poly(rel['rep']).degree()
            raise ValueError("doesn't have the variable")
        elif 'index' in args:
            return poly(self.rel_list[args['index']]['rep']).degree()
        else:
            raise ValueError("rel_deg() doesn't have the argument option")

    def element(self, i):
        if i > self.num:
            raise ValueError
        _rep = 0
        for d in range(self.exdeg)[::-1]:
            c = i // (self.mod ** d)
            if c > self.mod / 2:
                c -= self.mod
            _rep += c * self.gens[d]
            i %= self.mod ** d
        return _rep

    def elements_iter(self):
        for i in range(self.num):
            yield self.element(i)

    def elements_list(self):
        return [f for f in self.elements_iter()]

    def elements_with_args_iter(self, args):
        for i in range(self.num):
            yield (self.element(i), args)

    def rand(self):
        return self.element(random.randint(0, self.num - 1))

    def point(self, n, i):
        if i > self.num ** n:
            raise ValueError
        _list = []
        for d in range(n)[::-1]:
            _list.insert(0, self.element(int(i / (self.num ** d))))
            i %= self.num ** d
        return tuple(_list)

    def point_as_dict(self, i, gens):
        n = len(gens)
        return dict((gens[d], self.point(n, i)[d]) for d in range(n))

    def points_iter(self, n):
        for i in range(self.num ** n):
            yield self.point(n, i)

    def points_as_dict_iter(self, gens):
        n = len(gens)
        for i in range(self.num ** n):
            yield self.point_as_dict(i, gens)

    def points_as_dict_with_poly_and_queue_iter(self, gens, poly, queue):
        n = len(gens)
        for i in range(self.num ** n):
            yield (poly, self.point_as_dict(i, gens), queue)

    def extend(self, rep):
        rel_ = [rel['rep'] for rel in self.rel_list]
        rel_.append(rep)
        return sff(rel_, self.mod)

def sff(rel, mod):
    return SFF(rel, mod)