from basics.polytools import poly, Poly

def sing(f):
	"""

	return singular locus of f over algebraic closure of f.coeff_dom

	Example:
	>>> f = poly(x**2 - y**3, mod = 3)
	>>> sing(f)
	[(0, 0)]

	"""

	validate_type(f, Poly)
	mod = f.get_modulus()
	pd_dict = dict()
	for v in f.indet_vars:
		pd_dict[v] = diff(f, v)

	rel_list = []
	count, a, p_x, sol_x = 0, [], [], []
	for f_x_ in f_x:
		f_x_ = f_x_[0]
		if not _has_roots(f_x_, x, mod):
			a.append(symbols('a_' + str(count)))
			f_x_a = f_x_.subs({x: a[count]})
			rel_list.append(f_x_a)
			p_x_ = sffpoly(f_x_, sff(f_x_a, mod))
			p_x.append(p_x_)
			sol_x.extend([{x: reduce(a[count] ** (mod ** d), p_x_.dom)} for d in range(p_x_.degree())])
			count += 1
		else:
			p_x_ = sffpoly(f_x_, sff(0, mod))
			p_x.append(p_x_)
			sol_x.extend([{x: n} for n in range(mod) if p_x_.subs({x: n}) == 0])

	count, b, p_y, sol_y = 0, [], [], []
	for f_y_ in f_y:
		f_y_ = f_y_[0]
		if not _has_roots(f_y_, y, mod):
			b.append(symbols('b_' + str(count)))
			f_y_b = f_y_.subs({y: b[count]})
			rel_list.append(f_y_b)
			p_y_ = sffpoly(f_y_, sff(f_y_b, mod))
			p_y.append(p_y_)
			sol_y.extend([{y: reduce(b[count] ** (mod ** d), p_y_.dom)} for d in range(p_y_.degree())])
			count += 1
		else:
			p_y_ = sffpoly(f_y_, sff(0, mod))
			p_y.append(p_y_)
			sol_y.extend([{y: n} for n in range(mod) if p_y_.subs({y: n}) == 0])

	sff_ = sff(rel_list, mod)
	f = sffpoly(f, sff_)
	f_x = f.diff(x)
	f_y = f.diff(y)
	sol_f = []
	for point in product(sol_x, sol_y):
		point_ = {x: point[0][x], y: point[1][y]}
		if f_x.subs(point_) == 0 and f_y.subs(point_) == 0 and f.subs(point_) == 0:
			sol_f.append(point_)

	return sol_f, sff_.as_SFF()