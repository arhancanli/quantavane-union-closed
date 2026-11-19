# Symbolic check of Boppana's proof (arXiv 2301.09664): f(x) = h(x^2) - phi x h(x) on (0,1); claim f''' = p(x)/(x(1-x^2)^2),
# p(x) = -phi x^3 - 4x^2 + 3 phi x + 2 phi - 4; roots of p; f(0) = f'(0+) = 0, f(1/phi) = f'(1/phi) = 0, f(1) = 0.
import sympy as sp
x = sp.symbols('x', positive=True)
phi = (1 + sp.sqrt(5)) / 2
h = lambda u: -u * sp.log(u) - (1 - u) * sp.log(1 - u)
f = h(x**2) - phi * x * h(x)
f3 = sp.diff(f, x, 3)
p = -phi * x**3 - 4 * x**2 + 3 * phi * x + 2 * phi - 4
print("f''' - p/(x(1-x^2)^2) simplifies to:", sp.simplify(sp.together(f3 - p / (x * (1 - x**2)**2))))
print("roots of p (numeric):", [sp.N(r_, 12) for r_ in sp.Poly(p, x).nroots()])
print("p(0) =", sp.N(p.subs(x, 0), 12), "; leading coefficient -phi < 0")
xi = 1 / phi
print("f(1/phi) =", sp.N(sp.simplify(f.subs(x, xi)), 20), " f'(1/phi) =", sp.N(sp.diff(f, x).subs(x, xi), 20))
print("f(1-) limit =", sp.limit(f, x, 1, '-'), " f(0+) =", sp.limit(f, x, 0, '+'), " f'(0+) =", sp.limit(sp.diff(f, x), x, 0, '+'))
print("f''(1/phi) =", sp.N(sp.diff(f, x, 2).subs(x, xi), 15), "(> 0 means a double root that is a local minimum)")
