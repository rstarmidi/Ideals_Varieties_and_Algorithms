#myRedGroebTest2

import sympy as sym
from sympy.abc import x, y, z
from sympy import Poly
from sympy import LT, LM, LC
from sympy import lcm

def myDiv(lst, f, t, gen, dom):
    a = [0] * len(lst)
    r = 0
    p = f
    while p != 0:
        i = 0
        divisionoccured = False
        while i<len(lst) and divisionoccured == False:
            if sym.div(Poly(LT(p, order = t), gen, domain = dom), Poly(LT(lst[i], order = t), gen, domain = dom))[1] == 0 :
                a[i] = a[i] + sym.div(Poly(LT(p, order = t), gen, domain = dom), Poly(LT(lst[i], order = t), gen, domain = dom))[0]
                p = p - sym.div(Poly(LT(p, order = t), gen, domain = dom), Poly(LT(lst[i], order = t), gen, domain = dom))[0] * lst[i]
                divisionoccured = True
            else:
                i = i + 1
        if divisionoccured == False:
            r = r + LT(p, order = t)
            p = p - LT(p, order = t)
    return a, r


def calS(f, g, t):
    x_gamma = lcm(LM(f, order = t), LM(g, order = t))
    S = x_gamma/LT(f, order = t) * f - x_gamma/LT(g, order = t) * g
    return S

def myGroebner(f, t, gen, dom):
    g = f
    rolling = True
    while rolling == True:
        g_control = g.copy()
        for i in range (len(g_control)-1):
            for j in range (i+1, len(g_control)):
                S = calS(f[i], f[j], t = t)
                rem = myDiv(g_control, S, t = t, gen = gen, dom = dom)[1]
                if rem !=0:
                    g.append(Poly(rem, gen, domain = dom))
        if g == g_control:
            rolling = False

    for i in range (len(g)):
        print(f"g{i} = {sym.latex(g[i].as_expr())}")
        
    return g

def myRedGroeb(f, t, gen, dom):
    #HAVING GROEBNER BASIS
    g = myGroebner(f, t = t, gen = gen, dom = dom)

    #HAVING MONIC GROEBNER BASIS
    for i in range(len(g)):
        g[i] = sym.monic(g[i])

    #MAKING LIST OF LT(g):
    lstLT = []
    for i in range (len(g)):
        lstLT.append(Poly(LT(g[i], order = t), gen, domain = dom))

    #CHANGING THE NON REDUCED BASIS INTO ZERO
    h = g.copy()
    for i in range (len(g)):
        while h[i] != 0 :
            denum = lstLT[0:i] + lstLT[i+1:]
            rem = myDiv(denum, Poly(LT(h[i], order = t), gen, domain = dom), t = t, gen = gen, dom = dom)[1]
            if rem == 0 :
                g[i] = 0
                break
            else:
                h[i] = h[i] - LT(h[i], order = t)

    #REMOVING THE ZERO BASIS
    g = list( filter( lambda a: a != 0, g ) )
    print(f"Reduced Groebner = {g}")
    return g

f0 = Poly(x**2 + y, x, y, domain = 'QQ')
f1 = Poly(x**4 + 2*x**2*y + y**2 + 3, x, y, domain = 'QQ')

f = [f0, f1]

myRedGroeb(f, t = "lex", gen = (x, y), dom = 'QQ')
