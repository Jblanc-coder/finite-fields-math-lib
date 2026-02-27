"""PRÁCTICA JUAN BLANCO GUILLÉN"""
# Se deben implementar 4 clases: cuerpo_fp, anillo_fp_x, cuerpo_fq, anillo_fq_x
# respetando el esqueleto de abajo. Se pueden añadir funciones auxiliares, pero
# las indicadas deben estar. Los elementos de cada uno de estos cuerpos/anillos
# son objetos opacos.
import random
class cuerpo_fp:
    def __init__(self, p): # construye el cuerpo de p elementos Fp = Z/pZ
        if not isinstance(p, int) or p <= 1:
            raise ValueError("p debe ser un entero > 1 (idealmente primo).")
        self.p = p
    def cero(self):        # devuelve el elemento 0
        return 0
    def uno(self):         # devuelve el elemento 1
        return 1
    def elem_de_int(self, n): # fabrica el elemento dado por la clase de n
        if not isinstance(n, int):
            raise TypeError("n debe ser un entero.")
        return n % self.p
    def elem_de_str(self, s): # fabrica el elemento a partir de un string (parser)
        if not isinstance(s, str):
            raise TypeError("s debe ser un str.")
        return int(s) % self.p
    def conv_a_int(self, a):  # devuelve un entero entre 0 y p-1
        return int(a) % self.p
    def conv_a_str(self, a):  # pretty-printer
        return str(self.conv_a_int(a))
    def suma(self, a, b):     # a+b
        return (a + b) % self.p
    def inv_adit(self, a):    # -a
        return (-a) % self.p
    def neg(self, a):
        return self.inv_adit(a)
    def mult(self, a, b):     # a*b
        return (a*b) % self.p
    def pot(self, a, k):      # a^k (k entero)
        if not isinstance(k, int):
            raise TypeError("k debe ser un entero.")
        if k < 0:
            a = self.inv_mult(a)
            k = -k
        if k == 0:
            return 1
        b = self.pot(a, k//2)
        b = self.mult(b, b)
        if k % 2 == 1:
            b = self.mult(b, a)
        return b
    def inv_mult(self, a):    # a^(-1)
        if a%self.p == 0:
            raise ValueError(" 0 no es invertible")
        return self.pot(a, self.p-2) #Por el Pequeño Teorema de Fermat
    def es_cero(self, a):     # a == 0
        return a%self.p == 0
    def es_uno(self, a):      # a == 1
        return a%self.p == 1
    def es_igual(self, a, b): # a == b
        return (a-b)%self.p == 0
    def aleatorio(self):      # fabrica un elemento aleatorio con prob uniforme
        return random.randint(0,self.p-1)
    def tabla_suma(self):     # devuelve la matriz de pxp (lista de listas) de la suma
        p = self.p
        return [[(i + j) % p for j in range(p)] for i in range(p)]
    def tabla_mult(self):     # devuelve la matriz de pxp (lista de listas) de la mult
        p = self.p
        return [[(i * j) % p for j in range(p)] for i in range(p)]
    def tabla_inv_adit(self): # devuelve una lista de long p con los inv_adit
        p = self.p
        return [(-i) % p for i in range(p)]
    def tabla_inv_mult(self): # devuelve una lista de long p con los inv_mult (en el índice 0 pone un '*')
        p = self.p
        inv = [''*''] * p
        for i in range(1, p):
            inv[i] = self.inv_mult(i)
        return inv
    def cuadrado_latino(self, a): # cuadrado latino a*i+j (con a != 0)
        if self.es_cero(a):
            raise ValueError("Cuadrado latino no está definido para a=0")
        return [[self.suma(j, self.mult(a,i)) for j in range(self.p)] for i in range(self.p)]
    
def factorizar_entero(n):
    fact = []
    p=2
    while n > 1:
        e = 0
        while n%p == 0:
            n //= p
            e += 1
        if e != 0:
            fact.append((p,e))
        p+=1
    return fact

class anillo_fp_x:
    def __init__(self, fp, var='x'): # construye el anillo Fp[var]
        self.fp = fp
        self.var = var
        self.campo = fp
    def reducir(self, a):  #quita ceros a la izquierda y devuelve tupla
        a = tuple(a)
        i = len(a)
        while i > 0 and self.fp.es_cero(a[i-1]):
            i -= 1
        return tuple(a[:i])   
    def cero(self):                 # 0
        return ()
    def uno(self):                  # 1
        return (1,)
    def elem_de_tuple(self, a): # fabrica un polinomio a partir de la tupla (a0, a1, ...)
        return self.reducir(tuple(self.fp.elem_de_int(n) for n in a))
    def elem_de_int(self, a):  # fabrica un polinomio a partir de los dígitos de a en base p
        if a == 0:
            return self.cero()
        p = self.fp.p
        coeffs = []
        n = int(a)
        while n > 0:
            coeffs.append(n % p)
            n //= p
        return self.reducir(tuple(coeffs))
    def elem_de_str(self, s):       # fabrica un polinomio a partir de un string (parser)
        s = s.strip()
        if s == '':
            return self.cero()
        if ',' in s:
            parts = [p.strip() for p in s.split(',') if p.strip()!='']
        else:
            parts = [p for p in s.split() if p!='']
        nums = [int(t) for t in parts]
        return self.elem_de_tuple(nums)
    def conv_a_tuple(self, a):      # devuelve la tupla de coeficientes
        return tuple(a)
    def conv_a_int(self, a):        # devuelve el entero correspondiente al polinomio
        a = self.reducir(a)
        p = self.fp.p
        res = 0
        mul = 1
        for coeff in a:
            res += int(coeff) * mul
            mul *= p
        return res
    def conv_a_str(self, a):        # pretty-printer
        a = self.reducir(a)
        if len(a) == 0:
            return '0'
        terms = []
        for i, c in enumerate(a):
            c = int(c)
            if c == 0:
                continue
            if i == 0:
                terms.append(str(c))
            elif i == 1:
                if c == 1:
                    terms.append(self.var)
                else:
                    terms.append(f"{c}*{self.var}")
            else:
                if c == 1:
                    terms.append(f"{self.var}^{i}")
                else:
                    terms.append(f"{c}*{self.var}^{i}")
        return ' + '.join(reversed(terms)) if terms else '0'
    def suma(self, a, b):           # a+b
        la = len(a)
        lb = len(b)
        n = max(la, lb)
        res = [self.fp.cero()]*n 
        for i in range(n):
            ca = a[i] if i < la else 0
            cb = b[i] if i < lb else 0
            res[i] = self.fp.suma(ca, cb)
        return self.reducir(tuple(res))    
    def inv_adit(self, a):          # -a
        return self.reducir(tuple(self.fp.inv_adit(c) for c in a))
    def mult(self, a, b):           # a*b
        if self.es_cero(a) or self.es_cero(b):
            return self.cero()
        na = len(a)
        nb = len(b)
        c = [self.fp.cero()] * (na+nb-1)
        for i in range(na):
            for j in range(nb):
                c[i+j] = self.fp.suma(c[i+j], self.fp.mult(a[i], b[j]))
        return self.reducir(tuple(c))
    def mult_por_escalar(self, a, e): # a*e (con e en Z/pZ)
        return self.reducir(tuple(self.fp.mult(x, e) for x in a))
    def divmod(self, a, b):        # devuelve q,r tales que a=bq+r y deg(r)<deg(b)
        if self.es_cero(b):
            raise ZeroDivisionError("no se puede dividir por 0")
        q = [0]* max(len(a) - len(b)+1,1)
        r = a[:] 
        while len(r) >= len(b) and r != [0]:
            coef = self.fp.mult(r[-1],self.fp.inv_mult(b[-1]))
            deg = len(r) - len(b)
            coci = [0]*deg + [coef]
            q = self.suma(q, coci)
            r = self.suma(r, self.inv_adit(self.mult(b,coci)))
        return (self.reducir(q), self.reducir(r))
    def div(self, a, b):           # q
        q, r = self.divmod(a, b)
        return q
    def mod(self, a, b):           # r
        q, r = self.divmod(a, b)
        return r
    def grado(self, a):            # deg(a)
        a = self.reducir(a)
        return len(a) - 1 if len(a) > 0 else -1
    def monico(self, a): #pol mónico proporcional
        if self.es_cero(a):
            return a
        a_norm = self.elem_de_tuple(a)
        coef = self.fp.inv_mult(a_norm[-1])
        a_mon = self.mult_por_escalar(a_norm,coef)
        return a_mon
    def gcd(self, a, b):           
        a = self.reducir(a)
        b = self.reducir(b)
        while not self.es_cero(b):
            r = self.mod(a, b)
            a, b = b, r
        if self.es_cero(a):
            return self.cero()
        # normalizar a mónico
        return self.monico(a)
    def gcd_ext(self, a, b): # devuelve g,x,y tales que g=ax+by, g=gcd(a,b) mónico
       a = self.reducir(a)
       b = self.reducir(b)
       x0, x1 = self.uno(), self.cero()
       y0, y1 = self.cero(), self.uno()
       while not self.es_cero(b):
           q, r = self.divmod(a, b)
           a, b = b, r
           x0, x1 = x1, self.suma(x0, self.inv_adit(self.mult(q, x1)))
           y0, y1 = y1, self.suma(y0, self.inv_adit(self.mult(q, y1)))
       if self.es_cero(a):
           return self.cero(), self.cero(), self.cero()
       inv_lc = self.fp.inv_mult(a[-1])
       g = self.mult_por_escalar(a, inv_lc)
       x = self.mult_por_escalar(x0, inv_lc)
       y = self.mult_por_escalar(y0, inv_lc)
       return g, x, y
    def inv_mod(self, a, b):       # devuelve x tal que ax = 1 mod b
        g, x, y = self.gcd_ext(a, b)
        if not self.es_uno(g):
            raise ValueError("No existe inverso modular; gcd != 1")
        return self.mod(x, b)
    def pot_mod(self, a, k, b):    # a^k mod b
        if k < 0:
            a = self.inv_mod(a, b)
            k = -k
        res = self.uno()
        base = self.mod(a, b)
        while k > 0:
            if k & 1:
                res = self.mod(self.mult(res, base), b)
            base = self.mod(self.mult(base, base), b)
            k >>= 1
        return res
    def es_cero(self, a):          # a == 0
        return len(self.reducir(a)) == 0
    def es_uno(self, a):           # a == 1
        return a == (1,)
    def es_igual(self, a, b):      # a == b
        return self.reducir(a) == self.reducir(b)
    def es_irreducible(self, f):   # test de irreducibilidad de Rabin
        if self.es_cero(f):
            raise ValueError("El polinomio cero no es irreducible")
        n = self.grado(f)
        if n <= 0:
            return True  # polinomios constantes (grado 0)  son irreducibles
        if f[-1] != 1:
            inv_lider = self.fp.inv_mult(f[-1])
            f_monico = self.mult_por_escalar(f, inv_lider)
        else:
            f_monico = f
        p = self.fp.p
        factores = factorizar_entero(n)
        divisores_primos = [d for d, exp in factores]  
        x = (0, 1)  
        x_pn = self.pot_mod(x, p**n, f_monico)
        if not self.es_igual(x_pn, x):
            return False
        for q in divisores_primos:
            k = n // q
            x_pk = self.pot_mod(x, p**k, f_monico)
            x_pk_minus_x = self.suma(x_pk, self.inv_adit(x))
            g = self.gcd(x_pk_minus_x, f_monico)
            if not self.es_uno(g):
                return False
        return True

class cuerpo_fq:
    def __init__(self, fp, g, var='a'): # construye el cuerpo Fp[var]/<g(var)>
        self.fp = fp
        self.g = g
        self.anillo_base = anillo_fp_x(fp, var=var)  
        if not self.anillo_base.es_irreducible(g):
             raise ValueError("El polinomio g debe ser irreducible para formar un cuerpo")
        self.n = self.anillo_base.grado(g)  
        self.q = self.fp.p ** self.n        
        self.var = var
    def reducir(self, a):  #pol a mod g
        return self.anillo_base.mod(a, self.g)
    def cero(self):                       # 0
        return self.anillo_base.cero()
    def uno(self):                        # 1
        return self.anillo_base.uno()
    def elem_de_tuple(self, a):           # fabrica elemento a partir de tupla de coeficientes
        # a es una tupla de elementos de Fp
        return self.reducir(self.anillo_base.elem_de_tuple(a))
    def elem_de_int(self, a):             # fabrica elemento a partir de entero
        # Convertimos el entero a un polinomio en Fp[a] y luego reducimos
        return self.reducir(self.anillo_base.elem_de_int(a))
    def elem_de_str(self, s):             # fabrica elemento parseando string
        # Usamos el parser de anillo_fp_x y luego reducimos
        return self.reducir(self.anillo_base.elem_de_str(s))
    def conv_a_tuple(self, a):            # devuelve tupla de coeficientes sin ceros "extra"
        # La representación del elemento *es* ya la tupla reducida
        return self.anillo_base.conv_a_tuple(a)
    def conv_a_int(self, a):              # devuelve el entero correspondiente
        # Si el grado es muy alto, esto puede ser un número ENORME.
        return self.anillo_base.conv_a_int(a)
    def suma(self, a, b):                 # a+b
        return self.reducir(self.anillo_base.suma(a, b))
    def inv_adit(self, a):                # -a
        return self.reducir(self.anillo_base.inv_adit(a))
    def neg(self, a):
        return self.inv_adit(a)
    def mult(self, a, b):                 # a*b
        # Multiplicación en el anillo y luego reducción módulo g
        return self.reducir(self.anillo_base.mult(a, b))
    def pot(self, a, k):                  # a^k (k entero)
        if k == 0:
            return self.uno()
        # ExponenciaciÃ³n modular a^k mod g
        return self.anillo_base.pot_mod(a, k, self.g)
    def inv_mult(self, a):                # a^(-1)
        # Inverso modular en el anillo base: a^(-1) mod g
        if self.es_cero(a):
            raise ValueError("0 no es invertible")
        return self.anillo_base.inv_mod(a, self.g)
    def es_cero(self, a):                 # a == 0
        return self.anillo_base.es_cero(a)
    def es_uno(self, a):                  # a == 1
        return self.anillo_base.es_uno(a)
    def es_igual(self, a, b):             # a == b
        return self.anillo_base.es_igual(a, b)
    def aleatorio(self):                  # devuelve un elemento aleatorio con prob uniforme
        coeffs = [self.fp.aleatorio() for _ in range(self.n)]
        return self.anillo_base.reducir(tuple(coeffs))
    def _elem_de_idx(self, i):
        a = self.anillo_base.elem_de_int(i)
        return a
    def tabla_suma(self):  # matriz de qxq correspondiente a la suma 
        q = self.q
        tabla = []
        for i in range(q):
            row = []
            a = self._elem_de_idx(i)
            for j in range(q):
                b = self._elem_de_idx(j)
                s = self.suma(a, b)
                row.append(self.conv_a_int(s))
            tabla.append(row)
        return tabla
    def tabla_mult(self): 
        q = self.q
        tabla = []
        for i in range(q):
            row = []
            a = self.elem_de_int(i)
            for j in range(q):
                b = self.elem_de_int(j)
                m = self.mult(a, b)
                row.append(self.conv_a_int(m))
            tabla.append(row)
        return tabla
    def tabla_inv_adit(self):  
        q = self.q
        inv = []
        for i in range(q):
            a = self._elem_de_idx(i)
            inv_a = self.inv_adit(a)
            inv.append(self.conv_a_int(inv_a))
        return inv     
    def tabla_inv_mult(self):
        q = self.q
        inv = []
        for i in range(q):
            a = self._elem_de_idx(i)
            if self.es_cero(a):
                inv.append('*')
            else:
                inv_a = self.inv_mult(a)
                inv.append(self.conv_a_int(inv_a))
        return inv
    def cuadrado_latino(self, a):  # cuadrado latino para a != 0 
        if self.es_cero(a):
            raise ValueError("Cuadrado latino no está definido para a=0")
        q = self.q
        tabla = []
        for i in range(q):
            row = []
            elem_i = self._elem_de_idx(i)
            for j in range(q):
                elem_j = self._elem_de_idx(j)
                # a*i + j
                res = self.suma(self.mult(a, elem_i), elem_j)
                row.append(self.conv_a_int(res))
            tabla.append(row)
        return tabla

class anillo_fq_x:
    def __init__(self, fq, var='x'): # Fq[var], var debe ser distinta que la de fq
        self.fq = fq
        self.var = var
        self.campo = fq
    def reducir(self, a):
        a = tuple(a)
        i = len(a)
        while i > 0 and self.fq.es_cero(a[i-1]):
            i -= 1
        return tuple(a[:i])
    def cero(self):                     # 0
        return (self.fq.cero(),)
    def uno(self):                      # 1
        return (self.fq.uno(),)
    def elem_de_tuple(self, a):
        return self.reducir(tuple(a))
    def elem_de_int(self, a):
        if a == 0:
            return self.cero()
        q = self.fq.q
        coeffs = []
        n = int(a)
        while n > 0:
            digito = n % q
            elem_fq = self.fq.elem_de_int(digito)
            coeffs.append(elem_fq)
            n //= q
        return self.reducir(tuple(coeffs))
    def elem_de_str(self, s):
        s = s.strip()
        if s == '':
            return self.cero()
        if ',' in s:
            parts = [p.strip() for p in s.split(',') if p.strip() != '']
        else:
            parts = [p for p in s.split() if p != '']
        coeffs = []
        for part in parts:
            try:
                elem = self.fq.elem_de_str(part)
                coeffs.append(elem)
            except:
                try:
                    elem = self.fq.elem_de_int(int(part))
                    coeffs.append(elem)
                except:
                    coeffs.append(self.fq.cero())
        return self.elem_de_tuple(coeffs)
    def conv_a_tuple(self, a):
        return self.reducir(a)
    def conv_a_int(self, a):
        a = self.reducir(a)
        q = self.fq.q
        res = 0
        mul = 1
        for coeff in a:
            coeff_int = self.fq.conv_a_int(coeff)
            res += coeff_int * mul
            mul *= q
        return res
    def conv_a_str(self, a):
        a = self.reducir(a)
        if len(a) == 0:
            return '0'
        terms = []
        for i, c in enumerate(a):
            if self.fq.es_cero(c):
                continue
            coeff_str = self.fq.conv_a_str(c)
            if i == 0:
                terms.append(coeff_str)
            elif i == 1:
                if coeff_str == '1':
                    terms.append(self.var)
                else:
                    terms.append(f"{coeff_str}{self.var}")
            else:
                if coeff_str == '1':
                    terms.append(f"{self.var}^{i}")
                else:
                    terms.append(f"{coeff_str}{self.var}^{i}")
        return ' + '.join(reversed(terms)) if terms else '0'
    def suma(self, a, b):
        la = len(a)
        lb = len(b)
        n = max(la, lb)
        res = []
        for i in range(n):
            ca = a[i] if i < la else self.fq.cero()
            cb = b[i] if i < lb else self.fq.cero()
            res.append(self.fq.suma(ca, cb))
        return self.reducir(tuple(res))
    def inv_adit(self, a):
        return self.reducir(tuple(self.fq.inv_adit(c) for c in a))
    def mult(self, a, b):
        if self.es_cero(a) or self.es_cero(b):
            return self.cero()
        a = self.reducir(a)
        b = self.reducir(b)
        na = len(a)
        nb = len(b)
        c = [self.fq.cero()] * (na + nb - 1)
        for i in range(na):
            for j in range(nb):
                producto = self.fq.mult(a[i], b[j])
                c[i+j] = self.fq.suma(c[i+j], producto)
        return self.reducir(tuple(c))
    def mult_por_escalar(self, a, e):
        return self.reducir(tuple(self.fq.mult(x, e) for x in a))
    def divmod(self, a, b):
        if self.es_cero(b):
            raise ValueError("No se puede dividir por 0")
        a = self.reducir(a)
        b = self.reducir(b)
        if self.es_cero(a) or self.grado(a) < self.grado(b):
            return self.cero(), a
        q_coeffs = [self.fq.cero()] * (self.grado(a) - self.grado(b) + 1)
        r = list(a)  
        inv_lead_b = self.fq.inv_mult(b[-1])
        for i in range(self.grado(a) - self.grado(b), -1, -1):
            if self.fq.es_cero(r[i + len(b) - 1]):
                continue
            factor = self.fq.mult(r[i + len(b) - 1], inv_lead_b)
            q_coeffs[i] = factor
            for j in range(len(b)):
                if i + j < len(r):
                    termino = self.fq.mult(factor, b[j])
                    r[i + j] = self.fq.suma(r[i + j], self.fq.inv_adit(termino))  
        return self.reducir(tuple(q_coeffs)), self.reducir(tuple(r))
    def div(self, a, b):
        q, r = self.divmod(a, b)
        return q
    def mod(self, a, b):
        q, r = self.divmod(a, b)
        return r
    def grado(self, a):
        a = self.reducir(a)
        return len(a) - 1 if len(a) > 0 else -1
    def gcd(self, a, b):
        a = self.reducir(a)
        b = self.reducir(b)
        # Casos especiales
        if self.es_cero(a) and self.es_cero(b):
            return self.cero()
        if self.es_cero(a):
            return self.mult_por_escalar(b, self.fq.inv_mult(b[-1])) if not self.es_cero(b) else self.cero()
        if self.es_cero(b):
            return self.mult_por_escalar(a, self.fq.inv_mult(a[-1])) if not self.es_cero(a) else self.cero()
        # Algoritmo de Euclides
        while not self.es_cero(b):
            r = self.mod(a, b)
            a, b = b, r
        # Hacer mónico
        if self.es_cero(a):
            return self.cero()
        inv_lc = self.fq.inv_mult(a[-1])
        return self.mult_por_escalar(a, inv_lc)
    def gcd_ext(self, a, b):
        a = self.reducir(a)
        b = self.reducir(b)
        if self.es_cero(b):
            if self.es_cero(a):
                return self.cero(), self.cero(), self.cero()
            inv_lc = self.fq.inv_mult(a[-1])
            return (self.mult_por_escalar(a, inv_lc), 
                    self.elem_de_tuple([inv_lc]), 
                    self.cero())
        x0, x1 = self.uno(), self.cero()
        y0, y1 = self.cero(), self.uno()
        while not self.es_cero(b):
            q, r = self.divmod(a, b)
            a, b = b, r
            x0, x1 = x1, self.suma(x0, self.inv_adit(self.mult(q, x1)))
            y0, y1 = y1, self.suma(y0, self.inv_adit(self.mult(q, y1)))
        if self.es_cero(a):
            return self.cero(), self.cero(), self.cero()
        # Normalizar a mónico
        inv_lc = self.fq.inv_mult(a[-1])
        g = self.mult_por_escalar(a, inv_lc)
        x = self.mult_por_escalar(x0, inv_lc)
        y = self.mult_por_escalar(y0, inv_lc)
        return g, x, y
    def inv_mod(self, a, b):
        g, x, y = self.gcd_ext(a, b)
        if not self.es_uno(g):
            raise ValueError("No existe inverso modular; gcd != 1")
        return self.mod(x, b)
    def pot_mod(self, a, k, b):
        if k < 0:
            a = self.inv_mod(a, b)
            k = -k
        if k == 0:
            return self.uno()        
        res = self.uno()
        base = self.mod(a, b)
        while k > 0:
            if k & 1:
                res = self.mod(self.mult(res, base), b)
            base = self.mod(self.mult(base, base), b)
            k >>= 1     
        return res
    def es_cero(self, a):
        return len(self.reducir(a)) == 0
    def es_uno(self, a):
        a_reducido = self.reducir(a)
        return len(a_reducido) == 1 and self.fq.es_uno(a_reducido[0])
    def es_igual(self, a, b):
        return self.reducir(a) == self.reducir(b)
    def derivada(self, f):
        f = self.reducir(f)
        if self.grado(f) <= 0:
            return self.cero()

        deriv_coeffs = []
        p = self.fq.fp.p  # característica de F_q

        for i in range(1, len(f)):
            coef_i = f[i]
            k = i % p
            if k == 0:
                # el factor i es múltiplo de p -> coeficiente derivado 0 en característica p
                deriv_coeffs.append(self.fq.cero())
            else:
                # metemos k como ESCALAR en F_q (constante)
                mult_coef = self.fq.elem_de_int(k)
                deriv_coeff = self.fq.mult(coef_i, mult_coef)
                deriv_coeffs.append(deriv_coeff)
        return self.reducir(tuple(deriv_coeffs))
        return self.elem_de_tuple(tuple(deriv_coeffs))
    def es_irreducible(self, f):
        if self.es_cero(f):
            raise ValueError("El polinomio cero no es irreducible")
        n = self.grado(f)
        if n <= 0:
            return False  # Polinomios constantes no son irreducibles (excepto tal vez por definiciÃ³n)
        if n == 1:
            return True   # Todos los polinomios lineales son irreducibles
        # Hacer f mónico
        if not self.fq.es_uno(f[-1]):
            inv_lider = self.fq.inv_mult(f[-1])
            f_monico = self.mult_por_escalar(f, inv_lider)
        else:
            f_monico = f
        q = self.fq.q
        # Obtener divisores primos de n
        def div_primos(n):
            factores = factorizar_entero(n)
            return [d for d, _ in factores]
        divisores = div_primos(n)
        x = self.elem_de_tuple((self.fq.cero(), self.fq.uno()))  # polinomio x
        try:
            x_qn = self.pot_mod(x, q**n, f_monico)
            if not self.es_igual(x_qn, x):
                return False
        except:
            return False
        for d in divisores:
            k = n // d
            try:
                x_qk = self.pot_mod(x, q**k, f_monico)
                x_qk_minus_x = self.suma(x_qk, self.inv_adit(x))
                g = self.gcd(x_qk_minus_x, f_monico)
                if not self.es_uno(g):
                    return False
            except:
                return False
        return True
    def monico(self, a):
        # Si es cero, devolvemos cero
        if self.es_cero(a):
            return a
        
        # El coeficiente líder es el último elemento de la tupla reducida
        lc = a[-1]
        
        # Si ya es 1, devolvemos el polinomio tal cual
        if self.fq.es_uno(lc):
            return self.reducir(a)
            
        # Calculamos el inverso del coeficiente líder en Fq
        inv_lc = self.fq.inv_mult(lc)
        
        # Multiplicamos todo el polinomio por ese inverso para hacer el líder 1
        return self.mult_por_escalar(a, inv_lc)