# para esta actividad necesitamos la anterior (cuerpos_finitos.py)
# completa, pero quitando las funciones de factorización en anillo_fp_x
# y en anillo_fq_x, ya que las implementaremos aquí (no se olviden de
# quitarlas)
import cuerpos_finitos as cf

# square-free factorization
# input: fpx --> anillo_fp_x
# input: f --> polinomio fabricado por fpx (objeto opaco) no nulo
# output: g = producto de los factores irreducibles mónicos de f, es decir,
# si f = c * f1^e1 * f2^e2 * ... * fr^er con los fi irreducibles mónicos
# distintos entre si, ei >= 1, c en fp, entonces g = f1 * f2 * ... * fr
def sqfree_fact_fpx(fpx, f):
    # Calcula el radical de f (producto de factores irreducibles sin repetición).
    # Algoritmo robusto para característica p.
    
    f = fpx.elem_de_tuple(fpx.conv_a_tuple(f))
    # Normalizar f para que sea mónico
    f = fpx.monico(f)
    
    if fpx.es_cero(f):
        raise ValueError("f no debe ser el polinomio cero")
    if fpx.grado(f) <= 0:
        return fpx.monico(f) # Devuelve 1 si es constante no nula

    p = fpx.fp.p

    # 1. Calcular derivada
    coeffs = fpx.conv_a_tuple(f)
    deriv = []
    for i in range(1, len(coeffs)):
        ci = coeffs[i]
        factor = i % p
        if factor == 0:
            deriv.append(0)
        else:
            deriv.append(fpx.fp.mult(ci, factor))
    fprime = fpx.elem_de_tuple(tuple(deriv))

    # 2. Caso derivada nula: f es una potencia p-ésima perfecta.
    # f(x) = h(x^p) -> rad(f) = rad(h)
    if fpx.es_cero(fprime):
        coeffs_f = fpx.conv_a_tuple(f)
        h_coeffs = []
        # Extraemos la raiz p-ésima (en Fp, a^(1/p) = a)
        for j in range(0, len(coeffs_f), p):
            h_coeffs.append(coeffs_f[j])
        h = fpx.elem_de_tuple(tuple(h_coeffs))
        return sqfree_fact_fpx(fpx, h)

    # 3. Caso general
    d = fpx.gcd(f, fprime)
    w = fpx.div(f, d)  # w es el producto de factores con multiplicidad NO div por p
    
    # El gcd 'd' todavía puede contener factores repetidos o potencias p-ésimas.
    # Debemos eliminar de 'd' cualquier factor que ya esté en 'w' para aislar
    # las partes puramente p-ésimas.
    y = fpx.gcd(d, w)
    z = d
    while not fpx.es_uno(y):
        z = fpx.div(z, y)
        y = fpx.gcd(z, w)
        
    # Ahora z contiene solo factores cuyas multiplicidades originales eran múltiplos de p.
    # Recursivamente calculamos el radical de z y lo multiplicamos por w.
    rad_z = sqfree_fact_fpx(fpx, z)
    
    return fpx.mult(w, rad_z)

# distinct-degree factorization
# input: fpx --> anillo_fp_x
# input: g --> polinomio de fpx (objeto opaco) que es producto de factores
# irreducibles mónicos distintos cada uno con multiplicidad uno
# output: [h1, h2, ..., hr], donde hi = producto de los factores irreducibles
# mónicos de h de grado = i, el último hr debe ser no nulo y por supuesto
# g = h1 * h2 * ... * hr
def didegr_fact_fpx(fpx, g):
    # Supuesto: g es cuadrado-libre y mónico.
    g = fpx.monico(fpx.elem_de_tuple(fpx.conv_a_tuple(g)))
    if fpx.es_cero(g) or fpx.es_uno(g):
        return []
        
    p = fpx.fp.p
    x = (0, 1) # Polinomio x en Fp[x] (tupla de coeficientes)
    h = g
    salida = []
    i = 1
    
    # 1. Iterar mientras 2*i <= grado(h)
    while 2 * i <= fpx.grado(h):
        # t = x^{p^i} mod h
        t = fpx.pot_mod(x, p**i, h)
        t_minus_x = fpx.suma(t, fpx.inv_adit(x))
        u = fpx.gcd(t_minus_x, h)
        u = fpx.monico(u)
        
        salida.append(u)
        
        if not fpx.es_uno(u):
            h = fpx.div(h, u)
            h = fpx.monico(h)
        i += 1
        
    # 2. Gestionar el residuo
    # Si h != 1, lo que queda tiene grado 'deg_h'. 
    # Pertenece a la posición deg_h-1 de la lista.
    if not fpx.es_uno(h):
        deg_h = fpx.grado(h)
        # Rellenar con 1s hasta llegar al índice correcto
        while len(salida) < deg_h - 1:
            salida.append(fpx.uno())
        salida.append(h)
        
    # 3. Limpiar unos redundantes al final (trailing ones)
    while len(salida) > 0 and fpx.es_uno(salida[-1]):
        salida.pop()
        
    return salida

# equal-degree factorization
# input: fpx --> anillo_fp_x (supondremos p impar)
# input: r --> int
# input: h --> polinomio de fpx (objeto opaco) que es producto de factores
# irreducibles mónicos distintos de grado r con multiplicidad uno
# output: [u1, ..., us], donde h = u1 * u2* ... * us y los ui son irreducibles
# mónicos de grado = r
def eqdegr_fact_fpx(fpx, r, h):
    # Cantor–Zassenhaus para factores de grado fijo r (p impar).
    # h = producto de irreducibles mónicos de grado r, sin repetición.
    import random
    p = fpx.fp.p
    if p % 2 == 0:
        raise ValueError("eqdegr_fact_fpx asume p impar")
    h = fpx.monico(fpx.elem_de_tuple(fpx.conv_a_tuple(h)))
    if fpx.es_uno(h):
        return []
    if fpx.grado(h) == r:
        return [h]
    n = fpx.grado(h)
    exp = (p**r - 1) // 2  # (q-1)/2 con q = p^r

    while True:
        # g(x) aleatorio (no nulo), deg(g) < deg(h)
        g_coeffs = [random.randrange(p) for _ in range(n)]
        if all(c == 0 for c in g_coeffs):
            continue
        g = fpx.elem_de_tuple(tuple(g_coeffs))
        # a = g^{(p^r-1)/2} mod h
        a = fpx.pot_mod(g, exp, h)
        a_minus_1 = fpx.suma(a, fpx.inv_adit(fpx.uno()))
        u = fpx.gcd(a_minus_1, h)
        u = fpx.monico(u)
        if fpx.es_uno(u) or fpx.es_igual(u, h):
            continue  # intento fallido, probamos otra g
        # Split encontrado
        v = fpx.div(h, u)
        u = fpx.monico(u); v = fpx.monico(v)
        return eqdegr_fact_fpx(fpx, r, u) + eqdegr_fact_fpx(fpx, r, v)

# multiplicidad de factor irreducible mónico
# input: fpx --> anillo_fp_x
# input: f --> polinomio de fpx (objeto opaco) no nulo
# input: u --> polinomio irreducible mónico (objeto opaco) de grado >= 1
# output: multiplicidad de u como factor de f, es decir, el entero e >= 0
# mas grande tal que u^e | f
def multiplicidad_fpx(fpx, f, u):
    f = fpx.elem_de_tuple(fpx.conv_a_tuple(f))
    u = fpx.monico(fpx.elem_de_tuple(fpx.conv_a_tuple(u)))
    if fpx.es_cero(f):
        return 0
    if fpx.es_uno(u):
        return 0
    if fpx.es_cero(u):
        raise ValueError("u no debe ser el polinomio cero")
    e = 0
    resto = f
    while True:
        q, r = fpx.divmod(resto, u)
        if fpx.es_cero(r):
            e += 1
            resto = q
        else:
            break
    return e

# factorización de Cantor-Zassenhaus
# input: fpx --> anillo_fp_x (supondremos p impar)
# input: f --> polinomio de fpx (objeto opaco)
# output: [(f1,e1), ..., (fr,er)] donde f = c * f1^e1 * ... * fr^er es la
# factorización completa de f en irreducibles mónicos fi con multiplicidad
# ei >= 1 y los fi son distintos entre si y por supuesto c es el coeficiente
# principal de f
def fact_fpx(fpx, f):                     # mantener esta implementación
    g = sqfree_fact_fpx(fpx, f)
    h = didegr_fact_fpx(fpx, g)
    irreducibles = []
    for r in range(len(h)):
        if fpx.grado(h[r]) > 0:
            irreducibles += eqdegr_fact_fpx(fpx, r+1, h[r])
    factorizacion = []
    for u in irreducibles:
        e = multiplicidad_fpx(fpx, f, u)
        factorizacion += [(u,e)]
    return factorizacion

# esta linea es para añadir la función de factorización de Cantor-Zassenhaus
# como un método de la clase anillo_fp_x
cf.anillo_fp_x.factorizar = fact_fpx

# square-free factorization
# input: fqx --> anillo_fq_x
# input: f --> polinomio fabricado por fqx (objeto opaco) no nulo
# output: g = producto de los factores irreducibles mónicos de f, es decir,
# si f = c * f1^e1 * f2^e2 * ... * fr^er con los fi irreducibles mónicos
# distintos entre si, ei >= 1, c en fq, entonces g = f1 * f2 * ... * fr
def sqfree_fact_fqx(fqx, f):
    # Calcula el radical de f en Fq[x]
    
    f = fqx.elem_de_tuple(fqx.conv_a_tuple(f))
    f = fqx.monico(f)
    
    if fqx.es_cero(f):
        raise ValueError("El polinomio no puede ser cero")
    if fqx.grado(f) <= 0:
        return fqx.monico(f)

    # Usamos la derivada que ya implementaste en cuerpos_finitos.py
    fprime = fqx.derivada(f)
    
    # 2. Caso derivada nula
    if fqx.es_cero(fprime):
        coeffs_f = fqx.conv_a_tuple(f)
        h_coeffs = []
        
        # Datos del cuerpo
        p = fqx.fq.fp.p
        q = fqx.fq.q
        # El exponente para la raíz p-ésima en Fq es q/p (es decir, p^(n-1))
        # ya que x^q = x => x^(q/p) = x^(1/p) en el sentido del automorfismo inverso
        exp_root = q // p 
        
        for j in range(0, len(coeffs_f), p):
            val = coeffs_f[j]
            # Raíz p-ésima del coeficiente
            root_val = fqx.fq.pot(val, exp_root)
            h_coeffs.append(root_val)
            
        h = fqx.elem_de_tuple(tuple(h_coeffs))
        return sqfree_fact_fqx(fqx, h)

    # 3. Caso general (idéntico a Fp pero llamando a métodos de fqx)
    d = fqx.gcd(f, fprime)
    w = fqx.div(f, d) # Parte libre de cuadrados "fácil"
    
    # Limpiar d de los factores que ya están en w
    y = fqx.gcd(d, w)
    z = d
    while not fqx.es_uno(y):
        z = fqx.div(z, y)
        y = fqx.gcd(z, w)
        
    # Recursión sobre la parte restante (potencias p-ésimas ocultas)
    rad_z = sqfree_fact_fqx(fqx, z)
    
    return fqx.mult(w, rad_z)

# distinct-degree factorization
# input: fqx --> anillo_fq_x
# input: g --> polinomio de fqx (objeto opaco) que es producto de factores
# irreducibles mónicos distintos cada uno con multiplicidad uno
# output: [h1, h2, ..., hr], donde hi = producto de los factores irreducibles
# mónicos de h de grado = i, el último hr debe ser no nulo y por supuesto
# g = h1 * h2 * ... * hr
def didegr_fact_fqx(fqx, g):
    g = fqx.monico(fqx.elem_de_tuple(fqx.conv_a_tuple(g)))
    if fqx.es_cero(g) or fqx.es_uno(g):
        return []
        
    q = fqx.fq.q
    # Construimos el polinomio x: 0 + 1*x
    x = fqx.elem_de_tuple((fqx.fq.cero(), fqx.fq.uno()))
    h = g
    salida = []
    i = 1
    
    # 1. Iterar mientras 2*i <= grado(h)
    while 2 * i <= fqx.grado(h):
        # t = x^{q^i} mod h (ojo: q, no p)
        t = fqx.pot_mod(x, q**i, h)
        t_minus_x = fqx.suma(t, fqx.inv_adit(x))
        u = fqx.gcd(t_minus_x, h)
        u = fqx.monico(u)
        
        salida.append(u)
        
        if not fqx.es_uno(u):
            h = fqx.div(h, u)
            h = fqx.monico(h)
        i += 1

    # 2. Gestionar el residuo
    if not fqx.es_uno(h):
        deg_h = fqx.grado(h)
        # Rellenar con 1s hasta llegar al índice correcto
        while len(salida) < deg_h - 1:
            salida.append(fqx.uno())
        salida.append(h)
        
    # 3. Limpiar unos redundantes al final
    while len(salida) > 0 and fqx.es_uno(salida[-1]):
        salida.pop()
        
    return salida



# equal-degree factorization
# input: fqx --> anillo_fq_x (supondremos q impar)
# input: r --> int
# input: h --> polinomio de fqx (objeto opaco) que es producto de factores
# irreducibles mónicos distintos de grado r con multiplicidad uno
# output: [u1, ..., us], donde h = u1 * u2* ... * us y los ui son irreducibles
# mónicos de grado = r
def eqdegr_fact_fqx(fqx, r, h):
    if fqx.es_cero(h) or fqx.es_uno(h):
        return []

    if fqx.grado(h) == r:
        return [h]

    q = fqx.fq.q
    to_factor = [h]
    irreducibles = []
    max_attempts = 100

    while to_factor:
        current = to_factor.pop()
        current_deg = fqx.grado(current)
        
        if current_deg == r:
            irreducibles.append(current)
            continue
            
        if current_deg < r:
            continue

        found = False
        attempts = 0
        while not found and attempts < max_attempts:
            # Generar polinomio aleatorio de grado menor que current_deg
            a_coeffs = [fqx.fq.aleatorio() for _ in range(current_deg)]
            a = fqx.elem_de_tuple(tuple(a_coeffs))
            
            if fqx.es_cero(a):
                attempts += 1
                continue

            exp = (q**r - 1) // 2
            b = fqx.pot_mod(a, exp, current)
            b_minus_1 = fqx.suma(b, fqx.inv_adit(fqx.uno()))
            g = fqx.gcd(b_minus_1, current)

            if not fqx.es_uno(g) and not fqx.es_igual(g, current):
                to_factor.append(g)
                to_factor.append(fqx.div(current, g))
                found = True

            attempts += 1

        if not found:
            irreducibles.append(current)

    return irreducibles


# multiplicidad de factor irreducible mónico
# input: fqx --> anillo_fq_x
# input: f --> polinomio de fqx (objeto opaco) no nulo
# input: u --> polinomio irreducible mónico (objeto opaco) de grado >= 1
# output: multiplicidad de u como factor de f, es decir, el entero e >= 0
# mas grande tal que u^e | f
def multiplicidad_fqx(fqx, f, u):
    e = 0
    temp = f
    
    # Asegurar que trabajamos con el polinomio original y el factor en forma mónica
    u_monic = fqx.monico(u)
    
    while True:
        # Calcular el resto de la división
        remainder = fqx.mod(temp, u_monic)
        
        # Verificar si el resto es cero (u divide a temp)
        if fqx.es_cero(remainder):
            e += 1
            # Dividir temp por u_monic
            temp = fqx.div(temp, u_monic)
        else:
            break
            
        # Prevenir bucle infinito por seguridad
        if e > fqx.grado(f):
            break
    
    return e


# factorización de Cantor-Zassenhaus
# input: fqx --> anillo_fq_x (supondremos q impar)
# input: f --> polinomio de fqx (objeto opaco)
# output: [(f1,e1), ..., (fr,er)] donde f = c * f1^e1 * ... * fr^er es la
# factorización completa de f en irreducibles mónicos fi con multiplicidad
# ei >= 1 y los fi son distintos entre si y por supuesto c es el coeficiente
# principal de f
def fact_fqx(fqx, f):                     # mantener esta implementación
    g = sqfree_fact_fqx(fqx, f)
    h = didegr_fact_fqx(fqx, g)
    irreducibles = []
    for r in range(len(h)):
        if fqx.grado(h[r]) > 0:
            irreducibles += eqdegr_fact_fqx(fqx, r+1, h[r])
    factorizacion = []
    for u in irreducibles:
        e = multiplicidad_fqx(fqx, f, u)
        factorizacion += [(u,e)]
    return factorizacion

# esta linea es para añadir la función de factorización de Cantor-Zassenhaus
# como un método de la clase anillo_fq_x
cf.anillo_fq_x.factorizar = fact_fqx
