import cuerpos_finitos as cf

# ====================================================================
# FUNCIONES AUXILIARES (OPTIMIZADAS)
# ====================================================================

def _poly_reverso(anillo, A, n):
    """Calcula el reverso Rev_n(A) = x^(n-1) * A(1/x)."""
    if anillo.es_cero(A):
        return (anillo.campo.cero(),) * n
    
    A_list = list(A)
    # Ajustar longitud exactamente a n
    if len(A_list) < n:
        A_list.extend([anillo.campo.cero()] * (n - len(A_list)))
    else:
        A_list = A_list[:n]
        
    return anillo.elem_de_tuple(tuple(reversed(A_list)))

def _poly_inv_mod_xk(anillo, B, k):
    """Calcula la inversa de B modulo x^k usando Newton.
    B * I = 1 mod x^k."""
    campo = anillo.campo
    
    # Caso base
    if k == 1:
        return (campo.inv_mult(B[0]),)

    # Recursión: k0 = ceil(k/2)
    k0 = (k + 1) // 2
    I0 = _poly_inv_mod_xk(anillo, B, k0)

    # Iteración de Newton: I = I0 + I0 * (1 - B*I0) (mod x^k)
    # Solo necesitamos los primeros k coeficientes de B
    B_trunc = B[:k]
    
    # B * I0 (usando Karatsuba)
    B_I0 = anillo.mult_fast(B_trunc, I0)[:k]
    
    # 1 - B*I0
    # Calculamos término a término
    diff = [campo.inv_adit(c) for c in B_I0]
    # Sumar 1 al término constante
    diff[0] = campo.suma(campo.uno(), diff[0])
    
    # Rellenar con ceros hasta k si es necesario
    while len(diff) < k:
        diff.append(campo.cero())
    
    # Termino de corrección: I0 * (1 - B*I0)
    term_mult = anillo.mult_fast(I0, tuple(diff))[:k]
    
    # Resultado final
    res = anillo.suma(I0, term_mult)
    return res[:k]


# ====================================================================
# VERSIONES SOBRE Fp
# ====================================================================

def fp_x_mult_karatsuba(fpx, f, g):
    f = fpx.reducir(f)
    g = fpx.reducir(g)
    
    # Si alguno es cero o vacío
    if not f or not g:
        return fpx.cero()
        
    n = max(len(f), len(g))
    
    # Umbral empírico, 16 suele ser bueno
    if n < 16: 
        return fpx.mult(f, g) 

    k = n // 2
    
    f0 = f[:k]
    f1 = f[k:]
    g0 = g[:k]
    g1 = g[k:]
    
    # Recursión
    p1 = fp_x_mult_karatsuba(fpx, f0, g0)
    p2 = fp_x_mult_karatsuba(fpx, f1, g1)
    
    sum_f = fpx.suma(f0, f1)
    sum_g = fpx.suma(g0, g1)
    p3 = fp_x_mult_karatsuba(fpx, sum_f, sum_g)
    
    # p_middle = p3 - p1 - p2
    p_middle = fpx.suma(p3, fpx.inv_adit(fpx.suma(p1, p2)))
    
    # Combinar: p1 + x^k * p_middle + x^(2k) * p2
    # Usamos tuplas de ceros para el shift, es más eficiente que rellenar listas
    fp = fpx.campo
    zeros_k = (fp.cero(),) * k
    zeros_2k = (fp.cero(),) * (2*k)
    
    term_mid = zeros_k + p_middle
    term_high = zeros_2k + p2
    
    res = fpx.suma(p1, term_mid)
    res = fpx.suma(res, term_high)
    
    return fpx.reducir(res)

# Inyección en la clase
cf.anillo_fp_x.mult_fast = fp_x_mult_karatsuba


def fp_toep_inf_vec(fp, n, a, b):
    fpx = cf.anillo_fp_x(fp)
    # La multiplicación de Toeplitz inferior es equivalente a multiplicar polinomios
    # y quedarse con los primeros n términos.
    C = fpx.mult_fast(a, b)
    # Aseguramos longitud n con padding si hace falta
    C_padded = C + (fp.cero(),) * n
    return C_padded[:n]


def fp_toep_sup_vec(fp, n, a, b):
    fpx = cf.anillo_fp_x(fp)
    # Para Toeplitz superior, invertimos 'a' y tomamos una ventana superior
    A_rev = tuple(reversed(a))
    C = fpx.mult_fast(A_rev, b)
    # Padding extra para seguridad en el slice
    C_padded = C + (fp.cero(),) * (2*n)
    return C_padded[n-1 : 2*n-1]


def fp_toep_vec(fp, n, a, b):
    fpx = cf.anillo_fp_x(fp)
    # Construcción del vector que representa la matriz Toeplitz completa
    # a = [columna_inv... | diag | fila...]
    # Estructura: a[0]..a[n-2] es parte inferior invertida, a[n-1] diagonal, a[n].. fin es superior
    
    # Reconstruimos el polinomio A
    # Nota: a[:n-1] son los elementos t_{n-1,0}... t_{1,0}
    # a[n-1] es t_{0,0}
    # a[n:] son t_{0,1}... t_{0,n-1}
    
    A_list = list(reversed(a[:n-1])) + [a[n-1]] + list(a[n:])
    A_poly = tuple(A_list)
    
    C = fpx.mult_fast(A_poly, b)
    C_padded = C + (fp.cero(),) * (2*n)
    
    # El resultado válido está en la ventana central
    return C_padded[n-1 : 2*n-1]


def fp_toep_inf_inv(fp, n, a):
    # La inversa de T_inf(a) es T_inf(a^-1 mod x^n)
    fpx = cf.anillo_fp_x(fp)
    res = _poly_inv_mod_xk(fpx, a, n)
    # Rellenar con ceros si quedó corto
    return (res + (fp.cero(),) * n)[:n]


def fp_toep_sup_inv(fp, n, a):
    # Por isomorfismo, es igual al caso inferior
    fpx = cf.anillo_fp_x(fp)
    res = _poly_inv_mod_xk(fpx, a, n)
    return (res + (fp.cero(),) * n)[:n]


def fp_x_divmod(fpx, f, g):
    f = fpx.reducir(f)
    g = fpx.reducir(g)
    
    n = fpx.grado(f)
    m = fpx.grado(g)
    
    if n < m:
        return fpx.cero(), f
        
    k = n - m + 1
    
    # 1. Reversos
    f_rev = _poly_reverso(fpx, f, n + 1)
    g_rev = _poly_reverso(fpx, g, m + 1)
    
    # 2. Inversa de g_rev mod x^k
    g_rev_inv = _poly_inv_mod_xk(fpx, g_rev, k)
    
    # 3. Cociente reverso. Usamos toep_inf_vec para eficiencia y recorte automático
    # q_rev = f_rev * g_rev_inv mod x^k
    q_rev = fp_toep_inf_vec(fpx.fp, k, g_rev_inv, f_rev[:k])
    
    # 4. Recuperar q
    q = _poly_reverso(fpx, q_rev, k)
    
    # 5. Resto: r = f - g*q
    prod_qg = fpx.mult_fast(q, g)
    r = fpx.suma(f, fpx.inv_adit(prod_qg))
    
    return fpx.reducir(q), fpx.reducir(r)

cf.anillo_fp_x.divmod_fast = fp_x_divmod


def fp_fft(fp, g, k, a):
    target_n = 2**k
    a_list = list(a)
    if len(a_list) < target_n:
        a_list.extend([fp.cero()] * (target_n - len(a_list)))
    a = tuple(a_list)

    n = len(a) # Ahora n siempre es 2^k
    if n == 1:
        return a
    
    g_sq = fp.mult(g, g)
    
    a_even = a[::2]
    a_odd = a[1::2]
    
    U = fp_fft(fp, g_sq, k - 1, a_even)
    V = fp_fft(fp, g_sq, k - 1, a_odd)
    
    res = [None] * n
    w_j = fp.uno()
    
    mid = n // 2
    for j in range(mid):
        T = fp.mult(w_j, V[j])
        res[j] = fp.suma(U[j], T)
        res[j + mid] = fp.suma(U[j], fp.inv_adit(T))
        w_j = fp.mult(w_j, g)
        
    return tuple(res)


def fp_ifft(fp, g, k, b):
    # IDFT(a) = (1/n) * DFT(a, g^-1)
    g_inv = fp.inv_mult(g)
    a_rev = fp_fft(fp, g_inv, k, b)
    
    n_val = fp.elem_de_int(2**k)
    n_inv = fp.inv_mult(n_val)
    
    return tuple(fp.mult(x, n_inv) for x in a_rev)


# ====================================================================
# VERSIONES SOBRE Fq
# ====================================================================

def fq_x_mult_karatsuba(fqx, f, g):
    # Misma lógica que Fp
    f = fqx.reducir(f)
    g = fqx.reducir(g)
    if not f or not g:
        return fqx.cero()
        
    n = max(len(f), len(g))
    if n < 16:
        return fqx.mult(f, g)

    k = n // 2
    f0, f1 = f[:k], f[k:]
    g0, g1 = g[:k], g[k:]

    p1 = fq_x_mult_karatsuba(fqx, f0, g0)
    p2 = fq_x_mult_karatsuba(fqx, f1, g1)
    p3 = fq_x_mult_karatsuba(fqx, fqx.suma(f0, f1), fqx.suma(g0, g1))

    p_middle = fqx.suma(p3, fqx.inv_adit(fqx.suma(p1, p2)))

    fq = fqx.campo
    zeros_k = (fq.cero(),) * k
    zeros_2k = (fq.cero(),) * (2*k)
    
    res = fqx.suma(p1, zeros_k + p_middle)
    res = fqx.suma(res, zeros_2k + p2)
    return fqx.reducir(res)

cf.anillo_fq_x.mult_fast = fq_x_mult_karatsuba


def fq_toep_inf_vec(fq, n, a, b):
    fqx = cf.anillo_fq_x(fq)
    C = fqx.mult_fast(a, b)
    return (C + (fq.cero(),)*n)[:n]


def fq_toep_sup_vec(fq, n, a, b):
    fqx = cf.anillo_fq_x(fq)
    A_rev = tuple(reversed(a))
    C = fqx.mult_fast(A_rev, b)
    return (C + (fq.cero(),)*(2*n))[n-1 : 2*n-1]


def fq_toep_vec(fq, n, a, b):
    fqx = cf.anillo_fq_x(fq)
    A_list = list(reversed(a[:n-1])) + [a[n-1]] + list(a[n:])
    C = fqx.mult_fast(tuple(A_list), b)
    return (C + (fq.cero(),)*(2*n))[n-1 : 2*n-1]


def fq_toep_inf_inv(fq, n, a):
    fqx = cf.anillo_fq_x(fq)
    res = _poly_inv_mod_xk(fqx, a, n)
    return (res + (fq.cero(),)*n)[:n]


def fq_toep_sup_inv(fq, n, a):
    fqx = cf.anillo_fq_x(fq)
    res = _poly_inv_mod_xk(fqx, a, n)
    return (res + (fq.cero(),)*n)[:n]


def fq_x_divmod(fqx, f, g):
    f = fqx.reducir(f)
    g = fqx.reducir(g)
    
    n = fqx.grado(f)
    m = fqx.grado(g)
    
    if n < m:
        return fqx.cero(), f
    
    k = n - m + 1
    
    f_rev = _poly_reverso(fqx, f, n + 1)
    g_rev = _poly_reverso(fqx, g, m + 1)
    
    g_rev_inv = _poly_inv_mod_xk(fqx, g_rev, k)
    q_rev = fq_toep_inf_vec(fqx.fq, k, g_rev_inv, f_rev[:k])
    
    q = _poly_reverso(fqx, q_rev, k)
    r = fqx.suma(f, fqx.inv_adit(fqx.mult_fast(q, g)))
    
    return fqx.reducir(q), fqx.reducir(r)

cf.anillo_fq_x.divmod_fast = fq_x_divmod


def fq_fft(fq, g, k, a):
    target_n = 2**k
    a_list = list(a)
    if len(a_list) < target_n:
        a_list.extend([fq.cero()] * (target_n - len(a_list)))
    a = tuple(a_list)

    n = len(a)
    if n == 1:
        return a
    
    g_sq = fq.mult(g, g)
    U = fq_fft(fq, g_sq, k - 1, a[::2])
    V = fq_fft(fq, g_sq, k - 1, a[1::2])
    
    res = [None] * n
    w_j = fq.uno()
    mid = n // 2
    
    for j in range(mid):
        T = fq.mult(w_j, V[j])
        res[j] = fq.suma(U[j], T)
        res[j + mid] = fq.suma(U[j], fq.inv_adit(T))
        w_j = fq.mult(w_j, g)
        
    return tuple(res)


def fq_ifft(fq, g, k, b):
    g_inv = fq.inv_mult(g)
    a_rev = fq_fft(fq, g_inv, k, b)
    n_inv = fq.inv_mult(fq.elem_de_int(2**k))
    return tuple(fq.mult(x, n_inv) for x in a_rev)
