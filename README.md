# finite-fields-math-lib
# Librería de Álgebra Computacional: Cuerpos Finitos y Factorización

## Descripción
Este repositorio contiene una librería matemática avanzada implementada íntegramente en Python para el cálculo algebraico sobre **Cuerpos Finitos ($F_p$ y $F_q$)** y sus respectivos **Anillos de Polinomios ($F_p[x]$ y $F_q[x]$)**. 

El proyecto une el rigor del álgebra abstracta con la ingeniería de software orientada a la eficiencia computacional. Se ha diseñado estructurando el código de forma modular y aplicando algoritmos rápidos de vanguardia para reducir drásticamente la complejidad asintótica de las operaciones polinómicas y de factorización.

## Arquitectura del Proyecto

El ecosistema está dividido en tres módulos principales, construidos uno sobre otro:

### 1. Modelado Algebraico (`cuerpos_finitos.py`)
Implementación orientada a objetos de las estructuras algebraicas base:
* Clases para representar elementos de campos primos ($F_p$) y campos de extensión ($F_q$).
* Operaciones aritméticas elementales, inversos modulares (usando el Pequeño Teorema de Fermat y el Algoritmo de Euclides Extendido) y generación de tablas de operaciones.
* Clases para anillos de polinomios, incluyendo división euclídea, máximo común divisor (MCD) y tests de irreducibilidad de Rabin.

### 2. Algoritmos de Alta Eficiencia (`algoritmos_rapidos.py`)
Sobrecarga y optimización de las operaciones del anillo mediante técnicas algorítmicas avanzadas para manejar polinomios de alto grado:
* **Multiplicación Rápida:** Implementación del algoritmo de **Karatsuba** recursivo con umbrales de optimización empíricos.
* **Transformada Rápida de Fourier (FFT):** Implementación de la FFT y su inversa (IFFT) sobre cuerpos finitos para optimizar drásticamente la convolución y multiplicación de polinomios.
* **División y Matrices de Toeplitz:** Operaciones optimizadas de división e inversión polinómica modular utilizando el reverso de polinomios y multiplicación matriz-vector basada en matrices de Toeplitz.

### 3. Factorización Polinómica (`factorizacion.py`)
Implementación completa del algoritmo de **Cantor-Zassenhaus** para factorizar polinomios sobre $F_p$ y $F_q$ en sus componentes irreducibles. El pipeline incluye:
* **Square-Free Factorization** (Factorización libre de cuadrados) con manejo adaptado para la característica del cuerpo (derivadas nulas y cálculo de raíces p-ésimas).
* **Distinct-Degree Factorization** (Factorización por grados distintos).
* **Equal-Degree Factorization** (Factorización de igual grado).

## Complejidad y Optimización
El uso de algoritmos rápidos permite superar las limitaciones de la aritmética clásica $O(n^2)$. La integración de Karatsuba ($O(n^{1.58})$) y FFT ($O(n \log n)$) garantiza que las operaciones escalen de manera eficiente incluso en anillos de alta dimensión geométrica, una característica crucial para aplicaciones en criptografía y sistemas de corrección de errores.

## Tecnologías Utilizadas
* **Lenguaje:** Python 3 (Puro, sin dependencias externas pesadas).
* **Paradigmas:** Programación Orientada a Objetos (POO), Recursión y Divide y Vencerás.
* **Conceptos:** Álgebra Computacional, Teoría de Números, Optimización Algorítmica.

## Contacto
Estudiante de último año de Matemáticas y Ciencias de la Computación en la Universidad Complutense de Madrid. 
* **LinkedIn:** www.linkedin.com/in/juan-blanco-guillén 
