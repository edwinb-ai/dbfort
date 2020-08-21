# Dinámica Browniana con Interacciones Hidrodinámicas

Este es un código modular para realizar dinámica browniana con interacciones hidrodinámicas
usando el clásico algoritmo de Ermak y McCammon.

Esta versión también incluye interacciones hidrodinámicas empleando el tensor de
Rotne-Prager-Yamakawa (RPY). El tensor se descompone mediante el algoritmo de Cholesky
para obtener los desplazamientos estocásticos.

## Uso general

El código está disponible para uso general y se compila mediante el comando

```
make -j
```

en la terminal. Se espera que se tenga `gfortran` instalado, así como `LAPACK`
y `BLAS`, ya que son necesarios para compilar y ejecutar la descomposición, respectivamente.

Una vez compilado, debe aparecer un ejecutable `dbrown` que podrá ser ejecutado y la simulación
comenzará de inmediato.

## Observables

Por el momento, el código solamente general las siguientes observables:

- Función de distribución radial (RDF).