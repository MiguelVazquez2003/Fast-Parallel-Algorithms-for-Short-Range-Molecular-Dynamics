# Observaciones
Se tradujeron algunas descripciones, pero el codigo sigue la estructura original. Solamente se cambiaron algunos parametros para jugar con los algoritmos.

# Fast-Parallel-Algorithms-for-Short-Range-Molecular-Dynamics

S. J. Plimpton, J Comp Phys, 117, 1-19 (1995).  Three parallel algorithms for classical molecular dynamics are presented. The first assigns each processor a fixed subset of atoms; the second assigns each a fixed subset of inter--atomic forces to compute; the third assigns each a fixed spatial region.


## Pasos para ejecutar y visualizar la simulación:

### 1. Guardar el archivo `lj.in`
Asegúrate de guardar los cambios en el archivo `lj.in`.

### 2. Ejecutar LAMMPS:

1. Abre una terminal (PowerShell, cmd, o MSYS2).
2. Navega al directorio donde se encuentra tu archivo `lj.in`:
   ```sh
   cd ruta/al/directorio
3.Ejecuta LAMMPS con el archivo de entrada:
lmp_serial -in lj.in

### 3. Revisar los archivos de salida:
El archivo dump.lj se generará en el mismo directorio y contendrá las posiciones de los átomos.

### 4. Visualizar los resultados con OVITO:
Abre OVITO.
Ve a File > Load File y selecciona el archivo dump.lj.
Usa las herramientas de OVITO para visualizar la trayectoria de los átomos, analizar las propiedades de los materiales y crear animaciones.
