# Arquitectura computacional multiómica para el análisis del microbioma intestinal en Parkinson

Repositorio asociado al trabajo de grado **TG263006**, desarrollado en la Pontificia Universidad Javeriana.

## Objetivo general

Desarrollar una arquitectura computacional de integración multiómica para identificar rutas funcionales **microbio–gen–metabolito** asociadas a la enfermedad de Parkinson.

## Módulos

| Módulo | Descripción |
|---|---|
| [Aleatorización mendeliana](./MR/README.md) | Evaluación de posibles relaciones causales entre taxones intestinales y Parkinson mediante datos GWAS. |
| [Metagenómica](./Metagenomics/README.md) | Procesamiento de lecturas shotgun para reconstruir MAGs de alta calidad y caracterizarlos taxonómica y funcionalmente. |
| [Metatranscriptómica](./Metatranscriptomics/README.md) | Procesamiento de RNA para identificar transcritos, asociarlos con MAGs y evaluar expresión diferencial. |
| [16S rRNA](./16S/README.md) | Procesamiento de amplicones para generar ASVs, asignar taxonomía y comparar la composición microbiana. |
| [Metabolómica](./Metabolomics/README.md) | Análisis de perfiles metabolómicos para identificar metabolitos diferenciales asociados con Parkinson. |
| [Integración multiómica](./Integration/README.md) | Integración de MAGs, ASVs, genes, rutas y metabolitos para identificar asociaciones funcionales. |

Cada módulo incluye su propio `README.md`, donde se documentan los elementos necesarios para comprender y reproducir el análisis, incluyendo:

- la organización y estructura del módulo;
- los datos y archivos de entrada requeridos;
- las herramientas y versiones utilizadas;
- el procedimiento de ejecución;
- las salidas y resultados generados.
