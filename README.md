# Arquitectura computacional multiómica para el análisis del microbioma intestinal en Parkinson

Repositorio asociado al trabajo de grado **TG263006**, desarrollado en la Pontificia Universidad Javeriana.

## Objetivo general

Desarrollar una arquitectura computacional de integración multiómica para identificar rutas funcionales **microbio–gen–metabolito** asociadas a la enfermedad de Parkinson.

## Módulos

| Módulo | Descripción |
|---|---|
| [Aleatorización mendeliana](./MR) | Evaluación de asociaciones potencialmente causales entre taxones intestinales y enfermedad de Parkinson mediante datos GWAS. |
| [Metagenómica](./Metagenomics) | Procesamiento de datos shotgun, ensamblaje, binning, reconstrucción y caracterización de MAGs. |
| [Metatranscriptómica](./Metatranscriptomics) | Procesamiento de RNA, identificación de transcritos y análisis de expresión diferencial asociado a MAGs. |
| [16S rRNA](./16S) | Procesamiento de amplicones, generación de ASVs, asignación taxonómica y análisis de diversidad. |
| [Metabolómica](./Metabolomics) | Procesamiento y análisis estadístico para la identificación de metabolitos diferenciales. |
| [Integración multiómica](./Integration) | Integración taxonómica MAG–ASV y análisis funcional microbio–gen–metabolito. |

Cada módulo contiene su propio `README.md` con información sobre:

- estructura de archivos;
- datos de entrada;
- herramientas utilizadas;
- ejecución;
- parámetros principales;
- resultados generados.

## Objetivo general

Desarrollar una arquitectura computacional de integración multiómica para identificar rutas funcionales **microbio–gen–metabolito** asociadas a la enfermedad de Parkinson.
