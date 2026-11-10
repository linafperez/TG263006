# 16S rRNA

Pipeline modular implementado en **Nextflow DSL2** para el procesamiento y análisis de datos de secuenciación de amplicones del gen 16S rRNA del microbioma intestinal en pacientes con enfermedad de Parkinson (PD) y controles sanos (HC).

El flujo integra control de calidad, remoción de primers, inferencia de variantes de secuencia de amplicón (ASVs), clasificación taxonómica, preprocesamiento de abundancias, caracterización de la composición microbiana y análisis de diversidad alfa y beta.

## Esquema general

[`Diagrama general 16S.pdf`](https://github.com/linafperez/TG263006/blob/main/16S/diagrama_16s.pdf)

El pipeline comprende las etapas de **validación de muestras y metadatos**, **agrupamiento de corridas técnicas**, **control de calidad**, **remoción de primers**, **denoising e inferencia de ASVs**, **entrenamiento del clasificador taxonómico**, **clasificación y filtrado taxonómico**, **preprocesamiento de abundancias**, **composición taxonómica**, **diversidad alfa** y **diversidad beta**.

La arquitectura se organiza jerárquicamente mediante:

```text
main.nf
│
├── workflows/
│   └── rrna_16s_analysis/
│
├── subworkflows/
│   ├── core/
│   └── local/
│
├── modules/
│   ├── core/
│   └── local/
│
├── bin/
└── conf/
```

Los módulos `core` contienen los procesos asociados principalmente con QIIME 2 y las herramientas de procesamiento de secuencias, mientras que los módulos `local` contienen procesos específicos desarrollados para validación, agrupamiento de muestras y análisis estadístico en R.

---

## Software y requerimientos computacionales

El pipeline fue validado con:

- **Nextflow v25.04.8**
- **Singularity v3.7.1**
- **QIIME 2 v2025.7.0**
- **Cutadapt v5.1**
- **DADA2 v1.30.0**
- **scikit-learn v1.4.2**
- **SILVA v138.2 SSURef NR99**
- **FastQC v0.12.1**
- **MultiQC v1.33**
- **R v4.4**

Los siguientes valores corresponden al consumo observado durante el procesamiento de las cohortes EUR y LATAM.

| **Software y versión**                                                                 | **Módulo / etapa**                                          | **CPU** | **RAM asignada (GB)** | **RAM pico (GB)** | **Almacenamiento** | **Tiempo** |
| -------------------------------------------------------------------------------------- | ----------------------------------------------------------- | ------: | --------------------: | ----------------: | -----------------: | ---------: |
| Nextflow v25.04.8 / scripts locales                                                    | Preparación de muestras (EUR)                               |       1 |                     4 |              0.05 |             4.0 GB |   00:19:35 |
| Nextflow v25.04.8 / scripts locales                                                    | Preparación de muestras (LATAM)                             |       1 |                     4 |              0.06 |             4.1 GB |   00:02:05 |
| FastQC v0.12.1 / MultiQC v1.33                                                         | Control de calidad (EUR)                                    |       4 |                    16 |              1.30 |             423 MB |   00:15:50 |
| FastQC v0.12.1 / MultiQC v1.33                                                         | Control de calidad (LATAM)                                  |       4 |                    16 |              1.20 |             187 MB |   00:06:00 |
| QIIME 2 v2025.7.0 / Cutadapt v5.1 / DADA2 v1.30.0 / scikit-learn v1.4.2 / SILVA v138.2 | Procesamiento, denoising y clasificación taxonómica (EUR)   |       4 |                    64 |             32.90 |             9.4 GB |   03:15:12 |
| QIIME 2 v2025.7.0 / Cutadapt v5.1 / DADA2 v1.30.0 / scikit-learn v1.4.2 / SILVA v138.2 | Procesamiento, denoising y clasificación taxonómica (LATAM) |       4 |                    64 |             30.80 |             9.3 GB |   03:12:42 |
| R v4.4                                                                                 | Preprocesamiento de abundancias (EUR)                       |       1 |                     4 |             <0.01 |             2.1 MB |   00:00:21 |
| R v4.4                                                                                 | Preprocesamiento de abundancias (LATAM)                     |       1 |                     4 |              0.01 |            0.55 MB |   00:00:42 |
| R v4.4                                                                                 | Composición taxonómica (EUR)                                |       1 |                     4 |              0.22 |            0.63 MB |   00:01:32 |
| R v4.4                                                                                 | Composición taxonómica (LATAM)                              |       1 |                     4 |              0.11 |            0.55 MB |   00:01:36 |
| R v4.4                                                                                 | Diversidad alfa (EUR)                                       |       1 |                     4 |              0.15 |            0.19 MB |   00:01:42 |
| R v4.4                                                                                 | Diversidad alfa (LATAM)                                     |       1 |                     4 |              0.19 |            0.13 MB |   00:01:25 |
| R v4.4                                                                                 | Diversidad beta y PERMANOVA (EUR)                           |       1 |                     4 |              0.31 |             4.2 MB |   00:02:13 |
| R v4.4                                                                                 | Diversidad beta y PERMANOVA (LATAM)                         |       1 |                     4 |              0.25 |             1.3 MB |   00:02:00 |

El tiempo total de pared del workflow fue de aproximadamente **02:15:31 para EUR** y **01:53:16 para LATAM**. Los resultados publicados ocuparon aproximadamente **14 GB por cohorte**.

Las etapas de mayor demanda computacional correspondieron al **denoising mediante DADA2** y al **entrenamiento del clasificador taxonómico**. Este último utilizó 4 CPU y solicitó 64 GB de RAM.

> **Nota:** los valores corresponden a las ejecuciones de referencia utilizadas durante el desarrollo y validación del pipeline. El consumo puede variar según el número de muestras, profundidad de secuenciación y recursos disponibles.

---

## Datasets de referencia

El pipeline fue aplicado a dos cohortes independientes con regiones hipervariables diferentes del gen 16S rRNA.

| **Característica**        | **EUR**                                                          | **LATAM**                                                          |
| ------------------------- | ---------------------------------------------------------------- | ------------------------------------------------------------------ |
| BioProject                | [PRJEB27564](https://www.ncbi.nlm.nih.gov/bioproject/PRJEB27564) | [PRJNA975118](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA975118) |
| Población                 | Europea                                                          | Latinoamericana                                                    |
| Muestras biológicas       | 124                                                              | 50                                                                 |
| Parkinson (PD)            | 63                                                               | 25                                                                 |
| Controles (HC)            | 61                                                               | 25                                                                 |
| Corridas de secuenciación | 124                                                              | 150                                                                |
| Región 16S                | V3–V4                                                            | V4–V5                                                              |
| Lecturas                  | paired-end                                                       | paired-end                                                         |
| Cohorte de referencia     | Pereira et al.                                                   | Forero-Rodríguez et al.                                            |

En la cohorte LATAM, las **150 corridas técnicas corresponden a 50 muestras biológicas**, con tres corridas por individuo. El pipeline agrupa las corridas correspondientes a una misma muestra antes del procesamiento posterior.

Los datos de secuenciación **no se distribuyen dentro de este repositorio** y deben obtenerse de sus repositorios públicos originales.

---

## Entradas

El pipeline requiere principalmente:

1. lecturas paired-end en formato FASTQ;
2. un archivo `samplesheet.csv`;
3. metadatos de las muestras;
4. secuencias de primers;
5. parámetros específicos de la región 16S analizada.

### 1. Lecturas FASTQ

Las lecturas deben encontrarse disponibles externamente al repositorio.

Ejemplo:

```text
/ruta/datos/
├── sample_01_R1.fastq.gz
├── sample_01_R2.fastq.gz
├── sample_02_R1.fastq.gz
└── sample_02_R2.fastq.gz
```

Para cohortes con varias corridas técnicas correspondientes a una misma muestra biológica, el pipeline permite agrupar y combinar las lecturas antes de las etapas de QIIME 2.

---

### 2. Samplesheet

El archivo de entrada define la asociación entre las muestras y sus archivos FASTQ.

Ejemplo conceptual:

```csv
sample,fastq_1,fastq_2
sample_01,/ruta/datos/sample_01_R1.fastq.gz,/ruta/datos/sample_01_R2.fastq.gz
sample_02,/ruta/datos/sample_02_R1.fastq.gz,/ruta/datos/sample_02_R2.fastq.gz
```

El pipeline incluye módulos específicos para validar la estructura del samplesheet antes de iniciar el procesamiento.

---

### 3. Metadatos

Los metadatos deben incluir, como mínimo, un identificador compatible con el utilizado en el samplesheet y el grupo experimental.

Ejemplo:

```text
sample_id    group
sample_01    PD
sample_02    Healthy
```

El módulo `metadata_check` verifica la correspondencia entre las muestras incluidas en los datos de secuenciación y los metadatos.

---

### 4. Primers

Los primers dependen de la región hipervariable amplificada y deben proporcionarse para cada cohorte.

#### Cohorte EUR — V3–V4

| Primer | Secuencia               |
| ------ | ----------------------- |
| 341F   | `CCTACGGGNGGCWGCAG`     |
| 785R   | `GACTACHVGGGTATCTAATCC` |

#### Cohorte LATAM — V4–V5

| Primer | Secuencia              |
| ------ | ---------------------- |
| 515F   | `GTGYCAGCMGCCGCGGTAA`  |
| 926R   | `CCGYCAATTYMTTTRAGTTT` |

Cutadapt elimina estas secuencias antes del denoising.

---

## Parámetros específicos por cohorte

Debido a que las cohortes fueron generadas utilizando regiones hipervariables diferentes, los parámetros de procesamiento fueron adaptados independientemente.

| **Parámetro**             |       **EUR** |     **LATAM** |
| ------------------------- | ------------: | ------------: |
| Región                    |         V3–V4 |         V4–V5 |
| Truncación forward        |        300 nt |        250 nt |
| Truncación reverse        |        228 nt |        230 nt |
| Rango de extracción SILVA |    300–500 nt |    240–420 nt |
| Confianza taxonómica      |           0.7 |           0.7 |
| `max_ee_f`                |           2.0 |           2.0 |
| `max_ee_r`                |           2.0 |           2.0 |
| `trunc_q`                 |             2 |             2 |
| `min_overlap`             |            12 |            12 |
| `pooling_method`          | `independent` | `independent` |
| `chimera_method`          |   `consensus` |   `consensus` |
| `n_reads_learn`           |     1,000,000 |     1,000,000 |

Los parámetros generales se encuentran definidos en:

```text
conf/params.yaml
```

y pueden modificarse según las características de la región amplificada y la calidad de las lecturas.

---

## Procesamiento taxonómico

Después del control de calidad, las lecturas son procesadas mediante QIIME 2.

El flujo incluye:

1. importación de las lecturas;
2. resumen de la calidad;
3. eliminación de primers con Cutadapt;
4. denoising mediante DADA2;
5. unión de lecturas paired-end;
6. remoción de quimeras;
7. generación de ASVs;
8. preparación de la referencia SILVA;
9. extracción de la región específica amplificada;
10. entrenamiento de un clasificador Naive Bayes;
11. asignación taxonómica;
12. filtrado de secuencias mitocondriales, cloroplastos y Eukaryota.

La referencia empleada corresponde a **SILVA v138.2 SSURef NR99**.

El clasificador se entrena específicamente para la región amplificada en cada cohorte, evitando aplicar directamente un clasificador entrenado sobre una región 16S diferente.

---

## Preprocesamiento de abundancias

Después de la clasificación taxonómica, las tablas son exportadas para su procesamiento en R.

Se realizan:

- eliminación de ASVs con abundancia total igual a cero;
- filtrado por prevalencia;
- cálculo de abundancias relativas;
- generación de perfiles de composición taxonómica.

En los análisis de referencia se excluyeron los ASVs presentes en menos del **10% de las muestras**.

Para la visualización taxonómica, los géneros que no superaron el **1% de abundancia relativa máxima en ninguna muestra** se agruparon como `Rare taxa`.

---

## Diversidad alfa

La diversidad intramuestra se evalúa mediante métricas complementarias:

- **Chao1**, como estimador de riqueza;
- **Shannon**, que integra riqueza y equidad;
- **Simpson**, con mayor ponderación de los taxones dominantes.

Las comparaciones entre PD y HC se consideran estadísticamente significativas cuando:

```text
P < 0.05
```

El uso de métricas complementarias evita caracterizar la diversidad alfa mediante un único índice.

---

## Diversidad beta

Debido a la naturaleza composicional de los datos de microbioma, la diversidad beta se calcula mediante:

1. reemplazo multiplicativo de ceros;
2. transformación **Centered Log-Ratio (CLR)**;
3. cálculo de distancia **Aitchison** mediante distancia euclidiana sobre los datos CLR;
4. representación mediante **PCoA**;
5. evaluación de diferencias entre grupos mediante **PERMANOVA**;
6. evaluación complementaria de dispersión mediante **PERMDISP**.

PERMANOVA se ejecuta con **999 permutaciones**, considerando:

```text
P < 0.05
```

como criterio de significancia.

PERMDISP se utiliza para evaluar la homogeneidad de las dispersiones entre grupos. Un resultado no significativo (`P > 0.05`) indica ausencia de evidencia de diferencias en dispersión y permite interpretar con mayor confianza una PERMANOVA significativa como una diferencia en la estructura composicional entre grupos.

---

## Estructura esperada de los datos

Los datos de secuenciación y metadatos se mantienen externos al repositorio.

```text
16S/
├── main.nf
├── nextflow.config
├── README.md
├── bin/
├── conf/
├── modules/
│   ├── core/
│   └── local/
├── subworkflows/
│   ├── core/
│   └── local/
└── workflows/
    └── rrna_16s_analysis/
        └── main.nf

/ruta/datos/
├── samplesheet.csv
├── metadata.tsv
├── primers.fasta
└── fastq/
    ├── sample_01_R1.fastq.gz
    ├── sample_01_R2.fastq.gz
    ├── sample_02_R1.fastq.gz
    └── sample_02_R2.fastq.gz
```

Los archivos FASTQ, bases de datos descargadas, resultados, archivos temporales y directorios de trabajo de Nextflow no forman parte del repositorio.

---

## Salidas principales

El pipeline genera:

- reportes de control de calidad mediante FastQC;
- reporte integrado mediante MultiQC;
- artefactos de importación y demultiplexado de QIIME 2;
- reportes posteriores a la remoción de primers;
- tablas de ASVs inferidas mediante DADA2;
- secuencias representativas de ASVs;
- estadísticas de denoising;
- referencia SILVA curada y región 16S extraída;
- clasificador taxonómico específico para la región analizada;
- asignaciones taxonómicas;
- tablas filtradas sin secuencias mitocondriales, cloroplastos ni Eukaryota;
- tablas de abundancia procesadas;
- perfiles de composición taxonómica;
- métricas y figuras de diversidad alfa;
- matrices CLR y distancias de Aitchison;
- coordenadas PCoA;
- resultados de PERMANOVA y PERMDISP;
- figuras de diversidad beta.

---

## Resultados obtenidos

La aplicación del pipeline a las dos cohortes permitió conservar, después del filtrado por prevalencia:

- **358 ASVs en EUR**;
- **275 ASVs en LATAM**.

La diversidad alfa evaluada mediante Chao1, Shannon y Simpson no presentó diferencias significativas entre pacientes con Parkinson y controles sanos en ninguna de las cohortes (`P > 0.05`).

La diversidad beta basada en distancia de Aitchison mostró diferencias significativas entre PD y HC:

| **Cohorte** | **R² PERMANOVA** | **F** | **P** |
| ----------- | ---------------: | ----: | ----: |
| EUR         |           0.0177 |  2.20 | 0.001 |
| LATAM       |           0.0312 |  1.55 | 0.015 |

Los tamaños de efecto fueron reducidos, pero la separación significativa se reprodujo en dos poblaciones con regiones hipervariables y contextos geográficos diferentes.

Estos resultados son utilizados posteriormente junto con los resultados metagenómicos, metatranscriptómicos y metabolómicos en la etapa de integración multiómica.

---

## Ejecución

El pipeline fue validado en un entorno Linux utilizando:

```text
Nextflow v25.04.8
Singularity v3.7.1
```

Desde el directorio `16S/`, la ejecución general se realiza mediante:

```bash
nextflow run main.nf \
  -profile singularity \
  --input /ruta/datos/samplesheet.csv \
  --metadata /ruta/datos/metadata.tsv \
  --primers /ruta/datos/primers.fasta
```

Los parámetros adicionales pueden definirse en:

```text
conf/params.yaml
```

o proporcionarse mediante la línea de comandos.

### Ejemplo EUR

```bash
nextflow run main.nf \
  -profile singularity \
  --input /ruta/EUR/samplesheet.csv \
  --metadata /ruta/EUR/metadata.tsv \
  --primers /ruta/EUR/primers.fasta \
  --trunc_len_f 300 \
  --trunc_len_r 228 \
  --max_ee_f 2.0 \
  --max_ee_r 2.0 \
  --trunc_q 2 \
  --min_overlap 12 \
  --pooling_method independent \
  --chimera_method consensus \
  --n_reads_learn 1000000
```

Para la cohorte EUR debe utilizarse adicionalmente el rango de extracción de SILVA correspondiente a **300–500 nt**.

### Ejemplo LATAM

```bash
nextflow run main.nf \
  -profile singularity \
  --input /ruta/LATAM/samplesheet.csv \
  --metadata /ruta/LATAM/metadata.tsv \
  --primers /ruta/LATAM/primers.fasta \
  --trunc_len_f 250 \
  --trunc_len_r 230 \
  --max_ee_f 2.0 \
  --max_ee_r 2.0 \
  --trunc_q 2 \
  --min_overlap 12 \
  --pooling_method independent \
  --chimera_method consensus \
  --n_reads_learn 1000000
```

Para LATAM debe emplearse el rango de extracción de SILVA correspondiente a **240–420 nt**.

Una ejecución interrumpida puede reanudarse utilizando:

```bash
nextflow run main.nf \
  -profile singularity \
  --input /ruta/datos/samplesheet.csv \
  --metadata /ruta/datos/metadata.tsv \
  --primers /ruta/datos/primers.fasta \
  -resume
```

Para garantizar la reproducibilidad deben conservarse el mismo directorio `work/` y la caché de Nextflow utilizados en la ejecución inicial.

---

## Configuración

Los principales archivos de configuración son:

```text
conf/
├── base.config
├── modules.config
├── params.yaml
├── test.config
└── test_gut_to_soil.yaml
```

- `base.config`: configuración general de procesos y recursos.
- `modules.config`: configuración específica de módulos.
- `params.yaml`: parámetros generales del análisis.
- `test.config`: configuración para pruebas.
- `test_gut_to_soil.yaml`: configuración utilizada durante la validación del pipeline con datos de prueba.

Cada módulo que requiere dependencias específicas contiene su correspondiente archivo `environment.yml`.

---

## Documentación detallada

El flujo completo, incluyendo decisiones internas, parámetros, filtros y salidas intermedias, puede incorporarse en:

[`Flujo completo del análisis 16S rRNA`](https://github.com/linafperez/TG263006/blob/main/16S/Anexo_16s.pdf)

---

## Relación con otros módulos

Los perfiles taxonómicos y ASVs generados por este pipeline son utilizados posteriormente en:

[`../Integration`](https://github.com/linafperez/TG263006/blob/main/Integration)

para su integración con los resultados obtenidos mediante:

- metagenómica;
- metatranscriptómica;
- metabolómica.

En particular, la etapa de integración permite relacionar la información taxonómica derivada de 16S con los MAGs reconstruidos mediante metagenómica y posteriormente con señales funcionales y metabolómicas.

El pipeline 16S constituye así una capa taxonómica independiente dentro de la arquitectura multiómica desarrollada para el estudio del microbioma intestinal asociado con enfermedad de Parkinson.

