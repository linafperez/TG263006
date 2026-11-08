# Metatranscriptómica

Pipeline metatranscriptómico modular implementado en **Nextflow** para caracterizar la actividad funcional del microbioma intestinal en pacientes con enfermedad de Parkinson (PD) y controles sanos (HC).

## Esquema general

[`Diagrama general metatranscriptomica.pdf`](./diagrama_metatranscriptomica.pdf)

El flujo comprende las etapas de **preprocesamiento**, **remoción de secuencias no deseadas**, **ensamblaje y predicción de transcritos**, **asignación a MAGs y cuantificación**, **expresión diferencial** y **control de calidad global**.

## Software y requerimientos computacionales

Los siguientes valores corresponden al consumo observado durante la ejecución del pipeline sobre **95 muestras**.

| Software y versión | Módulo / etapa | CPU | RAM asignada (GB) | RAM pico (GB) | Almacenamiento | Tiempo |
|---|---|---:|---:|---:|---:|---:|
| FastQC v0.12.1 | Control de calidad | 2 | 4 | 1.03 | 0.34 GB | 00:14:29 |
| fastp v0.23.4 | Control de calidad | 6 | 4 | 2.95 | 63 GB | 00:05:08 |
| SortMeRNA v4.3.7 | Remoción de rRNA y lecturas humanas | 12 | 16 | 5.80 | 55 GB | 04:34:16 |
| STAR v2.7.11b | Remoción de rRNA y lecturas humanas | 12 | 48 | 40.13 | 51 GB | 04:05:21 |
| rnaSPAdes v4.2.0 | Ensamblaje | 16 | 500 | 304.65 | 1.9 GB | 18:03:04 |
| rnaQUAST v2.3.0 | Ensamblaje | 6 | 20 | 1.72 | 0.00024 GB | 00:00:44 |
| TransDecoder v5.7.1 | Predicción de transcritos | 16 | 500 | 129.45 | 38.9 GB | 10:08:03 |
| DIAMOND v2.1.10 | Asignación a MAGs | 24 | 40 | 11.78 | 0.06 GB | 00:04:21 |
| Salmon v1.10.3 (índice) | Cuantificación | 8 | 32 | 7.63 | 7.1 GB | 00:12:36 |
| Salmon v1.10.3 (cuantificación) | Cuantificación | 8 | 24 | 20.78 | 19 GB | 01:48:24 |
| tximport v1.38.2 | Análisis de expresión | 6 | 128 | 96.20 | 0.36 GB | 00:14:05 |
| DESeq2 v1.50.2 | Análisis de expresión | 4 | 64 | 5.02 | 0.26 GB | 00:04:19 |
| MultiQC v1.25.2 | Reporte y control de calidad | 2 | 8 | 0.17 | 0.12 GB | 00:01:25 |

La ejecución completa del pipeline sobre las 95 muestras requirió aproximadamente **56 h 59 min**. Las etapas con mayor demanda de memoria fueron **rnaSPAdes** y **TransDecoder**, por lo que se recomienda el uso de infraestructura HPC.

## Dataset

| Característica | Descripción |
|---|---|
| BioProject | PRJNA782492 |
| Muestras | 95 |
| Parkinson (PD) | 46 |
| Controles (HC) | 49 |
| Tipo de datos | RNA-seq paired-end |
| Fuente original | Sequence Read Archive (SRA) |

## Entradas

El pipeline requiere cuatro grupos principales de archivos de entrada.

### 1. Lecturas metatranscriptómicas

Las lecturas paired-end de las 95 muestras deben encontrarse en:

```text
data/fastq/
```

Cada muestra contiene dos archivos comprimidos:

```text
<sample_id>_R1.fastq.gz
<sample_id>_R2.fastq.gz
```

Por ejemplo:

```text
PD_01_R1.fastq.gz
PD_01_R2.fastq.gz
HC_01_R1.fastq.gz
HC_01_R2.fastq.gz
```

### 2. Samplesheet

La metadata y las rutas de los archivos de entrada se definen en:

```text
data/samplesheet_95samples.csv
```

El archivo utiliza la siguiente estructura:

```text
sample_id,sra_id,condition,fastq_1,fastq_2
PD_01,SRR20608121,PD,data/fastq/PD_01_R1.fastq.gz,data/fastq/PD_01_R2.fastq.gz
PD_02,SRR20608120,PD,data/fastq/PD_02_R1.fastq.gz,data/fastq/PD_02_R2.fastq.gz
...
HC_49,SRR20608188,HC,data/fastq/HC_49_R1.fastq.gz,data/fastq/HC_49_R2.fastq.gz
```

| Columna | Descripción |
|---|---|
| `sample_id` | Identificador interno de la muestra |
| `sra_id` | Accesión original en SRA |
| `condition` | Grupo experimental: `PD` o `HC` |
| `fastq_1` | Ruta a la lectura forward |
| `fastq_2` | Ruta a la lectura reverse |

### 3. Referencias

Las bases de referencia utilizadas por el pipeline se almacenan en:

```text
data/references/
```

con la siguiente estructura:

```text
data/references/
├── human_GRCh38_p14_GENCODE_v44/
├── mags/
└── rrna_database_release_v4/
```

#### Referencia humana

```text
data/references/human_GRCh38_p14_GENCODE_v44/
```

Contiene:

```text
GRCh38.p14.genome.fa
gencode.v44.annotation.gtf
STAR_IDX/
```

El índice `STAR_IDX/` es utilizado por **STAR** para identificar y remover lecturas de origen humano.

#### Base de datos de rRNA

```text
data/references/rrna_database_release_v4/
```

Contiene las bases de datos e índices utilizados por **SortMeRNA**.

Entre los archivos disponibles se encuentran:

```text
smr_v4.3_default_db.fasta
smr_v4.3_fast_db.fasta
smr_v4.3_sensitive_db.fasta
smr_v4.3_sensitive_db_rfam_seeds.fasta
idx/
```

#### MAGs

Los genomas ensamblados y desreplicados provenientes del análisis metagenómico se encuentran en:

```text
data/references/mags/dereplicated_genomes/
```

Este directorio contiene los archivos FASTA correspondientes a los MAGs utilizados posteriormente para relacionar la actividad transcripcional con genomas microbianos específicos.

### 4. Anotaciones de MAGs

Las anotaciones derivadas de los MAGs utilizadas durante las etapas de asociación y análisis funcional se almacenan en:

```text
data/mag_annotations/
```

## Estructura esperada de los datos

```text
data/
├── fastq/
│   ├── PD_01_R1.fastq.gz
│   ├── PD_01_R2.fastq.gz
│   ├── ...
│   ├── HC_49_R1.fastq.gz
│   └── HC_49_R2.fastq.gz
│
├── samplesheet_95samples.csv
│
├── mag_annotations/
│
└── references/
    ├── human_GRCh38_p14_GENCODE_v44/
    │   ├── GRCh38.p14.genome.fa
    │   ├── gencode.v44.annotation.gtf
    │   └── STAR_IDX/
    │
    ├── rrna_database_release_v4/
    │   ├── smr_v4.3_default_db.fasta
    │   ├── smr_v4.3_fast_db.fasta
    │   ├── smr_v4.3_sensitive_db.fasta
    │   ├── smr_v4.3_sensitive_db_rfam_seeds.fasta
    │   └── idx/
    │
    └── mags/
        └── dereplicated_genomes/
            ├── M_COMEBin_*.fa
            ├── M_MetaBAT2_*.fa
            ├── M_SemiBin2_*.fa
            ├── M_Vamb_*.fa
            └── ...
```

## Salidas principales

El pipeline genera:

- reportes FastQC de lecturas raw y procesadas;
- reportes de filtrado de fastp;
- logs de SortMeRNA;
- métricas de alineamiento de STAR;
- transcritos ensamblados mediante rnaSPAdes;
- reportes de evaluación del ensamblaje mediante rnaQUAST;
- ORFs y CDS predichos mediante TransDecoder;
- asociaciones entre ORFs y MAGs generadas con DIAMOND;
- cuantificaciones por muestra mediante Salmon;
- matrices de expresión a nivel de genes;
- resultados de expresión diferencial entre PD y HC;
- reporte integrado de calidad y resultados mediante MultiQC.

## Resultados obtenidos

La aplicación del pipeline a las 95 muestras permitió identificar:

- **517 genes codificados por MAGs con expresión diferencial**;
- **480 genes con mayor expresión en PD**;
- **37 genes con mayor expresión en HC**.

Estos resultados son utilizados posteriormente en la etapa de integración multiómica para relacionar actividad génica, MAGs, funciones biológicas y metabolitos.

## Ejecución

La ejecución del pipeline se gestiona mediante:

[`metatranscriptomics_pipeline.sh`](https://github.com/linafperez/TG263006/blob/main/Metatranscriptomics/metatranscriptomics_pipeline.sh)

Este script actúa como punto de entrada para enviar el workflow al planificador **SLURM** mediante `sbatch` y ejecutar Nextflow utilizando los perfiles `slurm` y `singularity`.

Desde el directorio `Metatranscriptomics/`, se encuentran disponibles tres modos de ejecución:

```bash
./metatranscriptomics_pipeline.sh stub
```

Ejecuta una validación aislada del pipeline utilizando el entorno de prueba.

```bash
./metatranscriptomics_pipeline.sh real
```

Envía una ejecución completa sobre las 95 muestras.

```bash
./metatranscriptomics_pipeline.sh resume
```

Reanuda una ejecución productiva interrumpida utilizando la opción `-resume` de Nextflow.

El script envía los trabajos al clúster mediante `sbatch` y configura automáticamente los parámetros principales de SLURM, el entorno de Singularity y la ejecución de `main.nf`.

Antes de ejecutar el pipeline deben estar disponibles las lecturas, referencias y anotaciones descritas en la sección **Entradas**, respetando la estructura de directorios esperada.

## Documentación detallada

El flujo completo, incluyendo decisiones internas, filtros, criterios de selección y salidas intermedias, se encuentra en:

[`Flujo completo de metatranscriptomica.pdf`](./Anexo_Metatranscriptomica.pdf)

## Relación con otros módulos

Los MAGs utilizados por este pipeline provienen del módulo de:

[`../Metagenomics`](../Metagenomics)

Los resultados de expresión diferencial son posteriormente utilizados en:

[`../Integration`](../Integration)

para la construcción de asociaciones funcionales **microbio–gen–metabolito**.
