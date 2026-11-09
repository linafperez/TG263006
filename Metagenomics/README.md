# Metagenómica

Pipeline metagenómico modular implementado en **Nextflow** para caracterizar la composición taxonómica y el potencial funcional del microbioma intestinal en pacientes con enfermedad de Parkinson (PD) y controles sanos (HC).

## Esquema general

[`Diagrama general metagenomica.pdf`](./diagrama_metagenomica.pdf)

El flujo comprende las etapas de **obtención de lecturas**, **preprocesamiento y control de calidad**, **remoción de secuencias humanas**, **ensamblaje metagenómico**, **binning y refinamiento de MAGs**, **clasificación taxonómica y análisis filogenómico**, **predicción de genes y anotación funcional**, **estimación de abundancia** y **evaluación global del procesamiento**.

## Software y requerimientos computacionales

Los siguientes valores reproducen los requerimientos y el consumo computacional estimados para **95 muestras**, presentados en la Tabla 2 del trabajo de grado.

| Software y versión | Módulo / etapa | CPU | RAM asignada (GB) | RAM pico (GB) | Almacenamiento | Tiempo |
|---|---|---:|---:|---:|---:|---:|
| SRA Toolkit v3.4.1 / pigz v2.8 | Obtención de lecturas | 12 | 19 | 1.60 | 7.40 GB | 00:15:08 |
| FastQC v0.12.1 | Control de calidad | 4 | 19 | 1.59 | 1.75 MB | 00:04:43 |
| fastp v1.0.1 | Filtrado de lecturas | 16 | 19 | 1.41 | 2.93 GB | 00:02:31 |
| Bowtie2 v2.5.4 | Eliminación del hospedador | 16 | 19 | 3.84 | 7.12 GB | 00:21:47 |
| MEGAHIT v1.2.9 | Ensamblaje | 16 | 189 | 150.80 | 2.37 GB | 1-10:25:34 |
| metaSPAdes v4.2.0 | Ensamblaje | 16 | 811 | 648.74 | 8.42 GB | 3-04:12:15 |
| MetaQUAST v5.3.0 | Evaluación del ensamblaje | 6 | 20 | 4.00 | 1.90 MB | 00:06:11 |
| CoverM v0.7.0 | Cobertura de contigs | 16 | 100 | 32.17 | 117.69 GB | 03:08:00 |
| COMEBin v1.0.4 | Binning | 16 | 194 | 155.14 | 2.18 GB | 5-05:05:47 |
| SemiBin2 v1.5.0 | Binning | 16 | 243 | 193.63 | 1.92 GB | 11:48:05 |
| MetaBAT2 v2.18 | Binning | 16 | 100 | 5.35 | 1.65 GB | 01:06:47 |
| Vamb v5.0.4 | Binning | 16 | 100 | 4.85 | 2.22 GB | 01:57:58 |
| DAS Tool v1.1.7 | Refinamiento de bins | 16 | 128 | 3.01 | 1.16 GB | 00:29:54 |
| CheckM2 v1.1.0 | Evaluación de MAGs | 16 | 25 | 18.12 | 3.08 GB | 00:25:47 |
| GUNC v1.0.6 | Coherencia taxonómica | 16 | 25 | 11.70 | 13.00 GB | 01:14:04 |
| dRep v3.6.2 | Desreplicación de MAGs | 16 | 100 | 7.80 | 784.63 MB | 00:55:58 |
| GTDB-Tk v2.6.1 | Clasificación taxonómica | 16 | 126 | 100.46 | 151.03 GB | 02:55:22 |
| PhyloPhlAn v3.1.1 / IQ-TREE v3.0.1 | Análisis filogenómico | 16 | 120 | 30.29 | 3.76 GB | 1-04:01:26 |
| GeneMarkS-2 v1.15 | Predicción de genes | 1 | 20 | 0.39 | 9.25 MB | 00:00:57 |
| eggNOG-mapper v2.1.13 | Anotación funcional | 5 | 30 | 6.87 | 48.10 GB | 00:41:35 |
| InterProScan v5.59-91.0 | Anotación funcional | 16 | 128 | 8.99 | 15.02 GB | 01:59:23 |
| CoverM v0.7.0 | Abundancia de MAGs | 16 | 100 | 25.83 | 1.45 MB | 02:21:05 |
| MultiQC v1.35 | Reporte global | 2 | 8 | 1.50 | 30 MB | 00:03:00 |

Los tiempos por etapa sumaron aproximadamente **293.72 horas**. Las etapas con mayor demanda de memoria fueron **metaSPAdes**, **SemiBin2** y **COMEBin**, por lo que se recomienda el uso de infraestructura HPC. Según la nota de la tabla, la RAM pico corresponde a MaxRSS de SLURM y el tiempo a la mayor duración por tarea; el almacenamiento excluye entradas enlazadas y archivos internos de Nextflow. Los tiempos utilizan el formato `días-hh:mm:ss` o `hh:mm:ss`.

## Dataset

| Característica | Descripción |
|---|---|
| BioProject | PRJNA782492 |
| Muestras | 95 |
| Parkinson (PD) | 46 |
| Controles (HC) | 49 |
| Tipo de datos | Secuenciación metagenómica shotgun paired-end |
| Fuente original | Sequence Read Archive (SRA) |

El manifiesto utiliza `CT` para identificar los controles, equivalente al grupo HC.

## Entradas

El pipeline requiere cuatro grupos principales de archivos de entrada.

### 1. Lecturas metagenómicas

El pipeline admite archivos FASTQ locales o su adquisición desde SRA mediante BioProject y manifiesto TSV.

Para el modo local, las lecturas paired-end pueden organizarse en:

```text
/ruta/datos/fastq/
```

Cada muestra contiene dos archivos:

```text
<sample_id>_R1.fastq.gz
<sample_id>_R2.fastq.gz
```

Las rutas se especifican mediante `--input`. En modo SRA se utilizan `--sra-project` y `--sra-samples`; las corridas del mismo BioSample se combinan antes del procesamiento.

### 2. Samplesheet y manifiesto de muestras

El manifiesto de las 95 muestras se encuentra en:

```text
assets/raw/PRJNA782492_metagenomic_samples.tsv
```

El archivo utiliza la siguiente estructura, con columnas separadas por tabuladores:

```text
biosample_accession	group
SAMN29531760	PD
SAMN29531788	CT
```

Para archivos FASTQ locales se proporciona la plantilla:

```text
assets/samplesheet.csv
```

Contiene un ejemplo que debe adaptarse a las lecturas disponibles:

```text
sample,fastq_1,fastq_2
sample_A,/path/to/sample_A_R1.fastq.gz,/path/to/sample_A_R2.fastq.gz
```

| Columna | Descripción |
|---|---|
| `biosample_accession` | Accesión BioSample del manifiesto SRA |
| `group` | Grupo experimental: `PD` o `CT` |
| `sample` | Identificador de la muestra en el samplesheet local |
| `fastq_1` | Ruta a la lectura forward |
| `fastq_2` | Ruta a la lectura reverse |

`assets/samplesheet.csv` es una plantilla local, no el manifiesto completo de 95 muestras. Las rutas relativas se resuelven respecto al CSV y puede añadirse la columna `group`.

### 3. Referencias y bases de datos

Las referencias se proporcionan en el entorno de ejecución. Las siguientes rutas ilustrativas corresponden a la organización de `bin/prepare_databases.sh`:

```text
/ruta/referencias/
```

#### Referencia humana

```text
/ruta/referencias/human/GRCh38.p14/GRCh38.p14.fa
/ruta/referencias/human/bowtie2/GRCh38_p14
```

Contiene el genoma **GRCh38.p14 de GENCODE v44** y su índice Bowtie2. El parámetro `host_bowtie2_index` recibe el prefijo del índice, acompañado de sus seis archivos `.bt2` o `.bt2l`.

#### Calidad, taxonomía y filogenómica

```text
/ruta/referencias/checkm2/1.1.0/
/ruta/referencias/gunc/1.0.6/progenomes_2.1/
/ruta/referencias/gtdbtk/release226/
/ruta/referencias/phylophlan/phylophlan/
```

Se configuran mediante `checkm2_db`, `gunc_db`, `gtdbtk_db` y `phylophlan_db`. CheckM2 requiere `uniref100.KO.1.dmnd`; GUNC recibe su archivo `.dmnd`. La configuración filogenómica incluida es `assets/phylophlan_iqtree.cfg`, indicada mediante `phylophlan_config`.

#### Anotación funcional

```text
/ruta/referencias/eggnog/5.0.2/
/ruta/referencias/interproscan/5.59-91.0/data/
```

`eggnog_db` contiene `eggnog.db`, `eggnog.taxa.db` y `eggnog_proteins.dmnd`. `interproscan_data` señala los datos compatibles con InterProScan. GeneMarkS-2 requiere una instalación v1.15 y su licencia, configuradas mediante `genemark_home` y `genemark_key`.

### 4. Archivos auxiliares y configuraciones requeridas

Las rutas de las referencias se definen en una configuración externa:

```text
/ruta/referencias/metagenomics_databases.config
```

Este archivo se carga con `--database-config`. El repositorio incluye `nextflow.config`, `nextflow_schema.json`, las configuraciones de `conf/` y `assets/multiqc_config.yml`. La ejecución requiere **Nextflow ≥26.04.6** y el entorno seleccionado: Docker, Singularity, Apptainer o Conda.

## Estructura esperada de los datos

`Metagenomics/` contiene los archivos incluidos en GitHub; `/ruta/` representa los datos y referencias proporcionados en el entorno de ejecución. Los FASTQ locales se utilizan en lugar de la adquisición desde SRA.

```text
Metagenomics/
└── assets/
    ├── samplesheet.csv
    ├── raw/
    │   └── PRJNA782492_metagenomic_samples.tsv
    ├── phylophlan_iqtree.cfg
    └── multiqc_config.yml

/ruta/
├── datos/
│   ├── samplesheet.csv
│   └── fastq/
│       ├── <sample_id>_R1.fastq.gz
│       └── <sample_id>_R2.fastq.gz
│
└── referencias/
    ├── metagenomics_databases.config
    ├── human/
    │   ├── GRCh38.p14/GRCh38.p14.fa
    │   └── bowtie2/GRCh38_p14.*.bt2*
    ├── checkm2/1.1.0/
    ├── gunc/1.0.6/progenomes_2.1/
    ├── gtdbtk/release226/
    ├── phylophlan/phylophlan/
    ├── eggnog/5.0.2/
    ├── interproscan/5.59-91.0/data/
    └── genemarks2/
```

## Salidas principales

El pipeline genera:

- reportes FastQC de lecturas raw y procesadas.
- reportes de filtrado de fastp y lecturas filtradas.
- pares sin alineamiento humano concordante y logs de Bowtie2.
- ensamblajes mediante MEGAHIT y metaSPAdes.
- reportes de evaluación del ensamblaje mediante MetaQUAST.
- tablas de cobertura de contigs y archivos BAM mediante CoverM.
- bins de COMEBin, SemiBin2, MetaBAT2 y Vamb, refinados mediante DAS Tool.
- tablas de calidad de CheckM2 y evaluaciones de GUNC.
- MAGs seleccionados y desreplicados mediante dRep.
- clasificación taxonómica mediante GTDB-Tk y filogenia mediante PhyloPhlAn/IQ-TREE.
- proteínas, CDS y anotaciones GFF predichas mediante GeneMarkS-2.
- anotaciones de eggNOG-mapper e InterProScan, integradas por MAG.
- tablas de abundancia de MAGs mediante CoverM.
- reporte integrado de calidad y resultados mediante MultiQC.

Las lecturas filtradas de fastp, los BAM y otros productos intermedios se publican mediante `--save_clean_reads true`, `--save_bam true` y `--save_intermediates true`, respectivamente.

## Resultados obtenidos

El análisis metagenómico reportado en el trabajo de grado permitió obtener:

- **356 MAGs de alta calidad**, con completitud >90% y contaminación <5%.
- **98.43% de completitud promedio** y **0.73% de contaminación promedio**.
- **Siete filos bacterianos y dos arqueanos**.
- **194 MAGs de Bacillota** y **78 MAGs de Bacteroidota**.

Estos resultados son utilizados posteriormente en la etapa de integración multiómica para relacionar MAGs, genes, funciones biológicas, actividad transcripcional y metabolitos.

## Ejecución

La ejecución del pipeline se gestiona mediante:

[`metagenomics_pipeline.sh`](./metagenomics_pipeline.sh)

Este script ejecuta **Nextflow** con los perfiles local o HPC y el entorno seleccionado. En HPC, Nextflow envía las tareas al planificador **SLURM**, utilizando Singularity o Apptainer.

Desde el directorio `Metagenomics/`, se encuentran disponibles los siguientes modos de ejecución. Las rutas `/ruta/...` deben sustituirse por las del entorno de trabajo.

```bash
./metagenomics_pipeline.sh --local --docker --run \
  --input /ruta/datos/samplesheet.csv \
  --database-config /ruta/referencias/metagenomics_databases.config
```

Ejecuta el pipeline localmente con Docker y lecturas FASTQ.

```bash
./metagenomics_pipeline.sh --hpc --singularity --run \
  --sra-project PRJNA782492 \
  --sra-samples assets/raw/PRJNA782492_metagenomic_samples.tsv \
  --database-config /ruta/referencias/metagenomics_databases.config
```

Ejecuta el análisis de las muestras seleccionadas desde SRA mediante SLURM y Singularity.

```bash
./metagenomics_pipeline.sh --hpc --singularity --run --resume \
  --sra-project PRJNA782492 \
  --sra-samples assets/raw/PRJNA782492_metagenomic_samples.tsv \
  --database-config /ruta/referencias/metagenomics_databases.config
```

Reanuda una ejecución interrumpida utilizando `-resume` de Nextflow, desde el mismo directorio de lanzamiento y conservando la caché y el directorio de trabajo.

La preparación de referencias se realiza con `./metagenomics_pipeline.sh --local --prepare-databases --db-root /ruta/referencias`, proporcionando la instalación y clave de GeneMarkS-2 mediante `--genemark-home` y `--genemark-key`.

Antes de ejecutar el pipeline deben estar disponibles las lecturas o el manifiesto SRA, las referencias y las configuraciones descritas en la sección **Entradas**, respetando la estructura de directorios esperada.

## Documentación detallada

El flujo completo, incluyendo decisiones internas, filtros, criterios de selección y salidas intermedias, se encuentra en:

[`Flujo completo de metagenomica.pdf`](./Anexo_Metagenomica.pdf)

## Relación con otros módulos

Los MAGs generados por este pipeline sirven como referencia para el módulo de:

[`../Metatranscriptomics`](../Metatranscriptomics)

Las anotaciones funcionales y abundancias son posteriormente utilizadas en:

[`../Integration`](../Integration)

para la construcción de asociaciones funcionales **microbio–gen–metabolito**.
