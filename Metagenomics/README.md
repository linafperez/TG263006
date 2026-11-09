# Metagenómica

Pipeline metagenómico modular implementado en **Nextflow** para caracterizar la composición taxonómica y el potencial funcional del microbioma intestinal en pacientes con enfermedad de Parkinson (PD) y controles sanos (HC).

## Esquema general

[`Diagrama general metagenomica.pdf`](./diagrama_metagenomica.pdf)

El flujo comprende las etapas de **obtención de lecturas**, **preprocesamiento y control de calidad**, **remoción de secuencias humanas**, **ensamblaje metagenómico**, **binning y refinamiento de MAGs**, **clasificación taxonómica y análisis filogenómico**, **predicción de genes y anotación funcional**, **estimación de abundancia** y **evaluación global del procesamiento**.

## Software y requerimientos computacionales

Los valores corresponden a **solicitudes por tarea** de [`conf/resources.config`](./conf/resources.config), con versiones fijadas en los módulos; no son consumos observados para 95 muestras. **N/D** indica ausencia de mediciones verificables. Las asignaciones históricas de [`docs/original_resource_mapping.tsv`](./docs/original_resource_mapping.tsv) tampoco constituyen mediciones.

| Software y versión | Módulo / etapa | CPU | RAM asignada (GB) | RAM pico (GB) | Almacenamiento | Tiempo |
|---|---|---:|---:|---:|---:|---:|
| SRA Toolkit v3.4.1 | Adquisición por corrida SRA | 12 | 19 | N/D | N/D | N/D |
| pigz v2.8 | Compresión dentro de la misma tarea SRA | 12 | 19 | N/D | N/D | N/D |
| FastQC v0.12.1 | Control de calidad raw y clean | 4 | 19 | N/D | N/D | N/D |
| fastp v1.0.1 | Filtrado y trimming | 16 | 19 | N/D | N/D | N/D |
| Bowtie2 v2.5.4 | Remoción de lecturas humanas | 24 | 19 | N/D | N/D | N/D |
| MEGAHIT v1.2.9 | Coensamblaje global | 32 | 500 | N/D | N/D | N/D |
| metaSPAdes v4.2.0 | Coensamblaje global o por condición | 32 | 1800 | N/D | N/D | N/D |
| MetaQUAST v5.3.0 | Evaluación de ambos ensamblajes | 6 | 20 | N/D | N/D | N/D |
| CoverM v0.7.0 | Cobertura de contigs | 24 | 100 | N/D | N/D | N/D |
| COMEBin v1.0.4 | Binning | 32 | 100 | N/D | N/D | N/D |
| SemiBin2 v1.5.0 | Binning | 25 | 100 | N/D | N/D | N/D |
| MetaBAT2 v2.18 | Binning | 32 | 100 | N/D | N/D | N/D |
| Vamb v5.0.4 | Binning de MEGAHIT | 32 | 100 | N/D | N/D | N/D |
| Vamb v5.0.4 | Binning de metaSPAdes | 12 | 20 | N/D | N/D | N/D |
| DAS Tool v1.1.7 | Integración y refinamiento de bins | 24 | 250 | N/D | N/D | N/D |
| CheckM2 v1.1.0 | Calidad de bins y MAGs | 24 | 25 | N/D | N/D | N/D |
| GUNC v1.0.6 | Evaluación de contaminación entre linajes | 24 | 25 | N/D | N/D | N/D |
| dRep v3.6.2 | Desreplicación de MEGAHIT y catálogo final | 24 | 100 | N/D | N/D | N/D |
| dRep v3.6.2 | Desreplicación de metaSPAdes | 24 | 50 | N/D | N/D | N/D |
| GTDB-Tk v2.6.1 | Clasificación taxonómica | 24 | 100 | N/D | N/D | N/D |
| PhyloPhlAn v3.1.1 / IQ-TREE v3.0.1 | Análisis filogenómico | 32 | 120 | N/D | N/D | N/D |
| GeneMarkS-2 v1.15 | Predicción de genes por MAG | 1 | 20 | N/D | N/D | N/D |
| eggNOG-mapper v2.1.13 | Anotación funcional por MAG | 5 | 30 | N/D | N/D | N/D |
| InterProScan v5.59-91.0 | Anotación funcional por MAG | 16 | 128 | N/D | N/D | N/D |
| CoverM v0.7.0 | Abundancia del catálogo final | 24 | 100 | N/D | N/D | N/D |
| MultiQC v1.35 | Reporte y control de calidad global | 2 | 8 | N/D | N/D | N/D |

Las mayores **solicitudes de memoria** corresponden a metaSPAdes, MEGAHIT y DAS Tool. Los límites máximos por tarea son `max_cpus = 32`, `max_memory = '1800 GB'` y `max_time = '14d'`. SRA Toolkit y pigz comparten asignación; fastp utiliza seis hilos aunque solicita 16 CPU. Se requiere **Nextflow ≥26.04.6**, Java compatible y el entorno de ejecución seleccionado. Estas asignaciones justifican utilizar infraestructura HPC. El perfil `hpc` desactiva reportes, trazas y cronogramas de recursos; los consumos deben verificarse con registros del planificador u otra monitorización explícita.

## Dataset

| Característica | Descripción |
|---|---|
| BioProject | PRJNA782492 |
| Muestras | 95 |
| Parkinson (PD) | 46 |
| Controles (HC) | 49 |
| Tipo de datos | Secuenciación metagenómica shotgun paired-end |
| Fuente original | Sequence Read Archive (SRA) |

El manifiesto metagenómico identifica los controles mediante **`CT`**, equivalente al grupo **HC** descrito aquí. Los valores originales `PD` y `CT` se conservan durante el procesamiento.

## Entradas

El pipeline requiere cuatro grupos principales de archivos de entrada.

### 1. Lecturas metagenómicas

El modo local recibe archivos FASTQ paired-end, comprimidos o sin comprimir, mediante `--input`. Las rutas relativas del samplesheet se resuelven respecto al directorio del CSV. Por ejemplo, el usuario puede proporcionar:

```text
/ruta/datos/fastq/
```

Cada muestra debe disponer de una lectura forward y otra reverse:

```text
<sample_id>_R1.fastq.gz
<sample_id>_R2.fastq.gz
```

El modo SRA utiliza `--sra-project PRJNA782492` y `--sra-samples` para resolver corridas metagenómicas paired-end del manifiesto, adquirir sus lecturas y combinar las del mismo BioSample. Los modos son excluyentes; los FASTQ no están incluidos en GitHub.

### 2. Samplesheet y manifiesto de muestras

El manifiesto completo de selección SRA se encuentra en:

```text
assets/raw/PRJNA782492_metagenomic_samples.tsv
```

Contiene 95 BioSamples y dos columnas separadas por tabuladores. Este fragmento reproduce identificadores y grupos reales:

```tsv
biosample_accession	group
SAMN29531760	PD
SAMN29531788	CT
```

El archivo `assets/raw/PRJNA782492_subset8.tsv` selecciona ocho BioSamples, cuatro por grupo, sin acreditar una ejecución terminada.

Para lecturas locales se proporciona la plantilla:

```text
assets/samplesheet.csv
```

Su contenido es ilustrativo y debe reemplazarse por rutas reales:

```csv
sample,fastq_1,fastq_2
sample_A,/path/to/sample_A_R1.fastq.gz,/path/to/sample_A_R2.fastq.gz
```

| Columna | Descripción |
|---|---|
| `biosample_accession` | TSV SRA: accesión BioSample para resolver las corridas |
| `group` | TSV SRA y CSV local opcional: grupo experimental `PD` o `CT` |
| `sample` | Identificador único de la muestra local |
| `fastq_1` | Ruta al archivo forward |
| `fastq_2` | Ruta al archivo reverse |

La plantilla contiene una muestra ficticia, no las 95 muestras. Admite una columna `group`, detectada automáticamente, u otra indicada con `--group-column`. metaSPAdes realiza coensamblaje global por defecto; por condición requiere `--spades-coassembly-mode condition --group-column group`. `--sra-run-info` admite un CSV RunInfo existente.

### 3. Referencias y bases de datos

Las referencias deben estar disponibles en el entorno de ejecución. Con una raíz ilustrativa `/ruta/referencias`, `bin/prepare_databases.sh` utiliza las siguientes rutas y genera `metagenomics_databases.config`.

#### Referencia humana

```text
/ruta/referencias/human/GRCh38.p14/GRCh38.p14.fa
/ruta/referencias/human/bowtie2/GRCh38_p14
```

El FASTA corresponde a **GRCh38.p14 de GENCODE v44**. `host_bowtie2_index` recibe el **prefijo** del índice con sus seis archivos `.bt2` o `.bt2l`. Bowtie2 conserva mediante `--un-conc-gz` los pares sin alineamiento concordante, considerados potencialmente microbianos.

#### Evaluación de MAGs

```text
/ruta/referencias/checkm2/1.1.0/
/ruta/referencias/gunc/1.0.6/progenomes_2.1/
```

`checkm2_db` admite el archivo DIAMOND o su directorio; el preparador valida `uniref100.KO.1.dmnd`. `gunc_db` recibe el archivo `.dmnd`. La preparación vigente selecciona **ProGenomes 2.1**, mientras que el anexo representa GTDB95; la reproducción requiere conservar la base realmente utilizada.

#### Taxonomía y filogenómica

```text
/ruta/referencias/gtdbtk/release226/
/ruta/referencias/phylophlan/phylophlan/
```

`gtdbtk_db` recibe **GTDB release 226** y `phylophlan_db` el directorio de marcadores. `phylophlan_config` utiliza por defecto `assets/phylophlan_iqtree.cfg`, junto con PhyloPhlAn e IQ-TREE en un mismo entorno.

#### Anotación funcional

```text
/ruta/referencias/eggnog/5.0.2/
/ruta/referencias/interproscan/5.59-91.0/data/
```

`eggnog_db` contiene `eggnog.db`, `eggnog.taxa.db` y `eggnog_proteins.dmnd`. El preparador utiliza **eggNOG 5.0.2** para eggNOG-mapper v2.1.13, frente a v6.0 del anexo. `interproscan_data` recibe los datos de **InterProScan v5.59-91.0**.

### 4. Archivos auxiliares y configuraciones requeridas

GeneMarkS-2 requiere una instalación externa v1.15 mediante `genemark_home` y su clave de licencia mediante `genemark_key`. Estos archivos se proporcionan en el entorno de ejecución. El repositorio incluye las configuraciones de recursos, módulos y planificador:

```text
nextflow.config
nextflow_schema.json
conf/resources.config
conf/modules.config
conf/hpc.config
assets/phylophlan_iqtree.cfg
assets/multiqc_config.yml
```

Las rutas de referencias se cargan con `--database-config /ruta/referencias/metagenomics_databases.config`, o mediante los parámetros Nextflow correspondientes. Los parámetros opcionales `interproscan_fingerprintscan_path`, `interproscan_cath_resolve_hits_path` e `interproscan_ps_path` permiten especificar ejecutables auxiliares según el entorno.

La configuración aplica fastp con Phred 30 y longitud mínima de 30 bases, y contigs de al menos `min_contig_length = 1500` bases para binning. La selección exige **completitud >90% y contaminación <5%** según CheckM2; GUNC evalúa los bins sin filtrarlos automáticamente. Se desreplica al **99% ANI**, con cobertura mínima de comparación del **30%**, por rama y tras combinar ambos ensambladores; se generan representantes al **95% ANI** por separado. El anexo representa 99.5% y desigualdades inclusivas, diferentes del código vigente.

## Estructura esperada de los datos

El siguiente árbol distingue los archivos incluidos en GitHub de los datos externos. Las rutas `/ruta/datos` y `/ruta/referencias` son ejemplos; en modo SRA las lecturas se adquieren dentro del directorio de trabajo de Nextflow.

```text
Metagenomics/                                  # Incluido en GitHub
├── metagenomics_pipeline.sh
├── main.nf
├── nextflow.config
├── nextflow_schema.json
├── conf/
├── assets/
│   ├── samplesheet.csv                        # Plantilla local
│   ├── raw/
│   │   ├── PRJNA782492_metagenomic_samples.tsv # 95 BioSamples
│   │   └── PRJNA782492_subset8.tsv             # 8 BioSamples
│   ├── phylophlan_iqtree.cfg
│   └── multiqc_config.yml
├── diagrama_metagenomica.pdf
└── Anexo_Metagenomica.pdf

/ruta/datos/                                  # Externo; modo local
├── samplesheet.csv                           # Rutas reales del usuario
└── fastq/
    ├── <sample_id>_R1.fastq.gz
    └── <sample_id>_R2.fastq.gz

/ruta/referencias/                            # Externo; ambos modos
├── metagenomics_databases.config
├── database_manifest.tsv
├── human/
│   ├── GRCh38.p14/GRCh38.p14.fa
│   └── bowtie2/GRCh38_p14.*.bt2*              # Seis archivos de índice
├── checkm2/1.1.0/                             # Incluye uniref100.KO.1.dmnd
├── gunc/1.0.6/progenomes_2.1/                 # Contiene la base .dmnd
├── gtdbtk/release226/
├── phylophlan/phylophlan/
├── eggnog/5.0.2/
│   ├── eggnog.db
│   ├── eggnog.taxa.db
│   └── eggnog_proteins.dmnd
├── interproscan/5.59-91.0/data/
└── genemarks2/                               # Instalación y clave externas
```

## Salidas principales

El pipeline genera:

- reportes FastQC de lecturas raw y procesadas, y reportes HTML/JSON de fastp.
- lecturas filtradas y pares sin alineamiento humano concordante, con logs de Bowtie2.
- coensamblajes de MEGAHIT y metaSPAdes, y métricas TSV/HTML de MetaQUAST.
- tablas de cobertura de contigs, abundancia para Vamb y archivos BAM de CoverM.
- bins de COMEBin, SemiBin2, MetaBAT2 y Vamb, asignaciones contig–bin y refinamiento de DAS Tool.
- tablas de calidad CheckM2, evaluaciones GUNC y selección de MAGs de alta calidad.
- representantes y clústeres de dRep al 99% y 95% ANI, y catálogo final con trazabilidad y calidad.
- clasificación taxonómica bacteriana y arqueana mediante GTDB-Tk.
- alineamiento concatenado y árbol filogenómico Newick mediante PhyloPhlAn/IQ-TREE.
- proteínas, CDS y anotaciones GFF predichas mediante GeneMarkS-2.
- anotaciones de eggNOG-mapper e InterProScan, e integración funcional TSV por MAG.
- tablas de abundancia de MAGs en formatos ancho y largo mediante CoverM.
- reporte MultiQC, versiones de software y metadata de muestras; en modo SRA, metadata de adquisición y exclusiones.

Los productos se organizan bajo `outdir` en `01_quality_control_and_filtering/`, `02_mag_construction/`, `03_taxonomic_classification_and_phylogenomics/`, `04_gene_prediction_and_functional_annotation/`, `05_mag_abundance_estimation/`, `06_global_processing_evaluation/` y `pipeline_info/`. Las lecturas de fastp se publican con `--save_clean_reads true`; los pares no concordantes con el hospedero se publican por defecto. Los BAM y otros productos intermedios requieren `--save_bam true` y `--save_intermediates true`, respectivamente. MultiQC integra los reportes compatibles disponibles.

## Resultados obtenidos

El trabajo de grado **TG263006**, en el apartado Metagenómica y la Figura 6, reporta para el estudio del BioProject PRJNA782492:

- **356 MAGs de alta calidad**, tras desreplicación al **99% ANI**.
- **98.43% de completitud promedio** y **0.73% de contaminación promedio**.
- **Siete filos bacterianos y dos arqueanos**, con predominio de **Bacillota (194 MAGs)** y **Bacteroidota (78 MAGs)**.

La tesis describe una selección de 95 muestras y criterios de completitud >90% y contaminación <5%. Estas cifras son **resultados previos reportados**, sin catálogo ni registros en el repositorio que los vinculen a una ejecución concreta de esta versión. Los manifiestos de 95 y ocho muestras describen entradas; no permiten certificar la finalización del piloto, del preprocesamiento de 95 muestras ni de su análisis completo. Los 356 MAGs no se atribuyen a una ejecución productiva nueva ni al subconjunto de ocho muestras.

Los MAGs, sus anotaciones funcionales y sus abundancias permiten relacionar genomas microbianos, genes y funciones biológicas con la actividad transcripcional y los metabolitos durante la integración multiómica posterior.

## Ejecución

La ejecución del pipeline se gestiona mediante:

[`metagenomics_pipeline.sh`](./metagenomics_pipeline.sh)

Este script lanza **una ejecución de Nextflow**: `--local` utiliza el ejecutor local y `--hpc` envía tareas a **SLURM**. `--docker` selecciona Docker para ejecución local y `--singularity` selecciona Singularity, también en HPC. Admite además `--conda` y `--apptainer`; Docker no es compatible con `--hpc`.

Desde el directorio `Metagenomics/`, pueden utilizarse los siguientes modos. Todas las rutas `/ruta/...` son ilustrativas y deben sustituirse por archivos y directorios del entorno de ejecución.

```bash
./metagenomics_pipeline.sh --local --docker --run \
  --input /ruta/datos/samplesheet.csv \
  --database-config /ruta/referencias/metagenomics_databases.config \
  --outdir /ruta/resultados/local --work-dir /ruta/trabajo/local
```

Ejecuta el pipeline con lecturas locales y un samplesheet real. La plantilla de GitHub debe adaptarse antes de utilizarla.

```bash
./metagenomics_pipeline.sh --hpc --singularity --run \
  --sra-project PRJNA782492 \
  --sra-samples assets/raw/PRJNA782492_metagenomic_samples.tsv \
  --database-config /ruta/referencias/metagenomics_databases.config \
  --slurm_queue particion_del_cluster \
  --outdir /ruta/resultados/PRJNA782492 --work-dir /ruta/trabajo/PRJNA782492
```

Ejecuta la selección SRA mediante SLURM y Singularity. La partición y, cuando correspondan, `--slurm_account` y `--slurm_qos` deben ajustarse al clúster. El proceso coordinador de Nextflow debe mantenerse activo o ejecutarse dentro de una asignación apropiada del clúster.

```bash
./metagenomics_pipeline.sh --hpc --singularity --run --resume \
  --sra-project PRJNA782492 \
  --sra-samples assets/raw/PRJNA782492_metagenomic_samples.tsv \
  --database-config /ruta/referencias/metagenomics_databases.config \
  --slurm_queue particion_del_cluster \
  --outdir /ruta/resultados/PRJNA782492 --work-dir /ruta/trabajo/PRJNA782492
```

Reanuda con `-resume` de Nextflow desde el mismo directorio de lanzamiento, conservando el directorio de trabajo, la caché y las entradas de la ejecución anterior.

```bash
./metagenomics_pipeline.sh --local --prepare-databases \
  --db-root /ruta/referencias \
  --genemark-home /ruta/referencias/genemarks2 \
  --genemark-key /ruta/licencias/clave_genemark
```

Prepara y valida las referencias mediante `bin/prepare_databases.sh`, y genera el manifiesto de bases de datos y su configuración Nextflow. Este modo no inicia el workflow, no recibe una opción de contenedor y requiere las herramientas de preparación en el entorno, además de la instalación y licencia externas de GeneMarkS-2. `--dry-run` permite revisar el comando del lanzador sin iniciar la ejecución.

Antes de ejecutar deben estar disponibles las lecturas o el manifiesto SRA, las referencias y las configuraciones descritas en **Entradas**.

## Documentación detallada

El flujo completo, incluyendo decisiones internas, herramientas, filtros, procedimientos y productos intermedios, se encuentra en:

[`Flujo completo de metagenomica.pdf`](./Anexo_Metagenomica.pdf)

Los diagramas se conservan en su versión original. Las diferencias respecto a la implementación vigente se explican en **Entradas**, especialmente para las bases GUNC y eggNOG y los criterios de selección y desreplicación de MAGs.

## Relación con otros módulos

Los MAGs generados por este pipeline sirven como referencia genómica para el módulo de:

[`../Metatranscriptomics`](../Metatranscriptomics)

Las anotaciones funcionales y abundancias son posteriormente utilizadas en:

[`../Integration`](../Integration)

para la construcción de asociaciones funcionales **microbio–gen–metabolito**.
