# Aleatorización mendeliana

Análisis de aleatorización mendeliana de dos muestras desarrollado para evaluar posibles relaciones causales entre la abundancia de taxones del microbioma intestinal y el riesgo de enfermedad de Parkinson (PD), utilizando asociaciones entre variantes genéticas humanas y estas características provenientes de estudios GWAS.

## Esquema general

[`Diagrama general de aleatorización mendeliana.pdf`](./diagrama_aleatorizacion_mendeliana.pdf)

El análisis comprende la **preparación de exposición y desenlace**, **selección y armonización de SNPs**, **anotación genética**, **estimación causal**, **análisis de sensibilidad** y **visualización de resultados**.

## Software y librerías

| Software / librería | Función |
|---|---|
| R v4.5.2 / RStudio v2025.09.2+418 | Ejecución de scripts y documentos R Markdown |
| data.table | Lectura y escritura de estadísticas resumidas de GWAS |
| tidyverse / dplyr | Estandarización, filtrado y combinación de datos |
| PLINK2 | Selección de variantes mediante desequilibrio de ligamiento |
| biomaRt / Ensembl | Anotación de SNPs, genes humanos y consecuencias funcionales |
| TwoSampleMR | Estimación causal, odds ratios y análisis de sensibilidad |
| MR-PRESSO | Evaluación de pleiotropía e identificación de variantes atípicas |
| ggplot2 | Generación y exportación de gráficos |
| forestploter | Visualización de estimaciones mediante forest plots |


## Fuentes de datos

El análisis utiliza estadísticas resumidas de fuentes GWAS independientes para la exposición y el desenlace.

| Componente | Fuente | Información |
|---|---|---|
| Exposición: microbioma intestinal | Consorcio MiBioGen ([Kurilshikov et al., 2021](https://doi.org/10.1038/s41588-020-00763-1)) | Asociaciones entre variantes humanas y abundancia relativa de 211 taxones; 18 340 participantes de 24 cohortes |
| Desenlace: Parkinson en EUR | Global Parkinson's Genetics Program, GP2 ([Leonard et al., 2025](https://doi.org/10.1101/2025.03.14.24319455)) | Asociaciones genéticas con Parkinson en poblaciones europeas |
| Desenlace: Parkinson en EUR | FinnGen ([Kurki et al., 2023](https://doi.org/10.1038/s41586-022-05473-8)) | Asociaciones genéticas con Parkinson en población finlandesa |
| Desenlace: Parkinson en LATAM | Latin American Research Consortium on the Genetics of Parkinson's Disease, LARGE-PD ([Loesch et al., 2021](https://doi.org/10.1002/ana.26153)) | Asociaciones genéticas con Parkinson en población latinoamericana |

## Entradas

### 1. Datos de exposición

Los archivos de MiBioGen deben estar disponibles en `unzipped_mibiogen/`, dentro del directorio local de exposición, con el patrón:

```text
*.summary.txt.gz
```

Contienen identificadores SNP (`rsID`), alelos (`ref.allele`, `eff.allele`), efecto (`beta`), error estándar (`SE`), valor p (`P.weightedSumZ`), tamaño muestral (`N`) y taxón (`bac`). El script procesa los archivos por género bacteriano.

### 2. Datos de desenlace

`SNP_outcome.R` utiliza los siguientes archivos locales:

```text
GWAS_EUR-GP2_rsid.txt.gz
GWAS_EUR-FINN.gz
GWAS_Latin_rsid.tsv.gz
```

Estos contienen las asociaciones entre SNPs y Parkinson, incluyendo beta, error estándar, alelos, valor p, posición genómica y frecuencia del alelo de efecto (EAF).

### 3. Paneles de referencia genética

PLINK2 utiliza paneles externos del **Proyecto 1000 Genomas** ([The 1000 Genomes Project Consortium, 2015](https://doi.org/10.1038/nature15393)) con archivos `.pgen`, `.pvar` y `.psam`, organizados bajo el directorio `Genotype/` de cada fuente:

```text
Genotype/<población>/1KG_<población>.pgen
Genotype/<población>/1KG_<población>.pvar
Genotype/<población>/1KG_<población>.psam
```

La exposición utiliza **EUR, SAS, EAS, AFR y AMR**. Para el desenlace se emplean **EUR** en los archivos europeos y **AMR** en el latinoamericano.

## Preparación de exposición y desenlace

Ambas fuentes se procesan por separado mediante **data.table** y **tidyverse**. Los scripts estandarizan los nombres de columnas, eliminan SNPs duplicados y conservan las variables de asociación en los formatos de exposición y desenlace utilizados posteriormente por TwoSampleMR.

## Selección y armonización de instrumentos genéticos

La selección aplica los siguientes criterios:

- **Significancia:** p < 1 × 10⁻⁵.
- **Desequilibrio de ligamiento:** r² < 0.001, con una ventana de 10 000 kb.
- **Variantes palindrómicas:** exclusión de pares A/T y C/G, en ambas orientaciones.
- **Fortaleza de instrumentos:** F ≥ 10, calculada como F = (β/SE)².
- **Frecuencia en el desenlace:** EAF disponible y superior a 0.01.

Estos criterios seleccionan variantes asociadas, reducen la redundancia genética y descartan instrumentos débiles o con ambigüedad alélica. El código aplica los filtros de significancia, LD y F por separado a exposición y desenlace.

La armonización alinea los alelos de efecto entre ambas fuentes. `MR.Rmd` requiere el archivo **`PD_harmonised.txt` preparado previamente**; su generación no está implementada en los tres archivos disponibles. Posteriormente, **biomaRt / Ensembl** permiten anotar posiciones, consecuencias funcionales e identificadores de genes humanos asociados a los SNPs.

## Estimación causal

Las estimaciones se realizan mediante **TwoSampleMR**, utilizando cuatro métodos:

- Varianza inversa ponderada (**IVW**).
- Mediana ponderada.
- Moda ponderada.
- Razón de Wald, cuando se dispone de un único SNP.

Los efectos se expresan como **odds ratios (OR)** e **intervalos de confianza del 95 %**, mediante `generate_odds_ratios()`.

## Análisis de sensibilidad

La robustez de las estimaciones se examina mediante:

- **Q de Cochran:** evaluación de heterogeneidad entre instrumentos.
- **MR-Egger:** evaluación de pleiotropía horizontal direccional mediante el intercepto.
- **MR-PRESSO:** identificación de variantes atípicas, con **1 000 simulaciones** y umbral de significancia de **0.05**.
- **Leave-one-out:** evaluación de la influencia de cada SNP al excluirlo individualmente.


## Salidas principales

Los scripts y el documento R Markdown definen las siguientes salidas:

- instrumentos intermedios por fuente: archivos `instruments_*` y `all_instruments_combined.txt` en los respectivos directorios `results/`.
- anotaciones intermedias de SNPs: `PD_data.tsv`.
- estimaciones con OR e intervalos de confianza: `MR_results/PD_MR_results.csv`.
- heterogeneidad: `Sensitivity_analysis/Heterogeneity/heterogeneity_results.csv`.
- pleiotropía horizontal: `Sensitivity_analysis/Pleiotropy/pleiotropy_results.csv`.
- MR-PRESSO: `Sensitivity_analysis/MR-PRESSO/mr-presso_results.csv`.
- leave-one-out: `Sensitivity_analysis/Leave-one-out/leave_one_out_results.csv` y `leave_one_out_plot.pdf` en ese mismo directorio.
- gráficos de dispersión, embudo y forest: `MR_results/MR_scatter_plot.pdf`, `MR_results/funnel_plot.pdf` y `MR_results/forest_plot.pdf`.
- forest plots por método y SNP: `Forest/forest_plot.pdf` y `Forest/pd_snps_forest_plot.pdf`.

## Resultados obtenidos

El trabajo de grado reportó evidencia compatible con un posible efecto protector del género **Bifidobacterium** frente al riesgo de Parkinson:

- **IVW, cinco SNPs:** OR = 0,51; IC del 95 %: 0,44–0,59; p = 1,40 × 10⁻¹⁸.
- **Mediana ponderada:** OR = 0,55; p = 2,70 × 10⁻⁸.
- **Moda ponderada:** OR = 0,56; p = 0,004.
- **Razón de Wald, un SNP:** OR = 0,46; p = 1,09 × 10⁻⁷.

Las pruebas de sensibilidad reportaron **MR-Egger p = 0,075** y **MR-PRESSO p = 0,507**. Estos resultados aportan evidencia compatible con una asociación protectora, condicionada a los supuestos de la aleatorización mendeliana; no demuestran causalidad. Los instrumentos corresponden a SNPs humanos asociados con taxones intestinales.

## Ejecución

El análisis se ejecuta localmente mediante **R y R Markdown**, en el siguiente orden:

1. [`SNPs_exposure.R`](./Exposure/SNPs_exposure.R): preparación y selección de instrumentos de exposición.
2. [`SNP_outcome.R`](./Outcome/SNP_outcome.R): preparación y filtrado de variantes del desenlace.
3. Preparación de `PD_harmonised.txt` con las asociaciones y los alelos alineados.
4. [`MR.Rmd`](./MR/MR.Rmd): anotación, estimación causal, sensibilidad y visualización.

Antes de ejecutar deben ajustarse las rutas absolutas `base_dir` y las referencias a la ubicación local de los archivos. El documento exporta `PD_data.tsv` y después lee **`PD_data-annotated.tsv`**, que debe estar disponible con ese nombre. También deben existir los directorios de salida indicados en el documento R Markdown.

## Documentación detallada

El flujo completo, incluyendo preparación de datos GWAS, selección de instrumentos, anotación genética, estimación causal, análisis de sensibilidad y visualización, se encuentra en:

[`Flujo completo de aleatorización mendeliana.pdf`](./Anexo_Aleatorizacion_Mendeliana.pdf)
