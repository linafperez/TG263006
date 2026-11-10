# Metabolómica

Workflow computacional para el procesamiento, selección y consolidación de señales de **metabolómica no dirigida** asociadas con enfermedad de Parkinson (PD) y controles sanos (HC).

El análisis combina **MS-DIAL v5**, **MetaboAnalyst v6.0** y un flujo reproducible desarrollado en **R Markdown**. Las tres plataformas analíticas —LC-MS en modo positivo, LC-MS en modo negativo y GC-MS— se procesan inicialmente de forma independiente y convergen posteriormente mediante una estrategia de **consenso algorítmico** entre ANOVA, sPLS-DA y Random Forest.

## Esquema general

[`Diagrama general metabolómica.pdf`](https://github.com/linafperez/TG263006/blob/main/Metabolomics/diagrama_metabolomica.pdf)

El workflow comprende las etapas de **procesamiento y anotación de espectros**, **control de calidad y preprocesamiento**, **análisis estadístico univariado**, **análisis multivariado**, **clasificación mediante aprendizaje automático**, **selección por consenso** y **consolidación de metabolitos entre plataformas**.

De forma general:

```text
Datos metabolómicos
│
├── LC-MS positivo
├── LC-MS negativo
└── GC-MS
        │
        ▼
     MS-DIAL v5
        │
        ├── detección de picos
        ├── deconvolución
        ├── alineamiento
        └── anotación
        │
        ▼
 MetaboAnalyst v6.0
        │
        ├── control de calidad
        ├── normalización
        ├── ANOVA
        ├── sPLS-DA
        └── Random Forest
        │
        ▼
      R Markdown
        │
        ├── selección top 5% de Random Forest
        ├── intersección ANOVA / sPLS-DA / RF
        ├── diagramas de Venn
        ├── clasificación PD / Healthy
        ├── consolidación entre plataformas
        └── UpSet plot
        │
        ▼
 Firma metabolómica consolidada
```

A diferencia de los módulos metagenómico, metatranscriptómico y 16S, este workflow **no fue implementado en Nextflow**. La etapa final de procesamiento y consolidación se encuentra automatizada mediante un [`Script de R Markdown`](https://github.com/linafperez/TG263006/blob/main/Metabolomics/Metabolomics_analysis.Rmd), mientras que MS-DIAL y MetaboAnalyst constituyen etapas externas previas.

---

## Software

El análisis fue desarrollado utilizando:

- **MS-DIAL v5**
- **MetaboAnalyst v6.0**
- **R**
- **R Markdown**

El código R utiliza los siguientes paquetes:

```text
tidyverse
ggplot2
ggrepel
readr
ggvenn
UpSetR
```

El archivo R Markdown comprueba la disponibilidad de los paquetes CRAN requeridos e instala automáticamente aquellos que no estén presentes.

La versión de R y las versiones de los paquetes utilizadas en cada ejecución pueden consultarse mediante la salida de:

```r
sessionInfo()
```

incluida al final del documento R Markdown.

---

## Alcance del código R Markdown

El código disponible en este repositorio (https://github.com/linafperez/TG263006/blob/main/Metabolomics/Metabolomics_analysis.Rmd) corresponde principalmente a la etapa de **posprocesamiento, consenso y consolidación** del workflow metabolómico.

El R Markdown:

- procesa la tabla de alineamiento de GC-MS exportada desde MS-DIAL;
- elimina anotaciones desconocidas;
- calcula el RSD de las muestras QC;
- conserva características con `RSD < 80%`;
- importa los resultados estadísticos generados previamente mediante MetaboAnalyst;
- selecciona las características con mayor importancia en Random Forest;
- intersecta los resultados de ANOVA, sPLS-DA y Random Forest;
- genera diagramas de Venn;
- asigna la dirección de enriquecimiento a PD o HC;
- integra los metabolitos consenso de LC-MS(+), LC-MS(-) y GC-MS;
- genera archivos finales separados por grupo;
- genera un UpSet plot para comparar las plataformas.

Las etapas de detección de picos, deconvolución, alineamiento y anotación mediante MS-DIAL, así como el cálculo principal de ANOVA, sPLS-DA, Random Forest y el preprocesamiento realizado mediante MetaboAnalyst, deben realizarse previamente.

---

## Requerimiento de selección por consenso

Para reducir la incorporación de potenciales falsos positivos, el módulo establece como requerimiento que una característica metabólica sea retenida únicamente cuando esté respaldada por **tres aproximaciones estadísticas independientes**.

### 1. ANOVA

Las características deben cumplir:

```text
|log2FC| > 1.5
FDR ajustado < 0.05
```

### 2. sPLS-DA

Las características deben estar respaldadas por el modelo multivariado **sPLS-DA**, cuya estabilidad se evalúa mediante:

```text
1,000 permutaciones
```

### 3. Random Forest

Las características deben encontrarse dentro del:

```text
5% de mayor importancia
```

según `MeanDecreaseAccuracy`.

Solo los metabolitos presentes simultáneamente en los tres conjuntos son incorporados a la firma metabolómica de cada plataforma.

Conceptualmente:

```text
Consenso =
ANOVA ∩ sPLS-DA ∩ Random Forest
```

El criterio se aplica independientemente a:

- LC-MS positivo;
- LC-MS negativo;
- GC-MS.

---

## Dataset de referencia

El workflow se aplicó a la cohorte de Pereira et al., conformada por muestras séricas de pacientes con enfermedad de Parkinson y controles sanos.

| **Característica** | **Descripción** |
| --- | --- |
| Muestras | 124 |
| Parkinson (PD) | 63 |
| Controles sanos (HC) | 61 |
| Tipo de muestra | Suero |
| Diseño metabolómico | No dirigido |
| Plataformas | UHPLC-MS(+), UHPLC-MS(-) y GC-MS |
| Publicación de referencia | [Pereira et al. (2022)](https://doi.org/10.1038/s41531-022-00300-3) |

Las tres plataformas fueron utilizadas de forma complementaria para ampliar la cobertura química del metaboloma.

Los archivos metabolómicos originales y las matrices de resultados **no se distribuyen dentro de este repositorio**.

---

## Procesamiento mediante MS-DIAL

Los datos obtenidos mediante cada plataforma se procesan inicialmente utilizando **MS-DIAL v5**.

El procesamiento comprende:

1. detección de picos;
2. deconvolución espectral;
3. alineamiento por tiempo de retención;
4. anotación de metabolitos.

La anotación utiliza librerías espectrales compatibles con la plataforma analítica, incluyendo recursos como:

- MassBank;
- GNPS;
- ReSpect;
- RIKEN PlaSMA;
- LipidBlast;
- Fiehn BinBase;
- MoNA.

Como resultado se obtiene una matriz de intensidades de características metabólicas por muestra para cada plataforma.

---

## Control de calidad y preprocesamiento

Las matrices obtenidas mediante MS-DIAL se procesan posteriormente mediante **MetaboAnalyst v6.0**.

El workflow metodológico incluye:

- evaluación de la variabilidad de las muestras QC mediante RSD;
- filtrado por IQR;
- normalización por suma;
- transformación log10;
- escalamiento Pareto.

En el código R Markdown se implementa explícitamente para GC-MS el cálculo:

```r
RSD = (SD / media) × 100
```

y se conservan las características que cumplen:

```text
RSD < 80%
```

Las columnas correspondientes a muestras QC son eliminadas posteriormente antes de construir la matriz utilizada para los análisis estadísticos.

---

## Convención de nombres de muestras

El código R Markdown identifica automáticamente los grupos a partir del prefijo de las columnas:

| **Prefijo** | **Grupo** |
| --- | --- |
| `CT` | Healthy / control sano |
| `PD` | Parkinson |
| `QC` | Quality Control |

Por tanto, las matrices utilizadas deben conservar esta convención para permitir la identificación automática de los grupos.

---

## Análisis estadísticos

Después del preprocesamiento, las tres plataformas se analizan independientemente.

### ANOVA

El análisis univariado identifica características con cambios significativos entre PD y HC.

Los archivos utilizados posteriormente por el R Markdown se esperan en forma de:

```text
volcano.csv
```

y deben corresponder a las características que cumplen los criterios definidos para el análisis.

### sPLS-DA

El análisis multivariado utiliza **sPLS-DA** para identificar combinaciones de características que contribuyen a separar los grupos.

El R Markdown utiliza:

```text
splsda_loadings.csv
```

y conserva las características con cargas distintas de cero.

### Random Forest

Random Forest se utiliza como estrategia complementaria de clasificación y priorización de variables.

El R Markdown importa:

```text
randomforests_sigfeatures.csv
```

ordena las características por:

```text
MeanDecreaseAccuracy
```

y conserva el **5% superior**.

---

## Entradas del R Markdown

El script presupone que las etapas previas de MS-DIAL y MetaboAnalyst ya fueron completadas.

### GC-MS: alineamiento

Se espera el archivo:

```text
GC-MS/
└── Alignment_results/
    └── Area_peak.txt
```

El código detecta automáticamente la fila que comienza con:

```text
Alignment ID
```

y utiliza esa fila como encabezado de la tabla.

Durante la limpieza:

- se conserva `Metabolite name`;
- se seleccionan las columnas `CT`, `PD` y `QC`;
- se elimina el prefijo `no MS2:`;
- se eliminan metabolitos cuyo nombre comienza por `Unknown`;
- los nombres duplicados se hacen únicos.

---

### LC-MS positivo

Se esperan:

```text
LC-MS/
└── Positive/
    └── Results/
        ├── volcano.csv
        ├── splsda_loadings.csv
        ├── randomforests_sigfeatures.csv
        └── data_normalized.csv
```

---

### LC-MS negativo

Se esperan:

```text
LC-MS/
└── Negative/
    └── Results/
        ├── volcano.csv
        ├── splsda_loadings.csv
        ├── randomforests_sigfeatures.csv
        └── data_normalized2.csv
```

---

### GC-MS

Después del preprocesamiento se esperan:

```text
GC-MS/
└── Results/
    ├── volcano.csv
    ├── splsda_loadings.csv
    ├── randomforests_sigfeatures.csv
    └── data_normalized2.csv
```

---

## Estructura esperada

El código utiliza rutas relativas, por lo que debe conservarse una organización compatible con:

```text
Metabolomics/
├── README.md
├── Metabolomics.Rmd
│
├── LC-MS/
│   ├── Positive/
│   │   └── Results/
│   │       ├── volcano.csv
│   │       ├── splsda_loadings.csv
│   │       ├── randomforests_sigfeatures.csv
│   │       └── data_normalized.csv
│   │
│   └── Negative/
│       └── Results/
│           ├── volcano.csv
│           ├── splsda_loadings.csv
│           ├── randomforests_sigfeatures.csv
│           └── data_normalized2.csv
│
└── GC-MS/
    ├── Alignment_results/
    │   └── Area_peak.txt
    │
    └── Results/
        ├── volcano.csv
        ├── splsda_loadings.csv
        ├── randomforests_sigfeatures.csv
        └── data_normalized2.csv
```

Los archivos de datos y resultados indicados representan la **estructura esperada de entrada** y no necesitan formar parte del repositorio GitHub.

---

## Selección por consenso

Para cada plataforma se construyen tres conjuntos:

```r
ANOVA
sPLS.DA
RF
```

La intersección se calcula mediante uniones internas sucesivas:

```r
ANOVA %>%
  inner_join(sPLS.DA, by = "Metabolite") %>%
  inner_join(RF, by = "Metabolite") %>%
  inner_join(norm_data, by = "Metabolite")
```

De esta forma, únicamente permanecen las características presentes simultáneamente en las tres estrategias de selección y disponibles en la matriz normalizada.

Para cada plataforma se genera además un **diagrama de Venn** que documenta el número de características identificadas por cada estrategia y su intersección.

---

## Asignación PD / Healthy

Después de obtener los metabolitos consenso, el código utiliza conjuntamente la dirección del `log2FC` y la carga del primer componente de sPLS-DA para asignar las señales a PD o Healthy.

La orientación de las cargas depende del modelo generado para cada plataforma, por lo que la lógica utilizada se encuentra definida explícitamente dentro del R Markdown.

El resultado se incorpora en la columna:

```text
Group
```

---

## Consolidación entre plataformas

Los resultados consenso de las tres plataformas se integran mediante:

```r
bind_rows(
  LC.MS_positive,
  LC.MS_negative,
  GC.MS
)
```

y se agrega una variable `Type` que identifica el origen de cada metabolito:

```text
Positive
Negative
GC-MS
```

Posteriormente se generan tres archivos principales:

```text
Final_metabolites.csv
Healthy_metabolites.csv
PD_metabolites.csv
```

`Final_metabolites.csv` constituye la salida metabolómica consolidada utilizada posteriormente en la integración multiómica.

---

## Visualización de intersecciones

### Diagramas de Venn

Para cada plataforma se genera un diagrama de Venn con:

- ANOVA;
- sPLS-DA;
- Random Forest.

Las figuras se exportan como:

```text
LC-MS/Positive/Venn.pdf
LC-MS/Negative/Venn.pdf
GC-MS/Venn.pdf
```

### UpSet plot

Después de consolidar las tres plataformas, se genera una matriz binaria metabolito × plataforma.

El UpSet plot permite visualizar:

- metabolitos exclusivos de LC-MS(+);
- metabolitos exclusivos de LC-MS(-);
- metabolitos exclusivos de GC-MS;
- metabolitos compartidos entre plataformas.

La figura se exporta como:

```text
UpSet_metabolites.pdf
```

---

## Salidas principales

El workflow genera como principales productos:

- `GC-MS_filtered.csv`: matriz GC-MS después del filtrado por RSD;
- `LC-MS/Positive/Merged_Metabolites.csv`: metabolitos consenso de LC-MS(+);
- `LC-MS/Negative/Merged_Metabolites.csv`: metabolitos consenso de LC-MS(-);
- `GC-MS/Merged_Metabolites.csv`: metabolitos consenso de GC-MS;
- `Final_metabolites.csv`: firma metabolómica consolidada;
- `Healthy_metabolites.csv`: metabolitos asignados al grupo control;
- `PD_metabolites.csv`: metabolitos asignados a Parkinson;
- diagramas de Venn por plataforma;
- `UpSet_metabolites.pdf`: intersección entre plataformas.

---

## Resultados obtenidos

Después del control de calidad y preprocesamiento se analizaron **3,046 características metabólicas**:

| **Plataforma** | **Características analizadas** |
| --- | ---: |
| LC-MS positivo | 1,798 |
| LC-MS negativo | 285 |
| GC-MS | 963 |
| **Total** | **3,046** |

La estrategia de consenso entre ANOVA, sPLS-DA y Random Forest permitió seleccionar:

| **Plataforma** | **Metabolitos consenso** |
| --- | ---: |
| LC-MS positivo | 9 |
| LC-MS negativo | 14 |
| GC-MS | 12 |
| **Firma consolidada** | **35** |

De los **35 metabolitos** incluidos en la firma final:

- **20** estuvieron enriquecidos en pacientes con Parkinson;
- **15** estuvieron enriquecidos en controles sanos.

Entre las señales asociadas con PD se identificaron metabolitos pertenecientes principalmente a clases lipídicas y polioles, incluyendo galactinol, manitol, triacilglicéridos con enlaces éter, plasmalógenos y ésteres de colesterol.

En los controles se observaron, entre otros:

- ácido indol-3-acético;
- N-isovalerylglycine;
- 3β-TMS-5α,6α-epoxycolestane.

La firma consolidada constituye la salida de este módulo para la posterior integración con las demás capas ómicas.

---

## Ejecución del R Markdown

Desde RStudio, el documento puede ejecutarse mediante **Knit**.

Alternativamente, desde R:

```r
rmarkdown::render("Metabolomics.Rmd")
```

o desde una terminal con R instalado:

```bash
Rscript -e 'rmarkdown::render("Metabolomics.Rmd")'
```

El directorio de trabajo debe mantener las rutas relativas esperadas para las carpetas `LC-MS/` y `GC-MS/`.

Al finalizar la ejecución se registra:

```r
Sys.time()
getwd()
sessionInfo()
```

para facilitar la trazabilidad del entorno utilizado.

---

## Requerimientos computacionales

A diferencia de los pipelines implementados mediante Nextflow y ejecutados en HPC, este módulo no cuenta actualmente con un registro estandarizado de CPU, RAM pico, almacenamiento y tiempo total equivalente a los demás módulos.

Por esta razón, no se reportan estimaciones computacionales no verificadas.

El consumo dependerá principalmente de:

- tamaño de las matrices metabolómicas;
- número de muestras;
- número de características;
- etapas ejecutadas dentro de MS-DIAL;
- análisis realizados mediante MetaboAnalyst;
- generación de figuras y consolidación en R.

El R Markdown encargado de la etapa de consenso trabaja sobre matrices y resultados previamente procesados, por lo que representa únicamente una fracción del costo computacional total del workflow metabolómico.

---

## Reproducibilidad

La reproducción completa del análisis requiere distinguir dos componentes:

### Procesamiento externo

Debe reproducirse utilizando:

```text
MS-DIAL v5
MetaboAnalyst v6.0
```

manteniendo los parámetros de procesamiento, filtrado, normalización y análisis estadístico establecidos para el estudio.

### Posprocesamiento automatizado

La selección por consenso, consolidación entre plataformas y generación de resultados finales se reproduce mediante:

```text
Metabolomics.Rmd
```

El archivo contiene el código utilizado para:

- selección de características;
- integración de resultados;
- generación de diagramas;
- clasificación PD/HC;
- construcción de la firma metabolómica final.

---

## Documentación detallada

El flujo completo, incluyendo procesamiento mediante MS-DIAL, preprocesamiento en MetaboAnalyst, criterios estadísticos y consolidación final, puede incorporarse en:

[`Flujo completo de metabolómica`](https://github.com/linafperez/TG263006/blob/main/Metabolomics/Anexo_Metabolomica.pdf)

---

## Relación con otros módulos

La salida principal del módulo:

```text
Final_metabolites.csv
```

es utilizada posteriormente en:

[`../Integration`](https://github.com/linafperez/TG263006/blob/main/Integration)

para relacionar las señales metabolómicas con los resultados obtenidos mediante:

- metagenómica;
- metatranscriptómica;
- 16S rRNA.

La integración permite construir asociaciones funcionales entre **microorganismos, genes y metabolitos** dentro de la arquitectura multiómica desarrollada para el estudio del microbioma intestinal asociado con enfermedad de Parkinson.
