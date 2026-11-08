# Integración multiómica

Módulo de integración multiómica desarrollado para relacionar la composición microbiana, la actividad génica y los perfiles metabolómicos asociados con la enfermedad de Parkinson (PD) y controles sanos (HC).

## Esquema general

[`Diagrama general de integración multiómica.pdf`](./diagrama_integracion.pdf)

El análisis integra resultados provenientes de los módulos de **metagenómica, 16S rRNA, metatranscriptómica y metabolómica**, con el fin de identificar asociaciones entre microorganismos, genes, funciones biológicas y metabolitos.

## Software y requerimientos computacionales

Los siguientes valores corresponden al consumo observado durante la ejecución del módulo de integración multiómica.

| **Software / entorno** | **Módulo / etapa** | **CPU** | **RAM asignada (GB)** | **RAM pico (GB)** | **Almacenamiento** | **Tiempo** |
|---|---|---:|---:|---:|---:|---:|
| R / R Markdown | Integración multiómica | 8 | 8 | 2.80 | 550 MB | 00:10:51 |

La ejecución completa del módulo requirió aproximadamente **10 min 51 s**, con un consumo máximo de memoria de **2.80 GB**.

## Entradas

El módulo utiliza como entradas principales:

- **MAGs** obtenidos mediante el análisis metagenómico.
- **ASVs y asignaciones taxonómicas** provenientes del análisis 16S rRNA.
- **genes diferencialmente expresados asociados a MAGs** obtenidos mediante metatranscriptómica.
- **metabolitos diferenciales** identificados mediante el análisis metabolómico.

## Integración taxonómica

La integración entre metagenómica y 16S rRNA se realiza mediante la comparación de la clasificación taxonómica de **MAGs y ASVs**.

La concordancia se establece principalmente a nivel de **género**, permitiendo relacionar los microorganismos identificados mediante ambas aproximaciones y seleccionar aquellos con evidencia consistente entre los conjuntos de datos.

## Integración funcional

Los microorganismos identificados se relacionan posteriormente con:

- genes diferencialmente expresados.
- funciones y rutas biológicas.
- metabolitos diferenciales.

Esto permite construir asociaciones funcionales del tipo:

```text
Microorganismo → Gen → Función/Ruta → Metabolito
```

## Salidas principales

El módulo genera:

- correspondencias taxonómicas entre MAGs y ASVs.
- asociaciones entre microorganismos y genes diferencialmente expresados.
- relaciones funcionales entre genes, rutas y metabolitos.
- matrices de correlación para PD y HC.
- enriquecimiento KEGG y GO.
- diagramas de integración multiómica.

## Resultados obtenidos

En el grupo PD se identificaron relaciones que involucraron microorganismos como **Parabacteroides goldsteinii** y diferentes clases de lípidos, mientras que en HC se observaron asociaciones relacionadas con **Roseburia intestinalis** y metabolitos como esteroles.

Estos resultados permitieron caracterizar cambios en las relaciones **microbio–gen–metabolito** entre pacientes con Parkinson y controles sanos.

## Ejecución

El análisis de integración se desarrolla en **R / R Markdown**, utilizando como entrada los resultados generados por los diferentes módulos ómicos.

Los scripts deben ejecutarse una vez se encuentren disponibles los resultados finales de metagenómica, 16S rRNA, metatranscriptómica y metabolómica.

## Documentación detallada

El flujo completo de integración, incluyendo las relaciones entre las distintas capas ómicas, criterios de asociación y etapas internas del análisis, se encuentra en:

[`Flujo completo de integración multiómica.pdf`](./Anexo_Integracion.pdf)

## Relación con otros módulos

Este módulo integra los resultados generados en:

- [`../Metagenomics`](https://github.com/linafperez/TG263006/blob/main/Metagenomics)
- [`../16S`](https://github.com/linafperez/TG263006/blob/main/16S)
- [`../Metatranscriptomics`](https://github.com/linafperez/TG263006/blob/main/Metatranscriptomics)
- [`../Metabolomics`](https://github.com/linafperez/TG263006/blob/main/Metabolomics)

y constituye la etapa final de integración de la arquitectura computacional multiómica.
