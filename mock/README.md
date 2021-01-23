La idea es testear el efecto de cuantizar los valores de calidad en el resultado del polishing de un ensamblado. Los datos son de una secuenciación de una comunidad mock de microbios, obtenidos con nanopore R10.

El paper en el que me guío es:

Nicholls, S. M., Quick, J. C., Tang, S., & Loman, N. J. (2019). Ultra-deep, long-read nanopore sequencing of mock microbial community standards. Gigascience, 8(5), giz043. <https://academic.oup.com/gigascience/article/8/5/giz043/5486468>

La idea es seguir los pasos del paper y ver que todo es reproducible, luego ensamblar con flye los últimos datos que publicaron (R10.3) y hacer el polishing siguiendo dinstintos protocolos que usan la calidad en algún punto.

Como algunos de estos programas usan heuristicas (al menos flye seguro), voy a repetir el mismo pipeline sobre cada versión de los datos 3 veces.

Los programas de polishing a usar son racon, medaka, marginPolish, y HELEN. De ellos, racon y marginPolish usan información de calidad por base. medaka y helen usan redes neurales recurrentes y fueron diseñados para corregir la salida de racon y marginPolish respectivamente. Notar que los que usan calidad no están al final del pipeline, sino que luego son corregidos por una red neural que no usa datos de calidad.

Del github de marginPolish: "If quality scores are present, the base likelihood is factored into alignment weight estimation".

Del paper de racon:  "Chunks of reads with an average quality lower than a predefined threshold are removed from the corresponding window. With this quality filter, we are able to use only high-quality parts of reads in the final consensus. Each window is then processed independently by constructing a POA graph and calling the consensus of the window. Quality values are used again during POA graph construction, where each edge is weighted by the sum of qualities of its source and destination nodes (bases; logarithmic domain)".

 <br/>

**Instalación de programas**

Instalé las versiones correspondientes al paper con conda:

```
conda create --prefix ./env_paper python=3.6 pilon=1.23 racon=1.3.2 minimap2=2.14 kraken2 checkm-genome quast=5.0.2
```
`wtdbg2` no estaba disponible en la versión 2.2, por lo que bajé los ejecutables de https://github.com/ruanjue/wtdbg2/releases/download/v2.2/wtdbg-2.2_x64_linux.tar.gz

`medaka` tampoco estaba disponible en la versión 0.5.0, por lo que lo instalé desde la fuente en un environment aparte.
Para eso descargue los archivos de https://pypi.org/project/medaka/0.5.0/#files. La versión correcta de medaka es importante porque la nueva no tiene el modelo correspondiente al basecaller que usaron en el paper, pero esta versión tiene problemas con dependencias que ya no están disponibles.

Antes de hacer make install a medaka hay que cambiar algunas dependencias del archivo `requirements.txt`:
- numpy==1.19.2 en lugar de numpy, y pasar esa linea arriba de la de tensorflow
- sacar la version a tensorflow (la que pide, 1.x.x, ya no está disponible como libreria de python)

 <br/>

Instalé las últimas versiones de los programas en otro environment:

```
conda create --prefix ./env_nuevos flye=2.8.2 racon=1.4.13 medaka=1.2.1 minimap2=2.17
```

MarginPolish/HELEN no está disponible vía `conda` pero sí está disponible vía `pip` bajo el paquete `helen`.

<br/>

**Datos y descarga**

Zymo Community Standards 2 (Even) Batch ZRC190633:
10 species (5 Gram-positive, 3 Gram-negative, 2 yeast): the bacteria are present at 12% and yeast at 2% (by genomic DNA)

Zymo Community Standards 2 (Log/Staggered) Batch ZRC190842
10 species (5 Gram-positive, 3 Gram-negative, 2 yeast) ranging from 10^2 - 10^8 genomic DNA abundance (total input 5 x 10^8 cells)

| Community | # of Spots    | # of Bases | Size  |
| ----------|:-------------:| ----------:|------:|
| Even      | 3,491,390     | 14.4G      | 12Gb  | 
| Log       | 3,667,480		| 16.5G      | 13.7Gb|


Esos son los datos del poro 9.41 secuenciados con la plataforma GridION (para el pipeline del paper). Descargué los fastQ de <https://github.com/LomanLab/mockcommunity>.

En el paper también usan los datos secuenciados con PromethION. Estos son muy pesados porque tienen mucha profundidad. En principio no los voy a usar.

También están los datos de señal, en principio no los descargué, para seguir los pasos del paper estoy trabajando con los datos tal cual los de ellos.

Descargué los datos de poro R10.3 de <https://nanopore.s3.climb.ac.uk/mock/Zymo-GridION-EVEN-3Peaks-R103-merged.fq.gz>. Son pocos reads (1.16M). También tienen publicados datos con el poro R10, son muchos mas reads (6.74M), disponibles en <https://nanopore.s3.climb.ac.uk/mock/Zymo-GridION-Even-3Peaks-Native-R10_hac_meth.fq.gz>.

Descargué las referencias (assemblies de spades con reads illumina) de <http://nanopore.s3.climb.ac.uk/mockcommunity/v2/Zymo-Isolates-SPAdes-Illumina.fasta>.

Descargué del SRA los reads de la secuenciación shotgun illumina de la comunidad even para el polishing con pilon: ERR2984773.

<br/>

**Algunos problemas con los programas**

Racon 1.3.2 carga muy lento los reads comprimidos, es mejor dárselos sin comprimir. 

minimap2 y racon no se pueden usar con nohup porque la salida estándar va a parar al archivo final. Se puede poner el comando dentro de un script y llamar al script con nohup.

El primer intento de correr racon 1.3.2 me dio el siguiente error:
```
[racon::Polisher::initialize] error: empty overlap set!
```
Hay que borrar el newline al final del archivo .paf de minimap2, y voila! Esto no me pasó con racon 1.4.3, pero con esa versión estuve usando `.sam` en lugar de `.paf`.

La versión 0.5.0 de medaka tiene problemas de incompatibilidad con la nueva versión de tensorflow. 
Tuve que reemplazar `tf` por `tf.compat.v1` en varios lugares de algunos scripts (todo lo que iba dando error), esos archivos fueron:

- medaka-0.5.0/venv/lib/python3.8/site-packages/keras/backend/tensorflow_backend.py
- medaka-0.5.0/venv/lib/python3.8/site-packages/medaka-0.5.0-py3.8-linux-x86_64.egg/medaka/inference.py

En el archivo medaka-0.5.0/venv/lib/python3.8/site-packages/medaka-0.5.0-py3.8-linux-x86_64.egg/medaka/common.py tuve que cambiar, en la linea 443, `raise StopIteration` por `return`. Esto es porque en las últimas versión de python (>3.7) cambió cómo se usan los generators: `StopIteration` ahora genera una excepción en lugar de simplemente terminar la vida útil del generator como antes (antes la excepción era silenciada).

La versión de `samtools view` que usa medakav1.2.1 (samtools 1.10) cuando llama a mini_align tiene un bug en el multithreading, por lo que edité la linea que llama a `samtools view` en el script mini_align (borré `-@ ${THREADS}`). 

También la parte de `medaka consensus` tiene un bug en el multithreading, si se indica 1 thread usa 16 threads, si se indica más de 1 thread toma todos los threads de la máquina. Edité el script medaka_consensus en la parte que llama a `medaka consensus`: en la linea 128 cambié `--threads ${THREADS}` por `--threads 1`, de esa forma puedo usar multithreading en otras partes.

HELEN falla antes de terminar a veces y no tengo claro por qué (no da ningún error, directamente termina). No es un tema de memoria. Estoy usando la versión de docker que no me da problemas.

<br/>

**Selección de contigs más largos**

Antes de hacer el polishing seleccioné los 8 contigs más largos:

```
infoseq --only  --length --name asm-flye-TODOS.fa | awk '{print $2, $1}' | sort -nr | head -n8 | cut -d' ' -f2 > lista_contigs_mayores
fastaUtils.pl -u asm-flye-TODOS.fa | grep --no-group-separator -A 1 -w -f lista_contigs_mayores > r10-flye-selected.fasta
```

<br/>

**Evaluación**

Para la evaluación estoy usando metaQUAST. TODO:  probar [fastmer.py](https://github.com/jts/assembly_accuracy) (per-genome accuracy).
La idea es usar las referencias de illumina para evaluar los missmatches antes y después de cada etapa polishing de los 8 contigs mas largos, tal como hicieron en el paper. Los missassemblies los evaluaron usando las referencias de pacbio, pero eso no es relevante para este trabajo.

<br/>

**Cuantizar los datos de calidad**

Hice un script en python que toma como entrada un fastq y devuelve el fastq con sólo 4 valores de calidad:

- si qscore <= 7 -> newqscore = 6
- si 8 <= qscore <= 13 -> newqscore = 12
- si 14 <= qscore <= 19 -> newqscore = 18
- si qscore >= 20 -> newqscore = 24

Es una primera idea de cómo definir los bins, no hay nada que respalde esta elección, pero me parece claro que los bins de illumina no funcionan para los datos de nanopore porque la calidad media de illumina es mucho mayor.

Hice el binning de los valores de calidad de los reads de la secuencicion con el poro r10, y a esos datos le apliqué el mismo pipeline que a los datos originales.
