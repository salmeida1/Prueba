# Análisis de datos metagenómicos de regiones ITS de hongos filamentosos asociados al suelo de la Hacienda "El Prado" - Sangolquí 

**Autores:** 
**Ing. Almeida Stefany**, **Ing. Córdova Daniel** y **Ing. Drouet Ariana**

**Fecha:**  marzo 2024

# 1.  Problema
#### Los hongos, estimados en 1,5 millones de especies y ubicados en casi todos los ecosistemas terrestres, desempeñan roles ecológicos clave, como la descomposición de materia orgánica y la simbiosis con plantas. En Ecuador, el estudio y la caracterización de la diversidad fúngica son especialmente relevantes debido a su ubicación geográfica diversa y sus ecosistemas únicos. Sin embargo, su estudio se ve obstaculizado por su naturaleza efímera y la dificultad para identificar especies a partir de sus estructuras visibles. En particular, en la Hacienda "El Prado"-IASA I, no se ha caracterizado la microbiota fúngica de los suelos, y debido a que los cuerpos fructíferos similares a menudo representan varias especies distintas, se subraya la necesidad de abordajes moleculares centrados en el estudio de la región espaciadora transcrita interna (ITS), para así, obtener datos en cuanto a su taxonomía a nivel de género y subgénero.

# 2.  Antecedentes
#### La ausencia de estudios previos sobre la microbiota fúngica en los suelos del sector de Horticultura y Fruticultura de la Hacienda “El Prado” - IASA I, ubicada en la provincia Pichincha, cantón Rumiñahui, parroquia Sangolquí, destaca la necesidad de realizar esta investigación. Un proyecto que implica la evaluación pionera de la diversidad de hongos filamentosos, principalmente en aquellos con potencial benéfico para la biorrecuperación de los suelos. 

# 3.  Objetivos

## 3.1  Objetivo general
#### Analizar los datos metagenómicos de las regiones ITS amplificadas de hongos filamentosos utilizando herramientas bioinformáticas.

## 3.2  Objetivos específicos

#### * Realizar el control de calidad de los datos crudos de secuencia de las regiones ITS por medio de la herramienta FASTQc en terminal de Linux.
#### * Identificar molecularmente a las especies de hongos filamentosos mediante un BLASTN del NCBI.
#### * Extraer la información taxonómica con la alineación de las secuencias consenso mediante MSA en T-coffee y analizarla en la plataforma Galaxy Europe. 

# 4.  Flujograma de trabajo

![Flujograma](https://github.com/Irondaniel34/Proyecto_G1/blob/main/Capturas_de_pantalla/Trabajo_final_flujograma1.jpg?raw=true)

# 5.  Análisis de datos

## 5.1  Importación y preprocesamiento de los datos 
#### Los datos crudos provienen de muestras de hongos filamentosos asociados al suelo de la Hacienda "El Prado", los cuales fueron aislados en laboratorio para posteriormente ser identificados molecularmente mediante la secuenciación por el método Sanger y el análisis de la región ITS amplificada con los primers ITS1 (TCCGTAGGTGAACCTGCGG) e ITS4 (TCCTCCGCTTATTGATATGC) reportados por White et al. (1990).
#### Se realizó la importación de la data procediente del secuenciador en formato *.ab1*. Por lo cual, previo al análisis de calidad es necesario realizar la conversión al formato *.fastq*, la misma que puede realizarse mediante el uso de BioPhyton (a través de Google Colab), usando la plataforma Galaxy Europe o mediante el uso de la terminal del sistema operativo Linux, de acuerdo con la siguiente metodología:

### 5.1.1  Conversión de secuencias mediante BioPhyton
#### La conversión de secuencias en BioPhyton se realizó en dos fases; la primera que consistió en la carga de los archivos *.ab1* dentro de una carpeta dada por el usuario, de acuerdo con el siguiente código:
```
# Creación de la carpeta de entrada dentro del directorio "content"
import os
carpeta_entrada = input("Ingrese el nombre de la nueva carpeta: ")
carpeta_entrada_path = os.path.join("/content", carpeta_entrada)
os.makedirs(carpeta_entrada_path)

# Impresión de la ruta al nuevo directorio
print(f"Directorio", carpeta_entrada, "creado en:", carpeta_entrada_path)

# Subida de archivos.ab1 al directorio de trabajo actual
from google.colab import files
ITS_files = files.upload()

# Mover archivos a la carpeta de entrada solicitada pro el usuario
import shutil
for filename in ITS_files:
    src = os.path.join("/content", filename)
    dst = os.path.join(carpeta_entrada_path, filename)
    shutil.move(src, dst)
```
#### La carga de archivos .ab1 se muestra a continuación:

![Imagen1_galaxy](https://github.com/Irondaniel34/Proyecto_G1/blob/dc9492430763dcff0a7347c0b3520ad188615f7b/Capturas_de_pantalla/Biophyton_carga%20de%20archivos%20ab1.jpg)


#### La segunda fase consistió en la conversión de las secuencias a formato *.fastq* para lo cual se utilizó el siguiente código:
```
# Instalación de Biophyton e importación de paquetes
!pip install biopython
from Bio import SeqIO
import os

# Ingreso de la ruta a la carpeta de entrada dada por el usuario.
input_folder = input("Ingrese el nombre de la carpeta que contiene las secuencias.ab1: ")
input_folder_path = os.path.join("/content", input_folder)

# Ingreso de la ruta a la carpeta de salida dada por el usuario.
output_folder = input("Ingrese el nombre de la carpeta donde se guardarán los archivos .fastq: ")
output_folder_path = os.path.join("/content", output_folder)

# Creación de la carpeta de salida en caso que esta no existiera en el directorio
os.makedirs(output_folder_path, exist_ok=True)


# Conversión de todos los archivos .ab1 de la carpeta de entrada por iteración con la función for
for filename in os.listdir(input_folder_path):
    # Compruebación si el archivo es un archivo .ab1
    if filename.endswith(".ab1"):
        # Determinación de la ruta del archivo .ab1
        input_file_path = os.path.join(input_folder_path, filename)

        # Determinación de la ruta de almancenamiento del archivo .fastq
        output_file_path = os.path.join(output_folder_path, filename.split(".")[0] + ".fastq")

        # Converción del archivo de .ab1 a .fastq usando Biopython
        SeqIO.convert(input_file_path, "abi", output_file_path, "fastq")

# Impresión del mensaje de confirmación
print(f"Archivos convertidos y guardados en el directorio: {output_folder_path}")

# Descarga del archivo comprimido con las secuencias en formato .fastq
!zip -r {output_folder}.zip {output_folder_path}
files.download(f"{output_folder}.zip")
```
#### La conversión realizada de los archivos *.ab1* a *.fastq* mientras el código se encuentra corriendo en Biophyton se muestra a continuación:

![Imagen1_galaxy](https://github.com/Irondaniel34/Proyecto_G1/blob/dc9492430763dcff0a7347c0b3520ad188615f7b/Capturas_de_pantalla/Biophyton_conversi%C3%B3n%20de%20archivos%20ab1%20a%20fastq.jpg)

#### El cuaderno de Google Colab se puede encontrar en el siguiente documento: [Proyecto_final_G1_ab1_to_fastq.ipynb](https://github.com/Irondaniel34/Proyecto_G1/blob/a27b121c371726325dae74519bc4dcacb4aa3c32/Proyecto_final_G1_ab1_to_fastq.ipynb)

### 5.1.2  Conversión de secuencias mediante Galaxy Europe

#### Se cargaron en modo de colección los datos en formato ab1 tanto para ITS1 e ITS4.

![Imagen1_galaxy](https://github.com/Irondaniel34/Proyecto_G1/blob/main/Capturas_de_pantalla/Galaxy1.png)

#### Se construyeron las secuencias 

![Imagen1_galaxy](https://github.com/Irondaniel34/Proyecto_G1/blob/main/Capturas_de_pantalla/Galaxy2.2.png)

#### Y posteriormente se transofrmaron de formato *.ab1* a *.fastq* 

![Imagen1_galaxy](https://github.com/Irondaniel34/Proyecto_G1/blob/main/Capturas_de_pantalla/Galaxy3.png)

### 5.1.3  Conversión de secuencias mediante terminal de Linux

#### Se descargan los datos crudos del secuenciador en la máquina virtual, estos vienen en un archivo comprimido *.zip*, para lo cual se usa el comando: 
```
unzip nombre_del_archivo.zip
```
![Imagen1_linux](https://github.com/Irondaniel34/Proyecto_G1/blob/main/Capturas_de_pantalla/Linux1.jpg)

#### Para continuar con el procesamiento de las secuencias se tiene que convertir de un archivo *.ab1* a *.fastq* para lo cual se utilizaron los siguientes comandos:
```
sudo apt install emboss 
```
#### Para instalar el paquete EMBOSS, que contienen la herramienta SEQRET que nos permite el paso de .ab1 a .fastq
![Imagen2_linux](https://github.com/Irondaniel34/Proyecto_G1/blob/main/Capturas_de_pantalla/instalaci%C3%B3n%20emboss.jpg)

**Nota:** *El paquete EMBOSS fue previamente instalado, es por eso que no sale el mensaje tradicional de su instalación*
```
mkdir ITS_fastq1 # Para crear la carpeta de salida de los archivos FASTQ
for file in ITS1/*.ab1;do  # Para realizar un bucle, en el que los archivos dentro de la carpeta ITS1 que tengan la extensión AB1, realicen la conversión
seqret -sequence $file -outseq ITS_fastq1/$(basename $file .ab1).fastq -osformat2 fastq; done # Comando seqret para la conversión de .ab1 a .fastq y su salida a la carpeta previamente creada. 
```
![Imagen3_linux](https://github.com/Irondaniel34/Proyecto_G1/blob/main/Capturas_de_pantalla/ab12fastq.jpg)

#### Esto fue realizado en las dos carpetas que contenian archivos *.ab1*, el directorio ITS1 e ITS4.

## 5.2  Control de calidad de las secuencias

### 5.2.1 Informe de calidad en FASTQc mediante terminal de Linux

#### Se efectuó el análisis de calidad de las secuencias mediante la herramienta FASTQc, que se ejecutó mediante línea de comando de Linux. 
#### Para ello, se tuvo que primero instalar el paquete *FASTQc* en la terminal utilizando el siguiente comando:
```
sudo apt install fastqc
```
![Imagen4_linux](https://github.com/Irondaniel34/Proyecto_G1/blob/main/Capturas_de_pantalla/instalacion%20fastqc.jpg?raw=true)
**Nota:** El paquete FASTQc ya había sido previamente instalado.

#### Posteriormente se realizó el análisis FASTQc, esto crea archivos con extensión *.html*, en la que nos permitirá visualizar las calidades, esto con los comandos:
```
mkdir ITS_fastqc1 # Para crear un directorio al cual irán los archivos .html
fastqc ITS_fastq1/*.fastq -o ITS_fastqc1 # Todos los archivos dentro de la carpeta ITS_fastq1 realizarán el control de calidad FASTQc y el argumento -o envía los resultados a la carpeta que se creó previamente.
```
![Imagen5_linux](https://github.com/Irondaniel34/Proyecto_G1/blob/main/Capturas_de_pantalla/Fastqc%20y%20salida%20a%20otra%20carpeta.jpg?raw=true)

![Imagen6_linux](https://github.com/Irondaniel34/Proyecto_G1/blob/main/Capturas_de_pantalla/Fastqc%20en%20carpeta.jpg?raw=true)
Archivos *.html* en carpeta creada.
![Imagen7_linux](https://github.com/Irondaniel34/Proyecto_G1/blob/main/Capturas_de_pantalla/fasqcWEB.jpg?raw=true)
Visualización de FASTQc mediante archivo *.html*

### 5.2.1    Recorte de secuencias en Trimmomatic mediante terminal de Linux

#### También se efectuó la herramienta *Trimmomatic* por medio de línea de comando de Linux, se tomó esta decisión al ver que mediante otras plataformas existía problemas al cargar los archivos.
#### Para la realización de los cortes de *Trimmomatic*, se consultó a la inteligencia artificial, la cuál indicó:

Descargar la aplicación de Trimmomatic de la página oficial.

**Nota:** Página web de la herramienta Trimmomatic
![Imagen8Linux](https://github.com/Irondaniel34/Proyecto_G1/blob/main/Capturas_de_pantalla/P%C3%A1gina%20de%20descarga%20de%20trimo.jpg?raw=true)
[USADELLAB.org](http://www.usadellab.org/cms/index.php?page=trimmomatic)

**Link de descarga**

Se ejecutó el comando de descarga y tras lo cual se efectuó un unzip del archivo descargado usando:
```
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
unzip Trimmomatic-0.39.zip
```
![Imagen9Linux](https://github.com/Irondaniel34/Proyecto_G1/blob/main/Capturas_de_pantalla/Descarga%20y%20unzip%20de%20trimo.jpg?raw=true)

La ejecución de la herramienta se llevó a cabo con línea de comando, usando:

```
mkdir F1T #CREA CARPETA DESTINO
for file in ITS_fastq1/*.fastq; do    
filename=$(basename -- "$file")
output="F1T/trimmed_${filename}"
java -jar Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads 4 "$file" "$output" \
ILLUMINACLIP:Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:10 LEADING:20 TRAILING:20 \
MINLEN:36 -phred33
done
```

#### Este comando creó un directorio (F1T) al cual salieron las nuevas secuencias cortadas con la herramienta, se especificó para el corte de tipo *Leading* y *Trailing* y se especifica que tipo de criterio de calidad de Phred utilizar, en este caso 33.
![Imagen10Linux](https://github.com/Irondaniel34/Proyecto_G1/blob/main/Capturas_de_pantalla/Trimme%20linux.jpg?raw=true)
![Imagen11Linux](https://github.com/Irondaniel34/Proyecto_G1/blob/main/Capturas_de_pantalla/ls%20trimo.jpg?raw=true)

###### Nota: A las secuencias resultantes se les volvió a aplicar un control de calidad con FASTQc y se descartaron 7 secuencias que no cumplieron  con los estándares de los autores, las cuales se muestran a continuación.

![Imagen1Fastqc](https://github.com/Irondaniel34/Proyecto_G1/blob/a34be7ad09072a06a515c51272a31b035ca81801/Capturas_de_pantalla/Fastqc_Secuencias%20baja%20calidad%20descartadas.jpg)

#### Las secuencias finales se utilizaron en la plataforma Galaxy Europe para la realización de la secuencia consenso.


## 5.3    Obtención de secuencias consenso en la Plataforma Galaxy Europe

### 5.3.1    Preparación de la colección de secuencias pareadas ITS1 e ITS4

#### Se importaron los datos de las secuencias ITS1 e ITS4 en el área de *Collection* seleccionando en los parámetros: de tipo de colección: *List of Pairs* y en tipo de archivo *fastqsanger*:
![Imagen1Colección](https://github.com/Irondaniel34/Proyecto_G1/blob/ad01d8cd195e91146fe2198a9e6b1327b32e73b0/Capturas_de_pantalla/Coleccion_1.jpg)

#### Se configuraron los parámetros para que las archivos con el prefijo trimmed_ITS1- y trimmed_ITS4- sean reconocidos como lecturas forward y reverse, respectivamente: 
![Imagen2Colección](https://github.com/Irondaniel34/Proyecto_G1/blob/ad01d8cd195e91146fe2198a9e6b1327b32e73b0/Capturas_de_pantalla/Coleccion_2.jpg)

#### Se editó los atributos de la colección ITS para el tipo de archivo *fastqsanger*
![Imagen3Colección](https://github.com/Irondaniel34/Proyecto_G1/blob/ad01d8cd195e91146fe2198a9e6b1327b32e73b0/Capturas_de_pantalla/Coleccion_3.jpg)

#### Colección de archivos en formato *.fastq* pareados en Galaxy 

![Imagen4Colección](https://github.com/Irondaniel34/Proyecto_G1/blob/adeb62d070629b957acd133be9966fa576264b04/Capturas_de_pantalla/Coleccion_4.jpg)


### 5.3.2    Obtención de secuencias consenso ITS

#### Se ejecutó la herramienta *PEAR* utilizando la colección pareada ITS como entrada con el objetivo de obtener las secuencias ensambladas.
![ImagenConsenso](https://github.com/Irondaniel34/Proyecto_G1/blob/32d209c86518699ad4fbfa3ccc51290bfb70a33b/Capturas_de_pantalla/Consenso.jpg)


### 5.3.3    Conversión de secuencias consenso ITS a formato *.fasta*

#### Con el objetivo de identificar molecularmente las secuencias utilizando la herramienta BLASTn del NCBI se ejecutó la herramienta *FASTQ to FASTA* utilizando la colección se secuencias consenso *Pear on collection: Assembled reads*.
![Imagen1FASTQtoFASTA](https://github.com/Irondaniel34/Proyecto_G1/blob/d61a5c7200c5181f54400531fe90bc10fa86adb2/Capturas_de_pantalla/Fastq_to_FASTA_1.jpg)

#### Las secuencias obtenidas en formato *.fasta* se muestran a continuación:
![Imagen1FASTQtoFASTA](https://github.com/Irondaniel34/Proyecto_G1/blob/d61a5c7200c5181f54400531fe90bc10fa86adb2/Capturas_de_pantalla/Fastq_to_FASTA_2.jpg)


## 5.4  Identificación molecular de las secuencias en BLASTn del NCBI 

### 5.4.1   Obtención archivo MultiFasta

#### Para el *blasteo* de multiples especies en el NCBI se necesitaba realizar un archivo *.fasta* concatenado con todas las secuencias. Motivo por el que se usó la terminal de Linux para unificar los numerosos archivos *.fasta*. 

#### Esto se realizó utilizando los siguientes comandos:

```
for file in Fasta/*.fasta; do # Realiza en bucle de todos los archivos con extensión .fasta
cat "$file" >> ITS_SC.fasta # Con cat permite obtener un concatenado que se nombrará ITS_SC.fasta
done
```

#### Obtención de FASTA concatenado con las secuencias consenso.

![Imagen12Linux](https://github.com/Irondaniel34/Proyecto_G1/blob/main/Capturas_de_pantalla/cat%20FASTA.jpg?raw=true)
![Imagen13Linux](https://github.com/Irondaniel34/Proyecto_G1/blob/main/Capturas_de_pantalla/concatenado%20FASTA.jpg?raw=true)


### 5.4.2   BLAST en NCBI 

#### El archivo consolidado de secuencias fue cargado en BLASTn del NCBI, el cuál entregó los resultados como se aprecia:

![Imagen14W](https://github.com/Irondaniel34/Proyecto_G1/blob/main/Capturas_de_pantalla/blasteo.jpg?raw=true)
![Imagen15W](https://github.com/Irondaniel34/Proyecto_G1/blob/main/Capturas_de_pantalla/blast2.jpg?raw=true)
![Imagen10Linux](https://github.com/Irondaniel34/Proyecto_G1/blob/main/Capturas_de_pantalla/Galaxy10.png)

#### Se hizo el BLAST de múltiples secuencias, para múltiples especies que se fueron clasificando en la siguiente tabla:

|Identificación | Género/Especie | # de Acceso |E value|   %|
|-------------- |--------------  |-------------|---------|---------|
|ITS1-10|	*Fusarium oxysporum*	|MT276045.1 |0.0	|100|
|ITS1-11|	*Clonostachys rosea*	|MN452687.1|0.0	|99.34|
|ITS1-12|	*Trichoderma atroviride*	|MT514373.1|0.0	|99.80|
|ITS1-13|	*Fusarium oxysporum*	|MT276045.1|0.0	|99.77|
|ITS1-14|	*Mucor moelleri*|	MN270302.1|0.0	|100|
|ITS1-1|	*Mucor circinelloides*|	JF723599.1|0.0|	98.70|
|ITS1-24|	*Trichoderma atroviride*|	MT514373.1|0.0|	99.80|
|ITS1-28|	*Trichoderma atroviride*|	MT514373.1|0.0|	100|
|ITS1-29|	*Fusarium ramigenum*|	MH980135.1|0.0|	100|
|ITS1-2|	*Mucor circinelloides*|	JF723599.1|0.0	|99.43|
|ITS1-32|	*Mucor hiemalis*	|MT366055.1|0.0|	99.28|
|ITS1-34|	*Fusarium ramigenum*|	MH980135.1|0.0|	100|
|ITS1-35|	*Trichoderma hamatum*|	PP464109.1|0.0|100|
|ITS1-36|	*Fusarium circinatum*|	MT464451.1|0.0|	99.77|
|ITS1-37|	*Fusarium oxysporum*|	PP453660.1|0.0	|100|
|ITS1-38|	*Trichoderma atroviride*|	MT514373.1|0.0	|100|
|ITS1-39|	*Clonostachys solani*|	AF358243.1|0.0	|100|
|ITS1-42|	*Xylariales* sp.|	KT269510.1|0.0|	98.36|
|ITS1-44|	*Mucor janssenii*|	MH855051.1|0.0	|99.64|
|ITS1-4	|*Clonostachys divergens*	|OQ910570.1|0.0	|100|
|ITS1-5|	*Fusarium equiseti*	|PP464187.1|0.0	|99.77|
|ITS1-7|	*Chaetomium cochliodes*|	MT561402.1|0.0|99.15|
|ITS1-8|	*Minimedusa polyspora*|	MH859968.1|0.0|99.30|
|ITS1-9|	*Trichoderma hamatum*	|MT271927.1|0.0|100|
|ITS1-A|	*Clonostachys rosea*|	OQ513907.1|0.0	|100|
|ITS1-BB|	*Trichoderma asperellum*|	LC123601.1|0.0|	99.79|
|ITS1-CC|	*Mucor* sp.|	ON209714.1|0.0|96.87|
|ITS1-C	|*Mucor* sp.	|MK164211.1|0.0	|99.81|
|ITS1-D|	*Fusarium culmorum*	|MT640271.1|0.0	|100|
|ITS1-E	|*Fusarium equiseti*	|MN135744.1	|0.0	|100|
|ITS1-F|	*Mortierella* sp.	|MF423582.1	|0.0	|100|
|ITS1-I	|*Trichoderma asperelloides*	|PP336494.1	|0.0	|100|
|ITS1-J	|*Fusarium equiseti*|	PP464187.1 |	0.0	|100|
|ITS1-K|	*Trichoderma* sp.	|OP999633.1	|0.0	|99.80|
|ITS1-L	|*Clonostachys rosea*	|MN452687.1	|0.0	|99.34|
|ITS1-M	|*Clonostachys* sp.	|OM106429.1|	0.0	|99.57|
|ITS1-N|	*Clonostachys solani*|	OQ910839.1	|0.0	|100|
|ITS1-NN|	*Trichoderma koningiopsis*|	OM574767.1	|0.0|	100|
|ITS1-O	|*Trichoderma* sp.	|PP464108.1	|0.0	|99.80|
|ITS1-Q	| *Trichoderma* sp.|	MT740343.1|	0.0	|100|
|ITS1-R10	|*Absidia* sp.	|KU923829.1 	|0.0	|98.35|
|ITS1-R12	|*Fusarium culmorum*	|MT640271.1 	|0.0	|100|
|ITS1-R13|	*Monocillium* sp.|	 LC433827.1| 	0.0|	99.8|
|ITS1-R14|	*Mucor* sp.	|MT762712.1	|0.0	|100|
|ITS1-R16|	*Fusarium* sp.|	MT276045.1| 	0.0	|99.77|
|ITS1-R18|	*Fusarium graminearum*|	OQ979826.1|	0.0	|100|
|ITS1-R19|	*Trichoderma* sp.|	MN522765.1 |	0.0|	99.6|
|ITS1-R1|	*Epicoccum nigrum*|	PP469495.1|	0.0|	100|
|ITS1-R20|	*Fusarium graminearum*|	OQ979826.1| 	0.0|	100|
|ITS1-R21|	*Fusarium culmorum*|	MT640271.1|	0.0|	100|
|ITS1-R22|	*Trichoderma* sp.|	MN522765.1| 	0.0|	100|
|ITS1-R24|	*Trichoderma atroviride*	|MZ568323.1	|0.0	|100.00|
|ITS1-R25|	*Trichoderma* sp.	|MK870863.1|	0.0|	99.8|
|ITS1-R26|	*Trichoderma atrobrunneum*|	PP385058.1|	0.0	|99.81|
|ITS1-R27|	*Mucor moelleri*	|MG214587.1	|0.0	|99.17|
|ITS1-R28|	*Poculum* sp.	|KF725734.1	|0.0	|99.14|
|ITS1-R29|	*Clonostachys rhizophaga*	|OQ910695.1|	0.0	|99.78|
|ITS1-R2|	*Fusarium oxysporum*	|MT276045.1	|0.0	|100.00|
|ITS1-R31|	*Fusarium equiseti*|	PP464187.1	|0.0	|99.55|
|ITS1-R33	|*Fusarium oxysporum*|	PP437857.1|	0.0|	100.00|
|ITS1-R34|	*Mucor hiemalis*	|LC413619.1	|0.0|	99.29|
|ITS1-R35|	*Mucor circinelloides*|	OP597934.1|	0.0|	100.00|
|ITS1-R36|	*Penicillium camemberti*|	MT529802.1	|	0.0|	100.00|
|ITS1-R3	|*Irpex laceratus*	|LC431580.1	|0.0	|97.81|
|ITS1-R4	|*Pseudopithomyces chartarum*|	MN077447.1|	0.0	|99.41|
|ITS1-R5	|*Trichoderma spirale*	|KX858890.1	|0.0	|99.61|
|ITS1-R6	|*Lecanicillium cf. psalliotae*	|AB517935.1	|0.0	|99.71|
|ITS1-R7|	*Clonostachys rosea*	|MK304207.1	|0.0	|100.00|
|ITS1-R8|	*Fusarium graminearum*|	ON024849.1	|0.0	|100.00|
|ITS1-R	|*Trichoderma* sp.	|MT514373.1	|0.0	|100.00|
|ITS1-T	|*Trichoderma asperellum*	|OL868970.1	|0.0	|100.00|
|ITS1-U	|*Trichoderma koningiopsis*|	OM574767.1	|0.0	|100.00|
|ITS1-W	|*Trichoderma asperellum*	|PP464108.1	|0.0	|100.00|
|ITS1-X|	*Trichoderma koningiopsis*	|KC428393.1	|0.0	|99.40|

 


## 5.6 Alineamiento Multiple en T-Coffee EMBL-EBI
Con el fin de analizar la biodiversidad y la filogenia de los hongos, se optó por realizar un alineamiento múltiple en la herramienta T-Coffee EMBL-EBI de las secuencias ITS.

![Imagen14W](https://github.com/Irondaniel34/Proyecto_G1/blob/main/Capturas_de_pantalla/Galaxy11.png)
![Imagen14W](https://github.com/Irondaniel34/Proyecto_G1/blob/main/Capturas_de_pantalla/Galaxy12.png)

De esta manera se obtuvo en formato FASTA las secuencias alineadas, mismo que se descargó para continuar con el proceso de análisis en la plataforma de Galaxy Europe
![Imagen16W](https://github.com/Irondaniel34/Proyecto_G1/blob/main/Capturas_de_pantalla/Galaxy13.png)

## 5.7 Clasificación taxonómica en Galaxy Europe 
Se cargó el archivo en formato FASTA obtenido del alineamiento múltiple en la plataforma de T-Coffee. 
![Imagen16W](https://github.com/Irondaniel34/Proyecto_G1/blob/main/Capturas_de_pantalla/Galaxy14.png)

Después con ayuda de la herramienta Kraken2, se obtuvo la clasificación taxonómica 
![Imagen16W](https://github.com/Irondaniel34/Proyecto_G1/blob/main/Capturas_de_pantalla/Galaxy15.png)
![Imagen16W](https://github.com/Irondaniel34/Proyecto_G1/blob/main/Capturas_de_pantalla/Galaxy16.png)
![Imagen16W](https://github.com/Irondaniel34/Proyecto_G1/blob/main/Capturas_de_pantalla/Galaxy17.png)
![Imagen16W](https://github.com/Irondaniel34/Proyecto_G1/blob/main/Capturas_de_pantalla/Galaxy18.png)
![Imagen16W](https://github.com/Irondaniel34/Proyecto_G1/blob/main/Capturas_de_pantalla/Galaxy19.png)
![Imagen16W](https://github.com/Irondaniel34/Proyecto_G1/blob/main/Capturas_de_pantalla/Galaxy20.png)
![Imagen16W](https://github.com/Irondaniel34/Proyecto_G1/blob/main/Capturas_de_pantalla/Galaxy21.png)
![Imagen16W](https://github.com/Irondaniel34/Proyecto_G1/blob/main/Capturas_de_pantalla/Galaxy22.png)

Posteriormente se realiza una transformación de los formatos de salida con la herramienta Krakentools: Convert kraken report file

![Imagen16W](https://github.com/Irondaniel34/Proyecto_G1/blob/main/Capturas_de_pantalla/Galaxy23.png)

## 5.8 Visualización 
Dado que los resultados obtenidos con la herramienta kraken2 se transformaron en un archivo reporte, se cargó la herramienta Krona para así visualizar los datos taxonómicos de forma interactiva.

![Imagen16W](https://github.com/Irondaniel34/Proyecto_G1/blob/main/Capturas_de_pantalla/Galaxy24.png)

![Imagen16W](https://github.com/Irondaniel34/Proyecto_G1/blob/main/Capturas_de_pantalla/Galaxy25.png)

## 6.Conclusiones

#### - Se realizó el control de calidad por medio de la herramienta FASTQc en terminal de Linux de 84 secuencias, de las cuales 77 (91.7%) secuencias cumplieron los parámetros de calidad de los autores para seguir con el proceso de análisis de calidad.
#### - Se identificó molecularmente las diferentes especies de hongos filamentosos con un porcentaje de identidad mayor al 99%, con un E value de 0.0, entre los cuales los géneros más abundantes fueron *Trichoderma*, *Fusarium* y *Mucor* fueron los más abundantes de entre la población.
#### - Se extrajo la información taxonómica de las secuencias mediante el uso de plataformas como Kraken2 y la herramienta Krona Pie Chart del Galaxy Europe.

## 7. Recomendaciones
#### - El uso de la terminal de Linux proporciona un control sobre el sistema operativo que facilita realizar tareas complejas de manera eficiente y rápida, a comparación del uso de la plataforma Galaxy Europe.
#### - Descartar las secuencias que presenten problemas durante la corrida del FastQC en el control de calidad, para así, garantizar la calidad óptima de los datos analizados en la plataforma Galaxy, permitiendo que avance la corrida de datos de manera adecuada en posteriores análisis.


## 8. Referencias

#### Batut et al., 2018 Community-Driven Data Analysis Training for Biology Cell Systems 10.1016/j.cels.2018.05.012
#### EMBL-EBI, Wellcome Genome Campus, Hinxton, Cambridgeshire, CB10 1SD, UK.
#### Hiltemann, Saskia, Rasche, Helena et al., 2023 Galaxy Training: A Powerful Framework for Teaching! PLOS Computational Biology 10.1371/journal.pcbi.1010752
#### Nilsson, R. H., Ryberg, M., Abarenkov, K., Sjökvist, E., & Kristiansson, E. (2009). The ITS region as a target for characterization of fungal communities using emerging sequencing technologies. FEMS Microbiology Letters, 296(1), 97-101. https://doi.org/10.1111/j.1574-6968.2009.01618.x
#### Sophia Hampe, Bérénice Batut, Paul Zierep, Taxonomic Profiling and Visualization of Metagenomic Data (Galaxy Training Materials). https://training.galaxyproject.org/training-material/topics/microbiome/tutorials/taxonomic-profiling/tutorial.html Online; accessed Sun Mar 24 2024
#### The Galaxy Community. The Galaxy platform for accessible, reproducible and collaborative biomedical analyses: 2022 update, Nucleic Acids Research, Volume 50, Issue W1, 5 July 2022, Pages W345–W351, doi:10.1093/nar/gkac247

