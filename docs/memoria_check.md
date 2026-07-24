

Las enfermedades oncológicas, inflamatorias y neurológicas son un
problema clínico de primera magnitud. Los sistemas de liberación de
fármacos (DDS, del inglés ) ofrecen una solución concreta: encapsular el
principio activo, protegerlo de la degradación y liberarlo donde hace
falta. Las nanopartículas poliméricas basadas en poli(ácido
láctico--glicólico) (PLGA), poli(ácido láctico) (PLA) o policaprolactona
llevan décadas en desarrollo clínico, y los dendrímeros de
poliamidoamina (PAMAM), que son moléculas muy ramificadas con cavidades
internas para cargar fármacos, representan una alternativa
complementaria \[@danhier2012; @kesharwani2018\].

Existen un cierto conjunto de proteínas que determinan qué le pasa a
cualquier molécula que entra en el organismo. La P-glicoproteína
(P-gp/MDR1) expulsa los fármacos que reconoce como sustratos, lo que
causa resistencia a la quimioterapia \[@waghray2018\]. El CYP3A4 los
metaboliza en el hígado antes de que lleguen al tumor, y es responsable
del metabolismo de más del 50% de los fármacos clínicos \[@zanger2013\].
El receptor de transferrina (TfR1) y el receptor de folato (FR$\alpha$)
están sobreexpresados en las células tumorales y se pueden usar como
puerta de entrada para llevar nanopartículas específicamente al tumor.
La lisozima es un indicador de compatibilidad tisular. Y la albúmina
sérica (HSA) transporta los fármacos en sangre y determina cuánto tiempo
circulan.

Evaluar experimentalmente cómo interacciona cada componente de un DDS
con estas seis proteínas es inviable a gran escala. El docking
computacional, que simula cómo una molécula pequeña se une a una
proteína, permite cribar cientos de moléculas antes de tocar un tubo de
ensayo \[@trott2010\]. Los estudios existentes tienen limitaciones
claras: se centran en una proteína, trabajan con conjuntos pequeños de
moléculas, y ninguno integra cálculos mecanocuánticos de estructura
electrónica para obtener geometrías precisas y propiedades como punto de
partida \[@maingi2012\].

Este trabajo construye una plataforma computacional para evaluar 58
moléculas representativas de nanopartículas poliméricas y dendrímeros
frente a las seis proteínas descritas anteriormente. Los objetivos
concretos son: obtener geometrías moleculares precisas y propiedades
electrónicas mediante cálculos DFT y perfiles de polaridad superficial
(COSMO-RS); realizar docking sistemático y analizar las interacciones
estructurales con PLIP 2025; calcular descriptores moleculares y
propiedades ADMET; entrenar modelos de aprendizaje automático para
predecir la afinidad de unión y clasificar moléculas para drug delivery;
e implementar una herramienta web que recomiende el carrier más adecuado
para una molécula nueva.

## Base de datos molecular

Se construyó una biblioteca de 60 moléculas en cuatro categorías:
monómeros y productos de degradación de polímeros biodegradables (ácido
láctico, glicólico, $\varepsilon$-caprolactona y derivados), oligómeros
y unidades de dendrímero PAMAM (poliaminas, oligómeros de quitosano),
fármacos modelo (doxorrubicina, paclitaxel, 5-fluorouracilo,
metotrexato, curcumina, cisplatino, gemcitabina, sirolimus,
dexametasona, ibuprofeno) y ligandos de targeting activo (ácido fólico,
biotina, galactosa, manosa, PEG-aminas y péptido RGD). Los SMILES,
representación textual de la estructura química de una molécula, se
obtuvieron de la base de datos PubChem \[@pubchem\]. Las coordenadas 3D
iniciales se generaron con el algoritmo ETKDGv3 y se preoptimizaron
mediante simulaciones de mecánica molecular de campo de fuerzas con
MMFF94 \[@rdkit\].

## Optimización de geometrías (DFT)

Los cálculos de estructura electrónica se realizaron con TURBOMOLE 7.8.1
\[@turbomole\] en el clúster HPC del Centro Extremeño de Investigación,
Innovación Tecnológica y Supercomputación (CENITS) (16 núcleos, 31,GB
RAM por nodo). El cálculo se llevó a cabo con el funcional híbrido B3LYP
con corrección de dispersión D4 \[@grimme2010\] y base de funciones
def2-TZVP \[@weigend2005\], al nivel de la función de densidad
electrónica (DFT), con aproximación RI-J para reducir el coste
computacional. Para cada molécula se realizaron tres cálculos: (i)
optimización de la geometría hasta encontrar la estructura de mínima
energía, (ii) cálculo de la distribución de carga superficial (COSMO-RS)
para obtener el perfil sigma, y (iii) cálculo de frecuencias de
vibración IR para confirmar que la estructura corresponde a un mínimo
real. El cisplatino requirió un tratamiento especial con el
pseudopotencial ECP-60-mwb de Stuttgart para el platino, empleando la
base def2-SVP para este átomo y def2-TZVP para el resto; el condroitín
sulfato, con carga $-2$, se recalculó especificando la carga formal.

## Acoplamiento molecular

Las seis proteínas de barrera (Tabla ) se obtuvieron del Protein Data
Bank (PDB), base de datos de acceso libre que almacena estructuras
tridimensionales de macromoléculas biológicas \[@pdb\]. Las proteínas se
prepararon con MGLTools, añadiendo hidrógenos y asignando cargas
atómicas Gasteiger para modelar de forma eficiente el potencial
electrostático de la macromolécula mediante la estimación de sus cargas
parciales.. Los ligandos se convirtieron a formato PDBQT, formato de
entrada requerido por AutoDock Vina que codifica las coordenadas
atómicas junto con las cargas parciales y los grados de libertad
rotacionales del ligando, con el objetivo de definir la flexibilidad
conformacional durante el docking, con OpenBabel. El docking, que simula
cómo una molécula pequeña se une a una proteína, se realizó con AutoDock
Vina 1.2.5 \[@eberhardt2021\] con parámetros de alta exhaustividad (32),
9 poses por cálculo y 8 CPUs. Las cajas de búsqueda se centraron en los
sitios de unión conocidos de cada proteína. De 348 dockings posibles se
obtuvieron 342 válidos.

## Análisis de interacciones proteína-ligando

Las interacciones no covalentes entre los cinco candidatos con mayor
Índice de Eficiencia de Transporte (ITI), métrica desarrollada en este
trabajo para cuantificar el perfil global de una molécula para drug
delivery, y las tres proteínas más relevantes (P-gp, TfR1 y FR$\alpha$)
se analizaron con PLIP 2025 (del inglés ) \[@plip2025\], una herramienta
que detecta automáticamente qué tipo de contactos establece cada ligando
con los residuos de la proteína: puentes de hidrógeno, contactos
hidrofóbicos, puentes salinos e interacciones entre anillos aromáticos.
Se procesaron 15 complejos en total.

## Descriptores moleculares

Se calcularon inicialmente 4179 descriptores numéricos brutos para
caracterizar cada molécula, agrupados en cinco familias: descriptores 2D
calculados con RDKit \[@rdkit\] (38 valores: peso molecular, lipofilia,
superficie polar, grupos funcionales, índices de complejidad
topológica); perfiles sigma COSMO-RS (100 valores que describen la
distribución de densidad electrónica superficial); fingerprints
Morgan/ECFP4 (2048 bits que codifican los patrones de subestructura
química del entorno de cada átomo hasta un radio definido) y MACCS keys
(167 bits que representan la presencia o ausencia de subestructuras
químicas predefinidas) con RDKit; y descriptores 3D calculados con
Mordred \[@moriwaki2018\] (1826 valores sobre las geometrías DFT, entre
los que destacan los momentos de inercia, el volumen molecular, los
índices de forma y los autovalores de la matriz de distancias). Tras
eliminar columnas con más del 50,% de valores ausentes y columnas de
varianza nula (bits de fingerprint constantes para todas las moléculas
del dataset, o descriptores 3D no calculables para alguna estructura),
el espacio final de descriptores efectivamente utilizado en los modelos
se reduce a 1097 variables: 100 perfiles sigma, 529 bits de
Morgan/ECFP4, 120 bits de MACCS, 38 descriptores 2D y 310 descriptores
3D de Mordred (incluyendo 4 que Mordred recalcula con el mismo nombre
que un descriptor 2D — BertzCT, TPSA, LabuteASA y MW — y que se
conservan como variables independientes).

## Propiedades ADMET

Las propiedades de ADMET (Absorción, Distribución, Metabolismo,
Excreción y Toxicidad) se calcularon con SwissADME \[@swissadme\] para
las 60 moléculas. Se obtuvieron 49 parámetros por compuesto, destacando
cinco modelos de lipofilia (LogP), solubilidad en agua, absorción
intestinal, permeabilidad hematoencefálica, condición de sustrato de
P-gp, inhibición de enzimas CYP450, *score* de biodisponibilidad oral,
alertas de toxicidad, el cumplimiento de las reglas de Lipinski (donde
incurrir en más de una violación de sus límites de peso molecular,
lipofilia y enlaces de hidrógeno predice una deficiente absorción por
vía oral) y accesibilidad sintética.

## Modelos de aprendizaje automático

El dataset final tiene 56 moléculas (de las 60 de la biblioteca
original, tras exigir que cada molécula tenga geometría DFT convergida,
docking válido frente a las seis proteínas y descriptores completos) y
1097 descriptores tras la limpieza descrita arriba. Antes del
entrenamiento se aplica escalado estadístico (RobustScaler) y se
seleccionan las 200 variables más informativas para cada proteína
(SelectKBest). Se desarrollaron tres modelos:

Regresión de afinidad: Random Forest (300 árboles), Gradient Boosting
(200 estimadores) y SVR (máquina de vectores soporte) para predecir la
energía de unión $\Delta G$ en cada proteína. Se evalúa con validación
cruzada de 5 particiones usando el Coeficiente de Determinación ($R^2$,
proporción de la varianza explicada por el modelo), el Error Absoluto
Medio (MAE, promedio de los errores en magnitud absoluta) y la Raíz del
Error Cuadrático Medio (RMSE, que mide la desviación cuadrática promedio
y penaliza en mayor medida los errores grandes).

Índice de Eficiencia de Transporte (ITI): puntuación de 0 a 100 que
combina las seis afinidades de docking en una métrica única. Premia la
afinidad por TfR1, FR$\alpha$ (40,% del índice), HSA (25,%) y lisozima
(20,%). Penaliza la afinidad por P-gp y CYP3A4 (15,%). Incluye un factor
corrector por peso molecular () que penaliza moléculas muy pequeñas o
muy grandes, con umbrales de 0,1 para MW,$<$,60,Da, 0,4 para 60–100,Da,
0,7 para 100–150,Da, 1,0 para 150–900,Da y 0,8 para MW,$>$,900,Da.

Clasificador binario: Random Forest, Gradient Boosting y SVC para
clasificar moléculas como favorables (ITI $\geq$ 50) o desfavorables.
Evaluado mediante el Área Bajo la Curva ROC (AUC-ROC, que mide la
capacidad del modelo para distinguir entre clases), la precisión
(proporción de predicciones positivas que son correctas) y la puntuación
F1 (media armónica entre precisión y exhaustividad que evalúa el balance
global del clasificador).

La interpretabilidad de los modelos se analizó con valores SHAP
\[@lundberg2017\], que miden cuánto contribuye cada descriptor a la
predicción final.

Se implementó una aplicación Streamlit que, a partir del SMILES de una
molécula, calcula descriptores, predice el perfil de afinidad frente a
las seis proteínas, calcula el ITI y genera un diagrama radar del perfil
de transporte con una recomendación del carrier más adecuado.

## Cálculos DFT y perfiles COSMO-RS

Las 58 estructuras moleculares convergieron sin frecuencias imaginarias,
lo que confirma que todas son mínimos de energía reales. Los cálculos
más exigentes fueron Paclitaxel (113 átomos pesados) y Sirolimus (144
átomos), que requirieron varios cientos de ciclos de optimización. Los
perfiles sigma COSMO-RS describen la distribución de densidad
electrónica en la superficie de cada molécula: valores de $\sigma$
alejados de cero indican regiones polares (negativas si hay exceso de
carga, positivas si hay defecto), mientras que la concentración cerca de
$\sigma \approx 0$ refleja carácter apolar. Los resultados mostraron
patrones claros por grupo: la doxorrubicina (fármaco modelo) presenta
una distribución ancha que se extiende por toda la región polar,
mientras que el ácido láctico (monómero de degradación) concentra su
densidad cerca de $\sigma \approx 0$ (figura ).

## Acoplamiento molecular

La figura muestra la matriz de afinidad completa, en la que cada celda
representa la energía de unión $\Delta G$ (kcal/mol) de la pose más
favorable (modo 1) obtenida por AutoDock Vina para un par
molécula-proteína determinado: valores más negativos indican mayor
afinidad de unión. Cabe señalar que, aunque se emplea la notación
$\Delta G$, los valores reportados corresponden a energías de *scoring*
de docking y no a energías de Gibbs termodinámicas rigurosas; su
interpretación es por tanto comparativa entre moléculas, no absoluta.
Las energías medias por proteína se calcularon promediando los valores
de modo 1 de las 57 moléculas frente a cada proteína, obteniendo:
FR$\alpha$ ($\overline{\Delta G} = -5{,}70$ kcal/mol), CYP3A4
($-5{,}21$), TfR1 ($-5{,}08$), P-gp ($-5{,}07$), HSA ($-4{,}31$) y
Lisozima ($-3{,}54$). Los valores individuales más negativos, que
indican mayor afinidad, fueron Paclitaxel en P-gp ($-10{,}27$ kcal/mol)
y Methotrexate tanto en CYP3A4 ($-9{,}91$) como en FR$\alpha$
($-9{,}91$); el ácido fólico, pese a ser el ligando natural de
FR$\alpha$, obtuvo una afinidad algo menor ($-9{,}67$) que el
metotrexato, un antifolato estructuralmente relacionado. Sirolimus en
CYP3A4 obtuvo $-0{,}15$ kcal/mol, la interacción más débil del dataset,
lo que sugiere que esta molécula elude el metabolismo hepático.

El mapa de correlaciones (figura ), que permite identificar si las
moléculas que se unen bien a una proteína tienden también a unirse bien
a las demás, muestra que P-gp, TfR1, FR$\alpha$ y HSA están muy
correlacionadas entre sí (correlación de Pearson $r$ entre 0,90 y 0,95).
La correlación más baja es entre CYP3A4 y Lisozima ($r = 0{,}73$), lo
que tiene sentido dado que sus mecanismos de reconocimiento son
distintos. La figura de dispersión $\Delta G$ frente a peso molecular
(figura ), que permite conocer en qué medida el tamaño molecular
determina la afinidad de unión, confirma la tendencia general: las
moléculas más grandes se unen con mayor afinidad a todas las proteínas
debido a la existencia de más sitios de interacción. Los fármacos modelo
aparecen sistemáticamente en la zona de alta afinidad.

## Análisis de interacciones proteína-ligando (PLIP)

El análisis de los 15 complejos generó 708 contactos no covalentes en
total (figura ). FR$\alpha$, el receptor de folato que se sobreexpresa
en tumores de ovario y mama, concentra la mayor parte: 515 contactos, de
los que 340 son puentes de hidrógeno. Los residuos SER174 (67 contactos)
y ARG103 (64) son los más implicados; ambos son aminoácidos hidrofílicos
con grupos laterales capaces de formar puentes de hidrógeno, lo que
explica su papel central en el reconocimiento del ácido fólico según
describen los datos cristalográficos. Los puentes salinos con ASP81 y
LYS136 en los oligómeros de quitosano, aminoácidos de carga opuesta (ASP
negativo, LYS positivo) que establecen interacciones electrostáticas
fuertes, sugieren que Chitotriose y Chitobiose se unen en la misma
cavidad de unión que el folato natural.

En P-gp, que es la principal bomba de expulsión de fármacos del
organismo, 96 de los 100 contactos totales son hidrofóbicos, repartidos
entre ILE730, LEU312, LEU757 y LEU760 (8 cada uno). Este patrón confirma
lo que predice el análisis SHAP: P-gp reconoce sus sustratos por su
forma y complejidad molecular, no por grupos químicos concretos. En
TfR1, el receptor de transferrina que se usa para llevar moléculas al
interior de células tumorales, aparecen 82 puentes de hidrógeno, con
THR305 (10 contactos) y GLU972 (6) a la cabeza.

## Modelos de regresión de afinidad

Los modelos de regresión de afinidad tienen como objetivo predecir la
energía de unión $\Delta G$ de cualquier molécula frente a cada una de
las seis proteínas de barrera, a partir exclusivamente de sus
descriptores moleculares, sin necesidad de realizar un nuevo cálculo de
docking. Esto permite evaluar moléculas nuevas de forma inmediata una
vez entrenado el modelo.

Cada molécula se representa mediante un vector de 1097 descriptores que
concatena las cinco familias calculadas tras la limpieza descrita en
Materiales y Métodos: 38 descriptores 2D, 100 valores del perfil sigma
COSMO-RS, 529 bits de fingerprints Morgan/ECFP4, 120 bits de MACCS keys
y 310 descriptores 3D de Mordred, todos ellos calculados sobre las
geometrías optimizadas a nivel DFT. Antes del entrenamiento, SelectKBest
selecciona los 200 descriptores más correlacionados con el $\Delta G$ de
cada proteína concreta, de modo que el subconjunto de variables
seleccionadas es distinto para P-gp, TfR1, FR$\alpha$ y el resto. El
modelo entrena sobre esos 200 descriptores seleccionados mediante un
pipeline independiente por proteína.

SVR y Gradient Boosting fueron los modelos más precisos en la mayoría de
las proteínas, con Random Forest competitivo pero no siempre mejor
(figura ). El mejor modelo alcanzó R$^2 = 0{,}836$ en P-gp (SVR, MAE =
$0{,}371$ kcal/mol), $0{,}857$ en TfR1 (Gradient Boosting) y $0{,}769$
en Lisozima (SVR), todos por encima del umbral de $0{,}7$ que se
considera aceptable para modelos cuantitativos de actividad. En CYP3A4
($0{,}543$, Random Forest), FR$\alpha$ ($0{,}669$, SVR) y HSA
($0{,}775$, SVR) el rendimiento es más modesto o limítrofe. Con 56
muestras y validación cruzada de 5 particiones, cada fold de test tiene
11 moléculas, lo que limita la capacidad predictiva y explica la
desviación estándar relativamente alta de estos R$^2$. El modelo MLP
obtuvo sistemáticamente R$^2$ negativo en las seis proteínas, indicio de
sobreajuste severo dado el tamaño muestral y se descarta como candidato
viable para este dataset.

Los puntos con mayor error de predicción en la figura de predicción
vs. real (figura ) son Paclitaxel y Condroitín sulfato disacárido en
P-gp, que se utiliza para modelar interacciones con glicosaminoglicanos
en la superficie celular, probablemente porque su complejidad
estructural está fuera del rango habitual del dataset. Ammonia es el
principal outlier en TfR1 por ser la molécula más pequeña del conjunto.

## Índice de Eficiencia de Transporte (ITI)

El ITI es un índice con valor entre 0 y 100 que resume mediante un único
número el perfil completo de una molécula para drug delivery: penaliza
la afinidad por P-gp (expulsión) y CYP3A4 (metabolismo), y premia la
afinidad por TfR1 y FR$\alpha$ (entrada en células tumorales), HSA
(transporte en sangre) y lisozima (compatibilidad con tejidos).

El ranking (figura ) lo encabeza Triethylene glycol (PEG, ITI,=,100),
seguido de Tetraethylenepentamine ($87{,}3$), Pentaethylenehexamine
($85{,}0$), Chitotriose ($84{,}7$), Chitobiose ($82{,}7$), D-Galactose
($81{,}7$) y D-Mannose ($81{,}7$). Los fármacos modelo quedan al final
del dataset porque tienen alta afinidad por P-gp. El podio (Triethylene
glycol, Tetraethylenepentamine, Pentaethylenehexamine) es estable frente
a la reconstrucción del descrita en la nota metodológica de la sección
anterior; el orden exacto a partir del cuarto puesto es más sensible a
esa reconstrucción y debe confirmarse frente al script original si se
recupera. La posición alta de Chitotriose y Chitobiose, en cualquier
caso, está respaldada de forma independiente por el análisis estructural
PLIP (sección 3.3), que no depende del ITI.

La distribución por grupos (figura ) muestra qué ligandos de targeting y
oligómeros tienen medianas de los valores ITI alrededor de 70–72, los
monómeros de degradación en torno a 60, y los fármacos modelo cerca de
42 con alta variabilidad. Las figuras y , que representan el perfil de
afinidad de las ocho moléculas con mayor ITI en forma de diagrama radar,
ilustran el patrón de Triethylene glycol: radio pequeño en P-gp y
CYP3A4, moderado en TfR1 y FR$\alpha$.

## Clasificador y valores SHAP

Las curvas ROC con bandas de confianza se muestran en la figura . Un AUC
cercano a 1 significa que el modelo separa casi perfectamente moléculas
favorables de desfavorables. Random Forest obtuvo AUC-ROC = $0{,}957$,
Accuracy = $0{,}927$ y F1 = $0{,}947$, donde el AUC-ROC mide la
capacidad discriminativa del modelo independientemente del umbral de
clasificación, la Accuracy indica el porcentaje de moléculas
clasificadas correctamente, y el F1 es la media armónica entre precisión
y sensibilidad, útil cuando las clases no están perfectamente
balanceadas. Gradient Boosting alcanzó AUC = $0{,}904$ y Accuracy =
$0{,}891$. SVC y MLP rindieron notablemente peor (AUC $0{,}636$ y
$0{,}626$ respectivamente), probablemente por el desbalance de clases
(37 moléculas favorables frente a 19 desfavorables con el ITI corregido)
combinado con el alto número de variables respecto al tamaño muestral.

el factor corrector por peso molecular del ITI (, Materiales y Métodos)
se recuperó a partir de los valores observados en los resultados
guardados, ya que el script que lo calculaba no quedó versionado junto
al resto del código; se documenta explícitamente en para que el Modelo 2
sea reproducible en el futuro.

El análisis SHAP muestra qué descriptores moleculares son los más
importantes para predecir la afinidad por P-gp (figura ). BertzCT, que
mide la complejidad topológica de una molécula cuantificando el número
de caminos distintos en su grafo molecular — representación abstracta de
la estructura química en la que los átomos son nodos y los enlaces son
aristas —, presenta en el dataset valores que van desde 0 (Ammonia)
hasta 2280 (Paclitaxel), con los fármacos modelo concentrados entre 940
y 2280 y los monómeros de degradación entre 6 y 221. Este descriptor
tiene un peso SHAP de $0{,}27$, muy por encima de NumBonds ($0{,}11$),
que cuenta el número total de enlaces de la molécula; MW ($0{,}11$), que
es el peso molecular; y ExactMW ($0{,}07$), que es el peso molecular
exacto calculado a partir de las masas isotópicas. Los tres miden de
distintas formas el tamaño y la complejidad de la molécula. En la
práctica, esto significa que P-gp expulsa moléculas complejas y grandes,
sin importar mucho su química concreta – el mismo patrón cualitativo que
en la versión anterior del análisis, ahora con pesos SHAP recalculados
sobre el espacio de descriptores completo (incluyendo sigma-profiles y
Mordred 3D).

## Propiedades ADMET

El análisis ADMET con SwissADME (figura ) mide propiedades
farmacológicas clave: qué porcentaje de moléculas de cada grupo se
absorbe bien en el intestino (GI absorption), cuántas son expulsadas por
P-gp, y qué índice de biodisponibilidad oral tienen.

Los monómeros de degradación tienen 0 violaciones de la regla de
Lipinski (que establece que una molécula con peso molecular $\leq$
500,Da, LogP $\leq$ 5, no más de 5 dadores y 10 aceptores de puentes de
hidrógeno tiene probabilidad de ser un fármaco oral) y el 67% con
absorción intestinal alta. Los oligómeros presentan una media de 1,5
violaciones y 65% de absorción. Los ligandos de targeting obtienen el
mejor índice de biodisponibilidad oral (BA $\approx$ 0,85) y 87% de
absorción. Los fármacos modelo tienen en media 1,75 violaciones, 20% de
absorción intestinal alta y el 80% son sustratos de P-gp. Los monómeros
de degradación no son sustratos de P-gp en ningún caso.

## Herramienta web de recomendación de carrier

Como resultado final de la plataforma computacional, se implementó una
herramienta web interactiva desarrollada con Streamlit que permite
evaluar cualquier molécula nueva sin necesidad de ejecutar ningún script
ni cálculo externo. El flujo de trabajo es el siguiente: el usuario
introduce el SMILES de la molécula, la herramienta calcula
automáticamente los descriptores 2D con RDKit, aplica los modelos Random
Forest entrenados para predecir el $\Delta G$ frente a las seis
proteínas de barrera, calcula el ITI y genera un diagrama radar con el
perfil de afinidad completo. A partir de ese perfil, la herramienta
emite una recomendación del carrier más adecuado entre las opciones
evaluadas en este trabajo (nanopartículas PLGA/PLA, dendrímeros PAMAM o
recubrimientos de quitosano), indicando además si la molécula se
clasifica como favorable o desfavorable para drug delivery según el
clasificador binario.

La herramienta está diseñada para ser utilizada en fases tempranas del
diseño de nanocarriers, como filtro rápido previo a cálculos más
costosos como el docking o la simulación DFT.

## Geometrías DFT y perfiles sigma

Usar estructuras moleculares optimizadas a nivel mecanocuántico (DFT) en
lugar de geometrías estimadas por métodos de campo de fuerza, que no
consideran explícitamente los electrones de las moléculas, tiene mayor
grado de verosimilitud. Para moléculas grandes con tensión
conformacional, como el sirolimus o el paclitaxel, las diferencias
geométricas son suficientes para cambiar los perfiles sigma. Estos
perfiles describen cómo se distribuye la polaridad en la superficie de
la molécula, algo que los descriptores globales como LogP o TPSA no
capturan con la misma precisión. Es importante resaltar que los estudios
previos que calculan descriptores directamente desde SMILES asumen una
geometría que puede estar lejos de aquella que se corresponde con el
mínimo real de energía.

## Interpretación de los modelos

BertzCT domina el análisis SHAP con un peso de $0{,}27$, frente a
$0{,}11$ de NumBonds y $0{,}11$ de MW. Los tres miden, de distintas
formas, el tamaño y la complejidad de la molécula. P-gp tiene un sitio
de unión de gran tamaño, de unos $5000$,Å$^3$, y reconoce sus sustratos
principalmente por su hidrofobicidad y complejidad topológica, no por
grupos químicos específicos \[@seelig1998\]. Los resultados coinciden
con lo que describe la literatura.

La correlación alta entre P-gp y TfR1 ($r = 0{,}95$) tiene una
implicación práctica importante: diseñar una molécula con alta afinidad
por TfR1 (para dirigirla a células tumorales) probablemente también
aumentará su afinidad por P-gp (que la expulsará). En este contexto el
valor ITI se hace necesario, dado que una métrica que mida sólo la
afinidad por TfR1 no detecta ese conflicto.

El rendimiento moderado en CYP3A4 (mejor R$^2 = 0{,}543$) y FR$\alpha$
(mejor R$^2 = 0{,}669$) es esperable para una muestra de 56 moléculas
con validación cruzada de 5 particiones. Con 11 muestras de test por
fold, cabe esperar un margen limitado para cualquier modelo. SVR y
Gradient Boosting funcionan mejor que Random Forest en la mayoría de las
proteínas bajo estas condiciones, aunque con desviaciones estándar altas
entre folds.

## Validación estructural con PLIP

Los 708 contactos identificados por PLIP dan contexto molecular a los
resultados de docking. En FR$\alpha$, los residuos SER174 (67 contactos)
y ARG103 (64) son exactamente los que la cristalografía identifica como
esenciales para reconocer el ácido fólico. Que los oligómeros de
quitosano se acoplen en la misma cavidad sugiere que podrían funcionar
como ligandos de targeting: FR$\alpha$ está sobreexpresado en tumores
ováricos y de mama, así que una molécula que compita con el folato por
ese sitio podría usarse para dirigir un nanocarrier específicamente a
esos tumores.

En P-gp, los 96 contactos hidrofóbicos con ILE730, LEU312, LEU757 y
LEU760 (8 cada uno) confirman lo que predice SHAP. Triethylene glycol y
Tetraethylenepentamine forman muy pocos de esos contactos, lo que
explica cuantitativamente por qué tienen ITI alto.

## Relevancia del ITI y propiedades ADMET

Triethylene glycol (PEG, ITI,=,100) tiene el perfil más limpio: 0
violaciones de Lipinski, biodisponibilidad oral alta y 0% de sustrato
P-gp. Las ramas PAMAM de baja generación (Tetraethylenepentamine con
ITI,=,$87{,}3$ y Pentaethylenehexamine con $85{,}0$) combinan buena
afinidad por TfR1 con poca interacción con P-gp.

La información más relevante proviene del análisis ADMET: el 80% de los
fármacos modelo son sustratos de P-gp, frente al 0% de los monómeros de
degradación. Eso cuantifica el problema que resuelven los nanocarriers:
si el fármaco está encapsulado en una partícula cuyo recubrimiento no es
sustrato de P-gp, la bomba no puede expulsarlo directamente. Es una
ventaja que los ensayos clínicos con formulaciones de nanopartículas
PLGA y PAMAM han confirmado experimentalmente \[@soma2000;
@danhier2012\].

Chitotriose (ITI,=,$84{,}7$) y Chitobiose ($82{,}7$) tienen un perfil
equilibrado: buena afinidad por FR$\alpha$ confirmada por PLIP, poca
interacción con P-gp y pocas violaciones de Lipinski. El quitosano ya se
usa ampliamente como recubrimiento de nanopartículas en formulaciones
clínicas, así que estos resultados refuerzan las evidencias
experimentales.

## Limitaciones

El tamaño del dataset es la limitación principal. Con 56 moléculas y
cálculos DFT completos, con un alto grado de exactitud, el coste
computacional es el factor que lo acota. La ausencia de validación
experimental es inherente al carácter computacional del trabajo. Los
valores de $\Delta G$ obtenidos son coherentes con los publicados en
estudios de docking frente a P-gp y FR$\alpha$ con moléculas de
características similares, donde se reportan rangos típicos de $-4$ a
$-10$ kcal/mol \[@dolghih2011; @maingi2012\]. Los residuos identificados
por PLIP coinciden con datos cristalográficos. Los pesos del ITI se
definieron con criterio experto, es decir, asignados manualmente en
función del conocimiento farmacológico sobre el papel de cada proteína
en el transporte de fármacos, sin ajuste empírico sobre datos
experimentales reales. Con datos de biodistribución real (mediciones de
concentración del fármaco en tejidos obtenidas en ensayos in vivo) se
podrían calibrar esos pesos de forma objetiva y cuantitativa, en lugar
de depender del juicio del investigador.

Durante la preparación final de este trabajo se detectó y corrigió un
error en el script de combinación de descriptores (): una discrepancia
de formato en los nombres de molécula entre (que usa guiones bajos) y el
resto de archivos (que usan espacios) hacía que el por nombre solo
emparejara una fracción de las moléculas, lo que a su vez provocaba que
las familias de perfiles sigma COSMO-RS y de descriptores 3D de Mordred
se descartaran por completo del conjunto de entrenamiento durante la
limpieza de columnas con datos ausentes. Todas las cifras de este
capítulo corresponden a la versión corregida del pipeline (nombres
normalizados, sin descriptores duplicados, con las cinco familias de
descriptores realmente presentes). Los resultados de regresión se
mantienen en el mismo orden de magnitud que con el pipeline original e
incluso mejoran para varias proteínas, lo que sugiere que las
conclusiones del trabajo no dependían de ese error, pero sí subraya la
importancia de una comprobación de cobertura del como parte del pipeline
(incluida ahora en ).

El AUC-ROC obtenido por el clasificador ($0{,}957$ Random Forest,
$0{,}904$ Gradient Boosting) debe interpretarse con cautela dado el
tamaño muestral: con 56 moléculas y validación cruzada de 5 particiones,
cada fold de test contiene únicamente 11 muestras, lo que produce curvas
ROC de forma escalonada y valores de AUC más extremos que los que cabría
esperar con un conjunto de test mayor. Además, la selección de 200
descriptores de un total de 1097 candidatos, aplicada dentro de cada
fold de entrenamiento (evitando fuga de información), opera en un
régimen de alta dimensionalidad relativa a la muestra ($p \gg n$) que
favorece un cierto optimismo en las métricas de validación cruzada. Para
descartar que la separación observada se deba al azar, se realizó un
test de permutación: se barajaron aleatoriamente las etiquetas
favorable/desfavorable 20 veces y se repitió el mismo procedimiento de
validación cruzada sobre cada barajado. El AUC medio obtenido con
etiquetas aleatorias fue $0{,}52$ (máximo $0{,}73$ en las 20
repeticiones), muy por debajo del AUC real, lo que respalda que la
capacidad discriminativa del modelo refleja una señal genuina en los
datos y no un ajuste espurio al ruido. Aun así, la clasificación es una
tarea más permisiva que la regresión de $\Delta G$, ya que solo requiere
acertar el lado de un umbral definido sobre una combinación de las seis
afinidades de docking, afinidades que a su vez ya se predicen
razonablemente bien desde los descriptores moleculares.

Por último, el factor corrector por peso molecular del ITI () se
reconstruyó a partir de los valores observados en los resultados
guardados, ya que el script que lo generaba originalmente no quedó
versionado en el repositorio. La reconstrucción reproduce con fidelidad
el podio del ranking ITI (Triethylene glycol, Tetraethylenepentamine,
Pentaethylenehexamine) pero introduce incertidumbre en el orden exacto
de las moléculas siguientes y, por extensión, en qué moléculas concretas
caen a cada lado del umbral favorable/desfavorable del clasificador. Se
recomienda confirmar estas cifras si el script original se recupera.

Se ha construido una plataforma computacional completa para evaluar
componentes moleculares de nanopartículas poliméricas y dendrímeros
frente a seis proteínas que determinan cómo se distribuye y metaboliza
un fármaco en el organismo. Las conclusiones principales son:

1.  Las 58 estructuras moleculares se optimizaron a nivel cuántico
    mediante cálculos DFT sin que ninguna presentara frecuencias
    imaginarias. Los perfiles de polaridad superficial (COSMO-RS) y las
    frecuencias de vibración IR se calcularon para todas las moléculas,
    incluidos los casos especiales como el cisplatino y el condroitín
    sulfato.

2.  El protocolo de docking generó 342 resultados válidos frente a seis
    proteínas. La matriz de energías de unión resultante separa
    claramente fármacos, oligómeros y ligandos de targeting, con valores
    medios entre $-3,54$ y $-5,70$ kcal/mol según la proteína.

3.  El análisis de interacciones con PLIP identificó 708 contactos no
    covalentes en 15 complejos. Los residuos clave de FR$\alpha$
    (SER174, ARG103) y P-gp (ILE730, LEU312, LEU757) coinciden con los
    datos cristalográficos publicados, lo que valida las poses de
    docking obtenidas.

4.  Los modelos de regresión superaron R$^2 = 0,7$ en tres proteínas:
    P-gp (R$^2 = 0,836$, SVR, MAE = $0,371$ kcal/mol), TfR1 ($0,857$,
    Gradient Boosting) y Lisozima ($0,769$, SVR). El conjunto final de
    1097 descriptores (tras eliminar columnas constantes o con datos
    ausentes de un total de 4179 calculados) es suficientemente
    informativo para predecir la afinidad de unión con este dataset.

5.  El ITI identifica Triethylene glycol (100), Tetraethylenepentamine
    ($87,3$) y Pentaethylenehexamine ($85,0$) como los candidatos más
    favorables, seguidos por Chitotriose ($84,7$) y Chitobiose ($82,7$),
    cuya posición alta está además respaldada de forma independiente por
    el análisis PLIP. Todos son componentes habituales en formulaciones
    de nanopartículas clínicas.

6.  El análisis ADMET confirma que el 80% de los fármacos modelo son
    sustratos de P-gp frente al 0% de los monómeros de degradación. Ese
    contraste cuantifica el problema que resuelven los nanocarriers con
    ITI alto: al encapsular el fármaco, se reduce la exposición directa
    a la bomba de expulsión.

7.  El descriptor BertzCT, que mide la complejidad topológica de una
    molécula, tiene el mayor peso en el análisis SHAP (valor medio
    $0,27$), seguido de NumBonds ($0,11$) y MW ($0,11$). P-gp selecciona
    sus sustratos principalmente por tamaño y complejidad estructural,
    no por grupos químicos concretos.

8.  El clasificador binario alcanza AUC-ROC = $0,957$ y Accuracy =
    $0,927$ (Random Forest); un test de permutación de etiquetas
    confirma que esta capacidad discriminativa no es producto del azar
    (AUC medio con etiquetas aleatorias $= 0,52$). El espacio de
    descriptores moleculares separa razonablemente bien las moléculas
    favorables de las desfavorables para drug delivery.

9.  La herramienta Streamlit permite evaluar cualquier molécula nueva
    introduciendo su SMILES: calcula descriptores, predice el perfil de
    afinidad, estima el ITI y recomienda el carrier más adecuado sin
    necesidad de ejecutar ningún script.

El trabajo más evidente es ampliar el dataset. Con 200 a 500 moléculas,
incluyendo generaciones superiores de PAMAM (G3-G5) y más variedad de
ligandos de targeting, los modelos de regresión mejorarían notablemente.
Con ese volumen de datos también sería posible probar arquitecturas de
aprendizaje profundo como redes neuronales de grafos (GNN), que trabajan
directamente con la estructura molecular en lugar de con descriptores
numéricos calculados a mano.

La validación experimental es el paso que más valor daría a las
predicciones. Ensayos de inhibición de P-gp, estudios de captación
celular mediada por TfR1 o experimentos de biodistribución en animales
permitirían comparar las energías de unión predichas con datos reales.
Con esos datos también se podrían ajustar los pesos del ITI de forma
empírica, en lugar de definirlos con criterio experto, es decir,
asignados manualmente en función del conocimiento farmacológico sobre el
papel de cada proteína, sin ajuste sobre datos experimentales reales.

En cuanto al docking, aplicar métodos más precisos como FEP (,
perturbación de energía libre) o MM-GBSA (, mecánica molecular con
modelo de solvatación implícita) a los candidatos top del ITI daría una
caracterización energética más fiable. A diferencia de AutoDock Vina,
que estima la afinidad mediante una función de scoring rápida, estos
métodos calculan la energía libre de unión considerando explícitamente
la flexibilidad del receptor y los efectos del disolvente, a costa de un
coste computacional mucho mayor. Para un conjunto reducido de candidatos
son abordables.

La herramienta Streamlit tiene margen de mejora. Añadir un módulo de
dibujo de estructuras moleculares directamente en el navegador,
exportación de informes en PDF y conexión con bases de datos como
PubChem o ChEMBL haría que fuera más útil para investigadores que no
estén familiarizados con códigos SMILES.
