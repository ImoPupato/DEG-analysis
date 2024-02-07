**En este apartado nos propusimos realizar un análisis de la expresión de ciertos genes de nuestro interés en personas con Melanoma Cutánemo utilizando la información disponible en el proyecto TCGA-SKCM, diferenciando si provenían de muestras de metástasis o del tumor primario**.  
*Vale aclarar que en todos los casos se trabajó con RNASeq y las personas todas fueron diagnosticadas con melanoma cutáneo como tumor primario*.  

**Para comparar las expresiones (RNASeq) de las familias Vav y Sinucleina en muestras provenientes de metástasis y tumor primario, asignamos a cada muestra su proveniencia, siendo 101 correspondendientes a tumores primarios y 359 a metastásicos. Las expresiones fueron convertidas a CPM (*count per million*) con el paquete edgeR**.  

##### Ejemplo del código utilizado para la comparación
```R
t.test(fenotipos$VAV1.cpm[
  fenotipos$Type =="Metastasis"],
  fenotipos$VAV1.cpm[
    fenotipos$Type=="Primario"]
  ,paired = FALSE)
```
##### Ejemplo de la salida
```R
	Welch Two Sample t-test
data:  fenotipos$VAV1.cpm[fenotipos$Type == "Metastasis"] and fenotipos$VAV1.cpm[fenotipos$Type == "Primario"]
t = 7.8152, df = 425.81, p-value = 4.339e-14
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
  7.939098 13.274379
sample estimates:
mean of x mean of y 
17.410673  6.803934
```
#### **Tabla de resultados**  
Los valores de expresión de RNASeq están en *cuentas por millón* (cpm)

<div align="center">
	
|Gen| IC 95% | $\overline{x}$ metastasis | $\overline{x}$ primario | p-value |
|:---------------------------------:|:---------------------------------:|:---------------------------------:|:---------------------------------:|:---------------------------------:|
|VAV1| 7.94 - 13.27 |17.41 | 6.80 |4.34 x $10^{-14}$|
|VAV2| 2.49 - 10.84 | 33.15 | 26.49 |0.0019 |
|VAV3| 4.66 - 18.64 | 27.59 | 15.94 | 0.0012 |
|SNCA| -30.90  - -3.47| 80.73 |  97.91 | 0.0143 |
|SNCB|  -1.19 - 2.78 | 3.70  | 2.91 | 0.4297 |
|SNCG|  -1.84 - 3.26 | 4.56 | 3.85 | 0.5839 |

</div>  

Luego, calculamos la mediana (*valor de corte*) para generar cada grupo de contraste. 

<div align="center">
	
|Grupo|VAV1| VAV2|VAV3| SNCA|SNCB|SNCG|
|:-----------:|:--------:|:---------:|:-----------:|:----------:|:---------:|:----------:|
|Metástasis | 10.19  | 28.18 | 10.34 | 66.52 | 0.47 | 1.62 | 
|Primario | 4.24 | 22.09 | 7.17 | 92.93 | 0.57 | 1.50 | 

</div>  
Tomando en cuenta los valores de corte, los grupos quedaron conformados según expresión del Gen:
<div align="center">
	
|Grupo|Alta|Baja|
|:-----------:|:--------:|:---------:|
|Metástasis | 180  | 179 |
|Primario | 51 | 50 |

</div> 
Para el caso del análisis *completo* cada grupo de expresión, alta y baja, quedó conformado por 230 muestras.  
  
**Se llevaron adelante cuatro tipos de análisis: Sobrevida, Correlación con perfil inmune, ORA y GSEA**. 

### Sobrevida
Se determinaron los valores de corte (*en cpm*) óptimos correspondientes a un análisis con la función *surv_cutpoint* del paquete 'survminer'  

**Ejemplo del código utilizado para la obtención de los valores de corte**
```R
library(survminer)  #carga de librería
FenotiposP<-subset(Fenotipos,Type=="Primario") #selección de expresiones por fenotipo primario
FenotiposM<-subset(Fenotipos,Type=="Metastasis")#selección de expresiones por fenotipo metastásico

surv_rnaseq<-surv_cutpoint(
  Fenotipos*, # * sería P o M según el análisis
  time = "times", # tiempo al seguimiento
  event = "status", # estado de sobrevida, 1 o 0
  variables = c("VAV1.cpm","VAV2.cpm", "VAV3.cpm","SNCA.cpm","SNCB.cpm","SNCG.cpm")
)
summary(surv_rnaseq) # para extraer los valores de corte
```

<div align="center">

 |Gen| Completo| Metastasis | Primario |
|:-----------:|:--------:|:---------:|:-----------:|
VAV1| 7.72| 7.18 | 12.64 |
VAV2| 57.48 | 57.55 | 6.27 |
VAV3|  9.79 | 9.14 | 2.40 |
SNCA| 94.80 |101.75|63.55 |
SNCB|  0.034 |0.033|0.062 |
SNCG|  3.39 |3.39|0.52 |
</div>

**Ejemplo del código utilizado para el gráfico**

```
library(survival)
surv_rnaseq <- surv_categorize(surv_rnaseq)
fit <- survfit(Surv(times, status) ~ *, # * sería cada gen de contraste
               data = Fenotipos*) # * sería P o M según el análisis
ggsurvplot(
  fit,                     
  risk.table = TRUE,       
  pval = TRUE,             
  conf.int = FALSE,                  
  xlim = c(0,5500),      
  break.time.by = 1000,    
  ggtheme = theme_bw(), 
  risk.table.y.text.col = T, 
  font.legend=8, 
  legend.title = "Expression",
  risk.table.y.text = FALSE,                          
)
```

![Gráfica de sobrevida VAVs](https://github.com/ImoPupato/DEG-analysis/blob/main/Survival%20Plots%20VAVs.jpg)

![Gráfica de sobrevida SNCs](https://github.com/ImoPupato/DEG-analysis/blob/main/Survival%20Plots%20SNCs.jpg)

Se repitió el análisis anterior utilizando como valores de corte la mediana de cada una de las expresiones obteniendo los siguientes p-value:

<div align="center">

 |Gen| Completo| Metastasis | Primario |
|:-----------:|:--------:|:---------:|:-----------:|
VAV1| 0.0056| 0.0056 | 0.61 |
VAV2| 0.9| 0.85 | 0.67 |
VAV3|  0.0027 | 0.0018 | 0.67 |
SNCA|  0.076 |0.064|0.22 |
SNCB|  0.76| 0.72| 0.9 |
SNCG|  0.57 |0.49|0.34 |
</div>

### Correlación con perfil inmune
Se realizó un análisis de correlación utilizando la matriz inmune obtenida de [TIMER2.0](http://timer.comp-genomics.org/timer/)  

**Ejemplo del código utilizado para la comparación**
```R
library("corrplot")
library("ggstatsplot")
cor.test(Gene Expression,xCell Score, method="kendall") #Gene Expression corresponde al valor de expresión (*en cpm*) y xCell Score corresponde al valor proporcionado por TIMER2.0
```
El coeficiente utilizado es $\tau$, que varía entre -1 y 1 dependiendo del tipo de asociación lineal (inversa o directa) y fuerza (débil, moderada, fuerte). Entre paréntesis Se indica el p-value correspondiente al análisis.  

#### Vav1
|Set de análisis| Completo| Metastasis | Primario |
|:-----------:|:--------:|:---------:|:-----------:|
|Microambiente tumoral | 0.78 (2.2 x $10^{-16}$) |0.80 (2.2 x $10^{-16}$) | 0.65 (2.2 x $10^{-16}$)|
|Puntuación Inmune| 0.76 (2.2 x $10^{-16}$) |  0.78 (2.2 x $10^{-16}$)| 0.65 (2.2 x $10^{-16}$)|

#### Vav2
|Set de análisis| Completo| Metastasis | Primario |
|:-----------:|:--------:|:---------:|:-----------:|
|Microambiente tumoral | 0.027 (0.39) |0.018 (0.61) | -0.074 (0.28)|
|Puntuación Inmune| -0.018 (0.56) | -0.033 (0.35)|-0.077 (0.26)|

#### Vav3
|Set de análisis| Completo| Metastasis | Primario |
|:-----------:|:--------:|:---------:|:-----------:|
|Microambiente tumoral | 0.20 (2.36 x $10^{-10}$) |0.23 (1.08 x $10^{-10}$) | -0.017 (0.80)|
|Puntuación Inmune| 0.19 (9.22 x $10^{-10}$) | 0.22 (2.31 x $10^{-10}$)|-0.021 (0.75)|

#### SNCA
|Set de análisis| Completo| Metastasis | Primario |
|:-----------:|:--------:|:---------:|:-----------:|
|Microambiente tumoral | -0.16 (2.97 x $10^{-5}$) |-0.15 (2.97 x $10^{-5}$) | -0.15 (0.03)|
|Puntuación Inmune| -0.11 (3.74 x $10^{-4}$) | -0.097 (6.07 x $10^{-3}$)|-0.095 (0.16)|

#### SNCB
|Set de análisis| Completo| Metastasis | Primario |
|:-----------:|:--------:|:---------:|:-----------:|
|Microambiente tumoral | -0.077 (0.013) |-0.070 (0.049) | -0.038 (0.58)|
|Puntuación Inmune| -0.084 (0.0069) | -0.073 (0.039)|-0.065 (0.34)|

#### SNCG
|Set de análisis| Completo| Metastasis | Primario |
|:-----------:|:--------:|:---------:|:-----------:|
|Microambiente tumoral | 0.19 (4.81 x $10^{-10}$) |0.22 (8.02 x $10^{-10}$) | -0.038 (0.58)|
|Puntuación Inmune| 0.13 (3.02 x $10^{-5}$) | 0.15 (1.18 x $10^{-5}$)|-0.065 (0.34)|

A través del siguiente gráfico de correlación, podemos visualizar la fuerza de asociación entre los tipos celulares y cada expresión génica según set de análisis. Se utilizó la matriz correspondiente al algoritmo xCell:  

  
![Gráfico de correlación](https://github.com/ImoPupato/DEG-analysis/blob/main/Correlacion.jpg)

### ORA
Primero se llevó adelante un análisis de la expresión diferencial para cada contraste, los resultados de los *exact tests* de comparación son los siguientes:
- [Metástasis vs Primario](https://drive.google.com/file/d/1C4Guj3ZnUNzbCbqUPLPuF3dJgMLkn83N/view?usp=drive_link)
- Expresión de Vav1: [Completo](https://drive.google.com/file/d/1CIADvX8z7qRZzBCFrVh-b1eSVgJxdLCg/view?usp=drive_link), [Metastasis](https://drive.google.com/file/d/19C5vIoan3DWFPivaBssyoBNwVNZlzAz8/view?usp=drive_link), [Primario](https://drive.google.com/file/d/19G1d4nvvGEkTxGfqmB-HXtMOHJ2_9nm_/view?usp=drive_link)
- Expresión de Vav2: [Completo](https://drive.google.com/file/d/1CJnIBh708jfYfYAEp2ecvF9f9DQ4Hu3W/view?usp=drive_link), [Metastasis](https://drive.google.com/file/d/197y1I_elcBIkaAxXii9ueV5IjvXQVqOK/view?usp=drive_link), [Primario](https://drive.google.com/file/d/19IlufgTtkd-ChjXTz9ENxIhezy4ODQdr/view?usp=drive_link)
- Expresión de Vav3: [Completo](https://drive.google.com/file/d/1ClII0RRga8zBmidsgOYzQGgst2zoXcqe/view?usp=drive_link), [Metastasis](https://drive.google.com/file/d/199BoIXRLJUFlVadnbh34NuwPiE3IluG2/view?usp=sharing), [Primario](https://drive.google.com/file/d/19N6iRuVa2MPQHX5W1Ewi_4eNaW4a0Qfa/view?usp=drive_link)
- Expresión de SNCA: [Completo](https://drive.google.com/file/d/1CRIZYRWOe7NnOh-P7CAExus5by4JwOTw/view?usp=drive_link),[Metastasis](https://drive.google.com/file/d/19Ap6hGRPotYPr6Bz3HcwTzawenqfZdwq/view?usp=drive_link),[Primario](https://drive.google.com/file/d/19KlevtqKbVq2Ozevu8DkI25CLKDRTSaH/view?usp=drive_link)
- Expresión de SNCB: [Completo](https://drive.google.com/file/d/1CohMTbu61uYgrk_tyx92lh9d8aBkSP3k/view?usp=drive_link),[Metastasis](https://drive.google.com/file/d/19BBj2_gRtg47vOKEX2-1MWBRzUjBQdHz/view?usp=drive_link),[Primario](https://drive.google.com/file/d/19NwmKVOBD6Qf3uuPmLNGROhbGP7oqTsk/view?usp=drive_link)
- Expresión de SNCG: [Completo](https://drive.google.com/file/d/1CohNMa9hfH1JrycAZat3fPUrlLRldjWv/view?usp=drive_link),[Metastasis](https://drive.google.com/file/d/19Bo3lc3lrY-Vj5-nbOd1VGhbrG4BBhRF/view?usp=drive_link),[Primario](https://drive.google.com/file/d/19Rp0SVsoTL4OBFxsuuaKbJ_vbM3nmVJh/view?usp=drive_link)
  
Los parámetros utilizados fueron: LogFC>1 y FDR<0.01.  

Los gráficos de enriquecimiento se encuentran en el siguiente [link](https://drive.google.com/drive/folders/1ENIiLyNCG1wzut2m4cEzzHYxRxlw7JL4?usp=sharing)

  
En el caso del contraste de expresión alta vs baja de SNCB para el grupo de muestras provenientes de tumores primarios, no se observaron vias eniquecidas por los genes diferencialmente expresados.

### GSEA
Se generaron las siguientes matrices:
- [Matriz de expresión análisis Completo o MvsP](https://drive.google.com/file/d/1TE0wsNRsJk0JiLgzSH_G1ka6et60kKB_/view?usp=drive_link)
- [Matriz de expresión análisis expresión de genes de interes en muestras provenientes de tejido tumoral metastásico](https://drive.google.com/file/d/1_ZK_z0lGFaCkXZgYe6DH6b2WYyKztE6t/view?usp=drive_link)
- [Matriz de expresión análisis expresión de genes de interes en muestras provenientes de tejido tumoral primario](https://drive.google.com/file/d/1_H80pHoUc4c9PdDOTFMcMm9LGmWRuTSx/view?usp=drive_link)

Los fenotipos de contraste son:
- [Metástasis vs Primario](https://drive.google.com/file/d/1TE0wsNRsJk0JiLgzSH_G1ka6et60kKB_/view?usp=drive_link)
- Expresión de Vav1: [Completo](https://drive.google.com/file/d/1_igTwl-hgYfYzTiNlcIN8p61Q8QhQTfb/view?usp=drive_link), [Metastasis](https://drive.google.com/file/d/1_aT2Mj34gWdk4tJut97PoR_OxVgspjVb/view?usp=drive_link), [Primario](https://drive.google.com/file/d/1_L1zmKlKLINVPx1EX9WGgWd9SrbAryE9/view?usp=drive_link)
- Expresión de Vav2: [Completo](https://drive.google.com/file/d/1_jUzpTvgPC3ESDEGV5UL9PONVpAINtHj/view?usp=drive_link), [Metastasis](https://drive.google.com/file/d/1_c7ZFpLVxqsdr-JKmUoc9DOObfUDCu5B/view?usp=drive_link), [Primario](https://drive.google.com/file/d/1_X2TNGNE109nRLOFoTqRwTXDqmmKguGo/view?usp=drive_link)
- Expresión de Vav3: [Completo](https://drive.google.com/file/d/1_mRCUAPjYrvzWccvJwwD5NcpUy87jFYk/view?usp=drive_link), [Metastasis](https://drive.google.com/file/d/1_daAS897b7kYBZXkUn19ZaovEgxsrDNL/view?usp=drive_link), [Primario](https://drive.google.com/file/d/1_LGZdctT86jekzC5B--kXM5zMr1qKZ1g/view?usp=drive_link)
- Expresión de SNCA: [Completo](https://drive.google.com/file/d/1xZgsqepJ-sGsg7U-WCc3R3NeD8PWFnpD/view?usp=drive_link), [Metastasis](https://drive.google.com/file/d/1_dzjdtkNzaw7lMbOY4kSsadKhFKiMCyE/view?usp=drive_link), [Primario](https://drive.google.com/file/d/1_LneeL2rS02Uuz1MEQu-RYlIM-Ht_N6z/view?usp=drive_link)
- Expresión de SNCB: [Completo](https://drive.google.com/file/d/1_qlTwn2nLClX35w3e_v-xKND5V6-sKvz/view?usp=drive_link), [Metastasis](https://drive.google.com/file/d/1_endX2CGNLx0A4nKGFP_eWP9Op4H243s/view?usp=drive_link), [Primario](https://drive.google.com/file/d/1_UAkkxjCGXzOxitAagC12vPVeNrwKPGv/view?usp=drive_link)
- Expresión de SNCG: [Completo](https://drive.google.com/file/d/1_qsL2sMC2HDmmERVtuaeC29tMoOxDu_B/view?usp=drive_link), [Metastasis](https://drive.google.com/file/d/1_ezJTaxVUrLI-UcsrhxcqXNVEF5dgHI3/view?usp=drive_link), [Primario](https://drive.google.com/file/d/1_WFZ50xb8bYjEMWG_xPfoSjmIGn0S5z5/view?usp=drive_link)
  
Las salidas de GSEA de este análisis están disponibles en:
- [Metástasis vs Primario](https://drive.google.com/file/d/1_vl1J-aVshIkYyvoVo1pP2dNJSXSFxo4/view?usp=sharing)
- Expresión de Vav1: [Completo](https://drive.google.com/file/d/1lD2mM6bwc-qy3rDLO_6LKzTxWuu7GJh0/view?usp=drive_link), [Metastasis](https://drive.google.com/file/d/192uFvj194FxrqMJVtehVyiK9BKR9Fl0h/view?usp=sharing), [Primario](https://drive.google.com/file/d/196rG9jW5C2dLNsPQaKZ0lHXTw0Ryqfte/view?usp=sharing)
- Expresión de Vav2: [Completo](https://drive.google.com/file/d/1qnTv26THS800Kk1iJ2pVtvGyAL8_VD0i/view?usp=drive_link), [Metastasis](https://drive.google.com/file/d/193WkZ8RZMG_E06MKXOIji40OYFCZJz8O/view?usp=sharing), [Primario](https://drive.google.com/file/d/197N6hMANlz0bSG6SmsfFwvDqRhvZ0tQO/view?usp=drive_link)
- Expresión de Vav3: [Completo](https://drive.google.com/file/d/1VgkTZshyHJdtKzdBIJCHxhHHsdACR4Wj/view?usp=drive_link), [Metastasis](https://drive.google.com/file/d/193Xno3LGx83ax70YP2L_j3S-px7pvqIR/view?usp=sharing), [Primario](https://drive.google.com/file/d/1940Y27GFJ8c0nYGSy90ceW-cANNd213F/view?usp=drive_link)
- Expresión de SNCA: [Completo](https://drive.google.com/file/d/1xZgsqepJ-sGsg7U-WCc3R3NeD8PWFnpD/view?usp=drive_link), [Metastasis](https://drive.google.com/file/d/191_a_KLDzVxIuZo_U4LEviVpC90601Fu/view?usp=sharing), [Primario](https://drive.google.com/file/d/19631KSQgKSQWgcFkfZjcMHY_sucblhEU/view?usp=drive_link)
- Expresión de SNCB: [Completo](https://drive.google.com/file/d/1wn2yCmyW5hb0RIRFRUaMs_NSq63_fOVB/view?usp=drive_link), [Metastasis](https://drive.google.com/file/d/191MqzDOEF_72bEun1-A8HTcl6V5Glc9G/view?usp=sharing), [Primario](https://drive.google.com/file/d/1967KRZD2NwVzci0as58eiMPe15RtlUqI/view?usp=drive_link)
- Expresión de SNCG: [Completo](https://drive.google.com/file/d/1Uyk55i3dQ0sI6rs19CpxlNQASmowH95i/view?usp=drive_link), [Metastasis](https://drive.google.com/file/d/190D_S5d3R1T0R7XMHVci7CnXQVhB5QnI/view?usp=sharing), [Primario](https://drive.google.com/file/d/196T3oDDSPgtTzFV1zwC4EcIPh26j0sQR/view?usp=drive_link)  

A continuación se muestra el resumen de las vías significaticas (FWER p-val<0.01) de cada contraste:  
#### Vav1 Completo:
<div align="center">
	
|Orden |Via Reactome |NES |FWER p-val |
|:-----------:|:--------:|:---------:|:---------:|
|	1	|	REACTOME_IMMUNOREGULATORY_INTERACTIONS_BETWEEN_A_LYMPHOID_AND_A_NON_LYMPHOID_CELL	|	2.90	 |	0	|
|	2	|	REACTOME_COMPLEMENT_CASCADE	|	2.64	 |	0	|
|	3	|	REACTOME_INITIAL_TRIGGERING_OF_COMPLEMENT	|	2.61	 |	0	|
|	4	|	REACTOME_CREATION_OF_C4_AND_C2_ACTIVATORS	|	2.59	 |	0	|
|	5	|	REACTOME_FCGR_ACTIVATION	|	2.58	 |	0	|
|	6	|	REACTOME_INTERFERON_GAMMA_SIGNALING	|	2.57	 |	0	|
|	7	|	REACTOME_ANTIGEN_ACTIVATES_B_CELL_RECEPTOR_BCR_LEADING_TO_GENERATION_OF_SECOND_MESSENGERS	|	2.56	 |	0	|
|	8	|	REACTOME_FCGR3A_MEDIATED_IL10_SYNTHESIS	|	2.53	 |	0	|
|	9	|	REACTOME_CHEMOKINE_RECEPTORS_BIND_CHEMOKINES	|	2.53	 |	0	|
|	10	|	REACTOME_PARASITE_INFECTION	|	2.51	 |	0	|
|	11	|	REACTOME_FCERI_MEDIATED_CA_2_MOBILIZATION	|	2.50	 |	0	|
|	12	|	REACTOME_SCAVENGING_OF_HEME_FROM_PLASMA	|	2.50	 |	0	|
|	13	|	REACTOME_FCERI_MEDIATED_MAPK_ACTIVATION	|	2.50	 |	0	|
|	14	|	REACTOME_CD22_MEDIATED_BCR_REGULATION	|	2.48	 |	0	|
|	15	|	REACTOME_BINDING_AND_UPTAKE_OF_LIGANDS_BY_SCAVENGER_RECEPTORS	|	2.46	 |	0	|
|	16	|	REACTOME_CELL_SURFACE_INTERACTIONS_AT_THE_VASCULAR_WALL	|	2.46	 |	0	|
|	17	|	REACTOME_ANTI_INFLAMMATORY_RESPONSE_FAVOURING_LEISHMANIA_PARASITE_INFECTION	|	2.45	 |	0	|
|	18	|	REACTOME_INTERFERON_ALPHA_BETA_SIGNALING	|	2.42	 |	0	|
|	19	|	REACTOME_GENERATION_OF_SECOND_MESSENGER_MOLECULES	|	2.42	 |	0	|
|	20	|	REACTOME_ROLE_OF_PHOSPHOLIPIDS_IN_PHAGOCYTOSIS	|	2.41	 |	0	|
|	21	|	REACTOME_ROLE_OF_LAT2_NTAL_LAB_ON_CALCIUM_MOBILIZATION	|	2.41	 |	0	|
|	22	|	REACTOME_FCGAMMA_RECEPTOR_FCGR_DEPENDENT_PHAGOCYTOSIS	|	2.41	 |	0	|
|	23	|	REACTOME_SIGNALING_BY_THE_B_CELL_RECEPTOR_BCR	|	2.40	 |	0	|
|	24	|	REACTOME_TCR_SIGNALING	|	2.40	 |	0	|
|	25	|	REACTOME_LEISHMANIA_INFECTION	|	2.39	 |	0	|
|	26	|	REACTOME_PD_1_SIGNALING	|	2.39	 |	0	|
|	27	|	REACTOME_INTERLEUKIN_10_SIGNALING	|	2.38	 |	0	|
|	28	|	REACTOME_INTERLEUKIN_2_FAMILY_SIGNALING	|	2.37	 |	0	|
|	29	|	REACTOME_COSTIMULATION_BY_THE_CD28_FAMILY	|	2.36	 |	0	|
|	30	|	REACTOME_SIGNALING_BY_INTERLEUKINS	|	2.34	 |	0	|
|	31	|	REACTOME_POTENTIAL_THERAPEUTICS_FOR_SARS	|	2.33	 |	0	|
|	32	|	REACTOME_CLASS_A_1_RHODOPSIN_LIKE_RECEPTORS	|	2.33	 |	0	|
|	33	|	REACTOME_INTERLEUKIN_4_AND_INTERLEUKIN_13_SIGNALING	|	2.31	 |	0	|
|	34	|	REACTOME_ANTIGEN_PROCESSING_CROSS_PRESENTATION	|	2.29	 |	0	|
|	35	|	REACTOME_FCERI_MEDIATED_NF_KB_ACTIVATION	|	2.29	 |	0	|
|	36	|	REACTOME_FC_EPSILON_RECEPTOR_FCERI_SIGNALING	|	2.28	 |	0	|
|	37	|	REACTOME_PEPTIDE_LIGAND_BINDING_RECEPTORS	|	2.27	 |	0	|
|	38	|	REACTOME_DAP12_INTERACTIONS	|	2.23	 |	0	|
|	39	|	REACTOME_INTEGRIN_CELL_SURFACE_INTERACTIONS	|	2.23	 |	0	|
|	40	|	REACTOME_NEUTROPHIL_DEGRANULATION	|	2.22	 |	0	|
|	41	|	REACTOME_INTERLEUKIN_3_INTERLEUKIN_5_AND_GM_CSF_SIGNALING	|	2.19	 |	0	|
|	42	|	REACTOME_INTERFERON_SIGNALING	|	2.15	 |	0	|
|	43	|	REACTOME_TOLL_LIKE_RECEPTOR_CASCADES	|	2.13	 |	0	|
|	44	|	REACTOME_TNFR2_NON_CANONICAL_NF_KB_PATHWAY	|	2.11	 |	0	|
|	45	|	REACTOME_GPCR_LIGAND_BINDING	|	02.09	 |	0	|
|	46	|	REACTOME_GPVI_MEDIATED_ACTIVATION_CASCADE	|	02.08	 |	0	|
|	47	|	REACTOME_INTERLEUKIN_12_FAMILY_SIGNALING	|	02.06	 |	0	|
|	48	|	REACTOME_THE_ROLE_OF_NEF_IN_HIV_1_REPLICATION_AND_DISEASE_PATHOGENESIS	|	02.05	 |	1	|
|	49	|	REACTOME_TNFS_BIND_THEIR_PHYSIOLOGICAL_RECEPTORS	|	02.05	 |	3	|
|	50	|	REACTOME_DAP12_SIGNALING	|	02.05	 |	3	|
|	51	|	REACTOME_INTERLEUKIN_RECEPTOR_SHC_SIGNALING	|	02.04	 |	3	|
|	52	|	REACTOME_SIGNAL_REGULATORY_PROTEIN_FAMILY_INTERACTIONS	|	02.03	 |	4	|
|	53	|	REACTOME_DECTIN_2_FAMILY	|	02.02	 |	4	|
|	54	|	REACTOME_G_ALPHA_I_SIGNALLING_EVENTS	|	02.02	 |	5	|
|	1	|	RRNA PROCESSING	|	-3.54	|	0	|
|	2	|	EUKARYOTIC TRANSLATION INITIATION	|	-3.43	|	0	|
|	3	|	EUKARYOTIC TRANSLATION ELONGATION	|	-3.43	|	0	|
|	4	|	SELENOAMINO ACID METABOLISM	|	-3.29	|	0	|
|	5	|	RESPONSE OF EIF2AK4 GCN2 TO AMINO ACID DEFICIENCY	|	-3.28	|	0	|
|	6	|	SRP DEPENDENT COTRANSLATIONAL PROTEIN TARGETING TO MEMBRANE	|	-3.25	|	0	|
|	7	|	NONSENSE MEDIATED DECAY NMD	|	-3.24	|	0	|
|	8	|	TRANSLATION	|	-3.21	|	0	|
|	9	|	ACTIVATION OF THE MRNA UPON BINDING OF THE CAP BINDING COMPLEX AND EIFS AND SUBSEQUENT BINDING TO 43S	|	-3.02	|	0	|
|	10	|	INFLUENZA INFECTION	|	-3	|	0	|
|	11	|	CELLULAR RESPONSE TO STARVATION	|	-2.99	|	0	|
|	12	|	MITOCHONDRIAL PROTEIN IMPORT	|	-2.8	|	0	|
|	13	|	SARS COV 1 MODULATES HOST TRANSLATION MACHINERY	|	-2.78	|	0	|
|	14	|	SARS COV 2 MODULATES HOST TRANSLATION MACHINERY	|	-2.74	|	0	|
|	15	|	MITOCHONDRIAL TRANSLATION	|	-2.73	|	0	|
|	16	|	REGULATION OF EXPRESSION OF SLITS AND ROBOS	|	-2.63	|	0	|
|	17	|	RESPIRATORY ELECTRON TRANSPORT	|	-2.62	|	0	|
|	18	|	RRNA MODIFICATION IN THE NUCLEUS AND CYTOSOL	|	-2.61	|	0	|
|	19	|	THE CITRIC ACID TCA CYCLE AND RESPIRATORY ELECTRON TRANSPORT	|	-2.51	|	0	|
|	20	|	RESPIRATORY ELECTRON TRANSPORT ATP SYNTHESIS BY CHEMIOSMOTIC COUPLING AND HEAT PRODUCTION BY UNCOUPLING PROTEINS	|	-2.43	|	0	|
|	21	|	TRNA PROCESSING	|	-2.37	|	0	|
|	22	|	RNA POLYMERASE I TRANSCRIPTION INITIATION	|	-2.34	|	0	|
|	23	|	RNA POLYMERASE I PROMOTER ESCAPE	|	-2.29	|	0	|
|	24	|	COMPLEX I BIOGENESIS	|	-2.27	|	0	|
|	25	|	TRNA AMINOACYLATION	|	-2.2	|	0.001	|
|	26	|	TRNA MODIFICATION IN THE NUCLEUS AND CYTOSOL	|	-2.19	|	0.001	|
|	27	|	RNA POLYMERASE I TRANSCRIPTION TERMINATION	|	-2.17	|	0.001	|
|	28	|	MITOCHONDRIAL TRNA AMINOACYLATION	|	-2.16	|	0.002	|
|	29	|	MRNA SPLICING MINOR PATHWAY	|	-2.15	|	0.002	|
|	30	|	PROCESSING OF CAPPED INTRON CONTAINING PRE MRNA	|	-2.14	|	0.002	|
|	31	|	POSITIVE EPIGENETIC REGULATION OF RRNA EXPRESSION	|	-2.14	|	0.002	|
|	32	|	RNA POLYMERASE I TRANSCRIPTION	|	-2.13	|	0.002	|
|	33	|	METABOLISM OF AMINO ACIDS AND DERIVATIVES	|	-2.11	|	0.002	|
|	34	|	TRANSCRIPTIONAL REGULATION BY SMALL RNAS	|	-2.11	|	0.002	|
|	35	|	RNA POLYMERASE III TRANSCRIPTION INITIATION FROM TYPE 1 PROMOTER	|	-2.1	|	0.002	|
|	36	|	CRISTAE FORMATION	|	-2.1	|	0.002	|
|	37	|	PYRUVATE METABOLISM AND CITRIC ACID TCA CYCLE	|	-2.09	|	0.002	|
|	38	|	RNA POLYMERASE III TRANSCRIPTION	|	-2.08	|	0.004	|
|	39	|	MRNA SPLICING	|	-2.08	|	0.004	|
|	40	|	B WICH COMPLEX POSITIVELY REGULATES RRNA EXPRESSION	|	-2.07	|	0.005	|
|	41	|	SIGNALING BY ROBO RECEPTORS	|	-2.07	|	0.005	|
|	42	|	TP53 REGULATES METABOLIC GENES	|	-2.07	|	0.005	|
|	43	|	GENE SILENCING BY RNA	|	-2.07	|	0.006	|
|	44	|	RNA POLYMERASE III TRANSCRIPTION INITIATION FROM TYPE 3 PROMOTER	|	-2.06	|	0.009	|
</div> 

[Heatmap Vav1 (C)] (https://drive.google.com/file/d/11JSiNKhYOnGe3KVDchBO4zjynxQgrpCv/view?usp=sharing)  


#### Vav2 Completo:
<div align="center">
	
|Orden |Via Reactome |NES |FWER p-val |
|:-----------:|:--------:|:---------:|:---------:|
1	|	NCAM1 INTERACTIONS	|	2.22	|	0	|
2	|	ELASTIC FIBRE FORMATION	|	2.21	|	0	|
3	|	MOLECULES ASSOCIATED WITH ELASTIC FIBRES	|	2.18	|	0	|
4	|	ECM PROTEOGLYCANS	|	2.12	|	0.001	|
5	|	NRAGE SIGNALS DEATH THROUGH JNK	|	2.12	|	0.001	|
6	|	NCAM SIGNALING FOR NEURITE OUT GROWTH	|	2.12	|	0.001	|
7	|	COLLAGEN BIOSYNTHESIS AND MODIFYING ENZYMES	|	2.1	|	0.001	|
8	|	COLLAGEN CHAIN TRIMERIZATION	|	2.1	|	0.001	|
9	|	INTEGRIN CELL SURFACE INTERACTIONS	|	2.07	|	0.001	|
10	|	RUNX2 REGULATES BONE DEVELOPMENT	|	2.07	|	0.001	|
11	|	O GLYCOSYLATION OF TSR DOMAIN CONTAINING PROTEINS	|	2.07	|	0.001	|
12	|	ASSEMBLY OF COLLAGEN FIBRILS AND OTHER MULTIMERIC STRUCTURES	|	2.05	|	0.001	|
13	|	SIGNAL TRANSDUCTION BY L1	|	2.04	|	0.001	|
14	|	NON INTEGRIN MEMBRANE ECM INTERACTIONS	|	2.03	|	0.001	|
15	|	EPHB MEDIATED FORWARD SIGNALING	|	2.01	|	0.001	|
16	|	CROSSLINKING OF COLLAGEN FIBRILS	|	2.01	|	0.001	|
17	|	G ALPHA 12 13 SIGNALLING EVENTS	|	2.01	|	0.001	|
18	|	LAMININ INTERACTIONS	|	1.99	|	0.002	|
19	|	EXTRACELLULAR MATRIX ORGANIZATION	|	1.98	|	0.003	|
20	|	COLLAGEN FORMATION	|	1.97	|	0.004	|
21	|	DISEASES ASSOCIATED WITH O GLYCOSYLATION OF PROTEINS	|	1.96	|	0.008	|
1	|	METABOLISM OF AMINO ACIDS AND DERIVATIVES	|	---	|	0	|
2	|	NEUTROPHIL DEGRANULATION	|	---	|	0	|
3	|	OLFACTORY SIGNALING PATHWAY	|	---	|	0	|
4	|	RESPIRATORY ELECTRON TRANSPORT	|	-3.7679663	|	0	|
5	|	RESPIRATORY ELECTRON TRANSPORT ATP SYNTHESIS BY CHEMIOSMOTIC COUPLING AND HEAT PRODUCTION BY UNCOUPLING PROTEINS	|	-3.635556	|	0	|
6	|	THE CITRIC ACID TCA CYCLE AND RESPIRATORY ELECTRON TRANSPORT	|	-3.3615346	|	0	|
7	|	KERATINIZATION	|	-3.2431023	|	0	|
8	|	COMPLEX I BIOGENESIS	|	-3.1194296	|	0	|
9	|	CRISTAE FORMATION	|	-3.0371633	|	0	|
10	|	AMINO ACIDS REGULATE MTORC1	|	-2.9092202	|	0	|
11	|	MITOCHONDRIAL FATTY ACID BETA OXIDATION	|	-2.8600173	|	0	|
12	|	FORMATION OF THE CORNIFIED ENVELOPE	|	-2.6595182	|	0	|
13	|	CHOLESTEROL BIOSYNTHESIS	|	-2.606474	|	0	|
14	|	INSULIN RECEPTOR RECYCLING	|	-2.571626	|	0	|
15	|	GLUCONEOGENESIS	|	-2.5386863	|	0	|
16	|	TRANSFERRIN ENDOCYTOSIS AND RECYCLING	|	-2.534894	|	0	|
17	|	GLYCOSPHINGOLIPID CATABOLISM	|	-2.4508731	|	0	|
18	|	PYRUVATE METABOLISM AND CITRIC ACID TCA CYCLE	|	-2.3947058	|	0.003	|
19	|	MITOCHONDRIAL PROTEIN IMPORT	|	-2.3777003	|	0.004	|
20	|	IRON UPTAKE AND TRANSPORT	|	-2.3737867	|	0.004	|
21	|	MITOCHONDRIAL TRANSLATION	|	-2.2982318	|	0.008	|
22	|	TP53 REGULATES METABOLIC GENES	|	-2.2718391	|	0.009	|
23	|	CITRIC ACID CYCLE TCA CYCLE	|	-2.2347806	|	0.01	|
24	|	ANTIGEN PRESENTATION FOLDING ASSEMBLY AND PEPTIDE LOADING OF CLASS I MHC	|	-2.1774726	|	0.018	|
25	|	GLYCOSPHINGOLIPID METABOLISM	|	-2.177357	|	0.018	|
26	|	PROTEIN LOCALIZATION	|	-2.172199	|	0.018	|
27	|	BRANCHED CHAIN AMINO ACID CATABOLISM	|	-2.1617498	|	0.02	|
28	|	ROS AND RNS PRODUCTION IN PHAGOCYTES	|	-2.0883758	|	0.039	|
29	|	PYRUVATE METABOLISM	|	-2.0186293	|	0.069	|
</div> 

[Heatmap Vav2 (C)] (https://drive.google.com/file/d/1K1DrEFhyb-BD8vVETTZu-Z0ZOPsu72p9/view?usp=drive_link)  


#### Vav3 Completo:
<div align="center">
	
|Orden |Via Reactome |NES |FWER p-val |
|:-----------:|:--------:|:---------:|:---------:|
|	1	|	GENERATION OF SECOND MESSENGER MOLECULES	|	2.54	|	0	|
|	2	|	IMMUNOREGULATORY INTERACTIONS BETWEEN A LYMPHOID AND A NON LYMPHOID CELL	|	2.35	|	0	|
|	3	|	PD 1 SIGNALING	|	2.14	|	0	|
|	4	|	CHEMOKINE RECEPTORS BIND CHEMOKINES	|	2.14	|	0	|
|	5	|	INTERFERON GAMMA SIGNALING	|	2.11	|	0	|
|	6	|	INTERLEUKIN 2 FAMILY SIGNALING	|	2.05	|	0	|
|	7	|	DAP12 INTERACTIONS	|	2.04	|	0	|
|	8	|	ANTIGEN ACTIVATES B CELL RECEPTOR BCR LEADING TO GENERATION OF SECOND MESSENGERS	|	2.03	|	0	|
|	9	|	INITIAL TRIGGERING OF COMPLEMENT	|	2.02	|	0	|
|	10	|	TCR SIGNALING	|	2.01	|	0	|
|	11	|	COSTIMULATION BY THE CD28 FAMILY	|	2	|	0	|
|	12	|	DAP12 SIGNALING	|	2	|	0	|
|	13	|	FCGR ACTIVATION	|	1.98	|	0.003	|
|	14	|	COMPLEMENT CASCADE	|	1.96	|	0.004	|
|	15	|	INTERLEUKIN 10 SIGNALING	|	1.96	|	0.004	|
|	1	|	NEURONAL SYSTEM	|	---	|	0	|
|	2	|	GLUCONEOGENESIS	|	-2.67777	|	0	|
|	3	|	FORMATION OF TUBULIN FOLDING INTERMEDIATES BY CCT TRIC	|	-2.3860688	|	0.002	|
|	4	|	MITOCHONDRIAL PROTEIN IMPORT	|	-2.3535988	|	0.002	|
|	5	|	MITOCHONDRIAL TRANSLATION	|	-2.2560556	|	0.007	|
|	6	|	DISEASES ASSOCIATED WITH N GLYCOSYLATION OF PROTEINS	|	-2.255748	|	0.007	|
|	7	|	ABORTIVE ELONGATION OF HIV 1 TRANSCRIPT IN THE ABSENCE OF TAT	|	-2.2518127	|	0.007	|
|	8	|	POST CHAPERONIN TUBULIN FOLDING PATHWAY	|	-2.2264192	|	0.01	|
|	9	|	RECYCLING PATHWAY OF L1	|	-2.1985188	|	0.011	|
|	10	|	CELL EXTRACELLULAR MATRIX INTERACTIONS	|	-2.196962	|	0.011	|
|	11	|	FORMATION OF THE BETA CATENIN TCF TRANSACTIVATING COMPLEX	|	-2.158657	|	0.021	|
|	12	|	GLUCOSE METABOLISM	|	-2.1556664	|	0.022	|
|	13	|	HIV ELONGATION ARREST AND RECOVERY	|	-2.1194372	|	0.032	|
|	14	|	DISEASES OF CARBOHYDRATE METABOLISM	|	-2.1001415	|	0.037	|
|	15	|	SYNAPTIC ADHESION LIKE MOLECULES	|	-2.0822177	|	0.044	|
|	16	|	TRANSLATION	|	-2.0541263	|	0.058	|
|	17	|	GLYCOSPHINGOLIPID METABOLISM	|	-2.0055656	|	0.097	|
</div> 

[Heatmap Vav3 (C)] (https://drive.google.com/file/d/1B384QgIamWDQTrwmukeWHmC494BAMsYn/view?usp=drive_link)  

#### SNCA Completo:
<div align="center">
	
|Orden |Via Reactome |NES |FWER p-val |
|:-----------:|:--------:|:---------:|:---------:|
|	1	|	THE CITRIC ACID TCA CYCLE AND RESPIRATORY ELECTRON TRANSPORT		|	3.11	|	0	|
|	2	|	RESPIRATORY ELECTRON TRANSPORT		|	2.92	|	0	|
|	3	|	RESPIRATORY ELECTRON TRANSPORT ATP SYNTHESIS BY CHEMIOSMOTIC COUPLING AND HEAT PRODUCTION BY UNCOUPLING PROTEINS		|	2.89	|	0	|
|	4	|	MITOCHONDRIAL TRANSLATION		|	2.84	|	0	|
|	5	|	AMINO ACIDS REGULATE MTORC1		|	2.84	|	0	|
|	6	|	MITOCHONDRIAL FATTY ACID BETA OXIDATION		|	2.82	|	0	|
|	7	|	CRISTAE FORMATION		|	2.79	|	0	|
|	8	|	MITOCHONDRIAL PROTEIN IMPORT		|	2.69	|	0	|
|	9	|	TP53 REGULATES METABOLIC GENES		|	2.66	|	0	|
|	10	|	PYRUVATE METABOLISM AND CITRIC ACID TCA CYCLE		|	2.64	|	0	|
|	11	|	HOMOLOGOUS DNA PAIRING AND STRAND EXCHANGE		|	2.55	|	0	|
|	12	|	PROTEIN LOCALIZATION		|	2.55	|	0	|
|	13	|	SYNTHESIS OF PIPS AT THE EARLY ENDOSOME MEMBRANE		|	2.51	|	0	|
|	14	|	GOLGI ASSOCIATED VESICLE BIOGENESIS		|	2.48	|	0	|
|	15	|	DISEASES OF DNA REPAIR		|	2.44	|	0	|
|	16	|	TBC RABGAPS		|	2.43	|	0	|
|	17	|	COMPLEX I BIOGENESIS		|	2.43	|	0	|
|	18	|	HDR THROUGH HOMOLOGOUS RECOMBINATION HRR		|	2.39	|	0	|
|	19	|	GLYCOSPHINGOLIPID METABOLISM		|	2.37	|	0	|
|	20	|	METABOLISM OF COFACTORS		|	2.37	|	0	|
|	21	|	GLYCOSPHINGOLIPID CATABOLISM		|	2.34	|	0	|
|	22	|	AUTOPHAGY		|	2.34	|	0	|
|	23	|	HOMOLOGY DIRECTED REPAIR		|	2.34	|	0	|
|	24	|	RETROGRADE TRANSPORT AT THE TRANS GOLGI NETWORK		|	2.33	|	0	|
|	25	|	RAB REGULATION OF TRAFFICKING		|	2.33	|	0	|
|	26	|	RESOLUTION OF D LOOP STRUCTURES THROUGH SYNTHESIS DEPENDENT STRAND ANNEALING SDSA		|	2.32	|	0	|
|	27	|	S PHASE		|	2.32	|	0	|
|	28	|	ACTIVATION OF ATR IN RESPONSE TO REPLICATION STRESS		|	2.31	|	0	|
|	29	|	RESOLUTION OF AP SITES VIA THE MULTIPLE NUCLEOTIDE PATCH REPLACEMENT PATHWAY		|	2.31	|	0	|
|	30	|	DNA REPAIR		|	2.31	|	0.001	|
|	31	|	MTOR SIGNALLING		|	2.3	|	0.001	|
|	32	|	INSULIN RECEPTOR RECYCLING		|	2.29	|	0.002	|
|	33	|	CITRIC ACID CYCLE TCA CYCLE		|	2.29	|	0.002	|
|	34	|	MITOCHONDRIAL TRNA AMINOACYLATION		|	2.27	|	0.004	|
|	35	|	G2 M CHECKPOINTS		|	2.25	|	0.005	|
|	36	|	TRANSLOCATION OF SLC2A4 GLUT4 TO THE PLASMA MEMBRANE		|	2.25	|	0.005	|
|	37	|	RESOLUTION OF D LOOP STRUCTURES		|	2.24	|	0.008	|
|	1	|	COLLAGEN BIOSYNTHESIS AND MODIFYING ENZYMES		|	-2.63	|	0	|
|	2	|	COLLAGEN CHAIN TRIMERIZATION		|	-2.6	|	0	|
|	3	|	EXTRACELLULAR MATRIX ORGANIZATION		|	-2.58	|	0	|
|	4	|	INTEGRIN CELL SURFACE INTERACTIONS		|	-2.55	|	0	|
|	5	|	IMMUNOREGULATORY INTERACTIONS BETWEEN A LYMPHOID AND A NON LYMPHOID CELL		|	-2.53	|	0	|
|	6	|	ECM PROTEOGLYCANS		|	-2.52	|	0	|
|	7	|	COLLAGEN FORMATION		|	-2.48	|	0	|
|	8	|	ELASTIC FIBRE FORMATION		|	-2.44	|	0	|
|	9	|	MOLECULES ASSOCIATED WITH ELASTIC FIBRES		|	-2.38	|	0	|
|	10	|	EUKARYOTIC TRANSLATION ELONGATION		|	-2.36	|	0	|
|	11	|	COLLAGEN DEGRADATION		|	-2.35	|	0	|
|	12	|	ASSEMBLY OF COLLAGEN FIBRILS AND OTHER MULTIMERIC STRUCTURES		|	-2.35	|	0	|
|	13	|	NCAM1 INTERACTIONS		|	-2.32	|	0	|
|	14	|	RESPONSE OF EIF2AK4 GCN2 TO AMINO ACID DEFICIENCY		|	-2.32	|	0	|
|	15	|	REGULATION OF INSULIN LIKE GROWTH FACTOR IGF TRANSPORT AND UPTAKE BY INSULIN LIKE GROWTH FACTOR BINDING PROTEINS IGFBPS		|	-2.26	|	0	|
|	16	|	GPCR LIGAND BINDING		|	-2.26	|	0	|
|	17	|	CLASS A 1 RHODOPSIN LIKE RECEPTORS		|	-2.24	|	0	|
|	18	|	GENERATION OF SECOND MESSENGER MOLECULES		|	-2.23	|	0	|
|	19	|	DEGRADATION OF THE EXTRACELLULAR MATRIX		|	-2.22	|	0	|
|	20	|	NON INTEGRIN MEMBRANE ECM INTERACTIONS		|	-2.21	|	0	|
|	21	|	CROSSLINKING OF COLLAGEN FIBRILS		|	-2.19	|	0	|
|	22	|	COMPLEMENT CASCADE		|	-2.18	|	0	|
|	23	|	INTERLEUKIN 10 SIGNALING		|	-2.16	|	0	|
|	24	|	SRP DEPENDENT COTRANSLATIONAL PROTEIN TARGETING TO MEMBRANE		|	-2.15	|	0	|
|	25	|	SARS COV 1 MODULATES HOST TRANSLATION MACHINERY		|	-2.12	|	0.002	|
|	26	|	O GLYCOSYLATION OF TSR DOMAIN CONTAINING PROTEINS		|	-2.11	|	0.002	|
|	27	|	EUKARYOTIC TRANSLATION INITIATION		|	-2.1	|	0.003	|
|	28	|	NONSENSE MEDIATED DECAY NMD		|	-2.1	|	0.003	|
|	29	|	BINDING AND UPTAKE OF LIGANDS BY SCAVENGER RECEPTORS		|	-2.1	|	0.003	|
|	30	|	SYNDECAN INTERACTIONS		|	-2.09	|	0.003	|
|	31	|	INITIAL TRIGGERING OF COMPLEMENT		|	-2.08	|	0.003	|
|	32	|	PD 1 SIGNALING		|	-2.08	|	0.003	|
</div> 

[Heatmap SNCA (C)] (https://drive.google.com/file/d/1BEQpdVt-4UHpX6qWw5QMyx9KKlAfU7wb/view?usp=drive_link)  

#### SNCB Completo:
<div align="center">
	
|Orden |Via Reactome|NES |FWER p-val |
|:-----------:|:--------:|:---------:|:---------:|
|	1	|	GLUCONEOGENESIS	|	2.77	|	0	|
|	2	|	RESPIRATORY ELECTRON TRANSPORT ATP SYNTHESIS BY CHEMIOSMOTIC COUPLING AND HEAT PRODUCTION BY UNCOUPLING PROTEINS	|	2.77	|	0	|
|	3	|	FORMATION OF THE CORNIFIED ENVELOPE	|	2.73	|	0	|
|	4	|	THE CITRIC ACID TCA CYCLE AND RESPIRATORY ELECTRON TRANSPORT	|	2.72	|	0	|
|	5	|	RESPIRATORY ELECTRON TRANSPORT	|	2.7	|	0	|
|	6	|	COMPLEX I BIOGENESIS	|	2.35	|	0.004	|
|	1	|	TRANSPORT OF MATURE MRNAS DERIVED FROM INTRONLESS TRANSCRIPTS	|	-2	|	0.009	|
|	2	|	SYNTHESIS OF PIPS AT THE EARLY ENDOSOME MEMBRANE	|	-2	|	0.009	|
|	3	|	SUMOYLATION OF DNA DAMAGE RESPONSE AND REPAIR PROTEINS	|	-1.98	|	0.016	|
|	4	|	TRANSPORT OF MATURE TRANSCRIPT TO CYTOPLASM	|	-1.95	|	0.037	|
|	5	|	CILIUM ASSEMBLY	|	-1.92	|	0.068	|
|	6	|	SARS COV 2 ACTIVATES MODULATES INNATE AND ADAPTIVE IMMUNE RESPONSES	|	-1.91	|	0.074	|
|	7	|	CARGO TRAFFICKING TO THE PERICILIARY MEMBRANE	|	-1.91	|	0.074	|
</div> 

[Heatmap SNCB (C)] (https://drive.google.com/file/d/1rFKhz7GTGC5rhykqE6vO7mtMjNXj18Gs/view?usp=drive_link)  

#### SNCG Completo:
<div align="center">
	
|Orden |Via Reactome|NES |FWER p-val |
|:-----------:|:--------:|:---------:|:---------:|
|	1	|	IMMUNOREGULATORY INTERACTIONS BETWEEN A LYMPHOID AND A NON LYMPHOID CELL	|	2.73	|	0|
|	2	|	KERATINIZATION	|	2.7	|	0|
|	3	|	FORMATION OF THE CORNIFIED ENVELOPE	|	2.66	|	0|
|	4	|	CLASS A 1 RHODOPSIN LIKE RECEPTORS	|	2.5	|	0|
|	5	|	COMPLEMENT CASCADE	|	2.42	|	0|
|	6	|	GPCR LIGAND BINDING	|	2.38	|	0|
|	7	|	GENERATION OF SECOND MESSENGER MOLECULES	|	2.36	|	0|
|	8	|	CHEMOKINE RECEPTORS BIND CHEMOKINES	|	2.34	|	0|
|	9	|	ASSEMBLY OF COLLAGEN FIBRILS AND OTHER MULTIMERIC STRUCTURES	|	2.32	|	0|
|	10	|	EXTRACELLULAR MATRIX ORGANIZATION	|	2.3	|	0|
|	11	|	INTERLEUKIN 10 SIGNALING	|	2.28	|	0|
|	12	|	INITIAL TRIGGERING OF COMPLEMENT	|	2.28	|	0|
|	13	|	ECM PROTEOGLYCANS	|	2.27	|	0
|	14	|	INTEGRIN CELL SURFACE INTERACTIONS	|	2.26	|	0|
|	15	|	COLLAGEN CHAIN TRIMERIZATION	|	2.25	|	0|
|	16	|	PEPTIDE LIGAND BINDING RECEPTORS	|	2.24	|	0|
|	17	|	BINDING AND UPTAKE OF LIGANDS BY SCAVENGER RECEPTORS	|	2.21	|	0|
|	18	|	COLLAGEN FORMATION	|	2.2	|	0|
|	19	|	ELASTIC FIBRE FORMATION	|	2.2	|	0|
|	20	|	COLLAGEN BIOSYNTHESIS AND MODIFYING ENZYMES	|	2.17	|	0.001|
|	21	|	FORMATION OF FIBRIN CLOT CLOTTING CASCADE	|	2.17	|	0.001|
|	22	|	CROSSLINKING OF COLLAGEN FIBRILS	|	2.16	|	0.002|
|	23	|	MOLECULES ASSOCIATED WITH ELASTIC FIBRES	|	2.15	|	0.002|
|	24	|	INTERLEUKIN 4 AND INTERLEUKIN 13 SIGNALING	|	2.14	|	0.002|
|	25	|	NCAM1 INTERACTIONS	|	2.12	|	0.003|
|	26	|	TNFS BIND THEIR PHYSIOLOGICAL RECEPTORS	|	2.12	|	0.003|
|	27	|	COLLAGEN DEGRADATION	|	2.11	|	0.003|
|	28	|	CREATION OF C4 AND C2 ACTIVATORS	|	2.11	|	0.004|
|	29	|	TNF RECEPTOR SUPERFAMILY TNFSF MEMBERS MEDIATING NON CANONICAL NF KB PATHWAY	|	2.11	|	0.004|
|	30	|	DEGRADATION OF THE EXTRACELLULAR MATRIX	|	2.09	|	0.007|
|	31	|	ANTIMICROBIAL PEPTIDES	|	2.08	|	0.007|
|	32	|	ANCHORING FIBRIL FORMATION	|	2.06	|	0.008
|	1	|	HOMOLOGOUS DNA PAIRING AND STRAND EXCHANGE	|	-2.65	 |	0	|
|	2	|	ORGANELLE BIOGENESIS AND MAINTENANCE	|	-2.61	 |	0	|
|	3	|	TRNA AMINOACYLATION	|	-2.59	 |	0	|
|	4	|	DISEASES OF DNA REPAIR	|	-2.58	 |	0	|
|	5	|	TRNA PROCESSING IN THE NUCLEUS	|	-2.57	 |	0	|
|	6	|	MITOTIC SPINDLE CHECKPOINT	|	-2.56	 |	0	|
|	7	|	CELL CYCLE CHECKPOINTS	|	-2.56	 |	0	|
|	8	|	RETROGRADE TRANSPORT AT THE TRANS GOLGI NETWORK	|	-2.53	 |	0	|
|	9	|	RESOLUTION OF SISTER CHROMATID COHESION	|	-2.52	 |	0	|
|	10	|	NS1 MEDIATED EFFECTS ON HOST PATHWAYS	|	-2.52	 |	0	|
|	11	|	M PHASE	|	-2.51	 |	0	|
|	12	|	CILIUM ASSEMBLY	|	-2.50	 |	0	|
|	13	|	TRANSLATION	|	-2.49	 |	0	|
|	14	|	PYRUVATE METABOLISM AND CITRIC ACID TCA CYCLE	|	-2.48	 |	0	|
|	15	|	INTERACTIONS OF REV WITH HOST CELLULAR PROTEINS	|	-2.48	 |	0	|
|	16	|	HDR THROUGH HOMOLOGOUS RECOMBINATION HRR	|	-2.47	 |	0	|
|	17	|	DNA REPAIR	|	-2.47	 |	0	|
|	18	|	MITOTIC PROMETAPHASE	|	-2.46	 |	0	|
|	19	|	HDR THROUGH SINGLE STRAND ANNEALING SSA	|	-2.46	 |	0	|
|	20	|	NUCLEAR IMPORT OF REV PROTEIN	|	-2.44	 |	0	|
|	21	|	MITOCHONDRIAL TRANSLATION	|	-2.44	 |	0	|
|	22	|	TRANSPORT OF MATURE MRNAS DERIVED FROM INTRONLESS TRANSCRIPTS	|	-2.44	 |	0	|
|	23	|	SUMOYLATION OF RNA BINDING PROTEINS	|	-2.44	 |	0	|
|	24	|	SUMOYLATION OF DNA REPLICATION PROTEINS	|	-2.44	 |	0	|
|	25	|	SUMOYLATION OF DNA DAMAGE RESPONSE AND REPAIR PROTEINS	|	-2.43	 |	0	|
|	26	|	NUCLEAR PORE COMPLEX NPC DISASSEMBLY	|	-2.42	 |	0	|
|	27	|	EXPORT OF VIRAL RIBONUCLEOPROTEINS FROM NUCLEUS	|	-2.41	 |	0	|
|	28	|	DNA DOUBLE STRAND BREAK REPAIR	|	-2.41	 |	0	|
|	29	|	INTERACTIONS OF VPR WITH HOST CELLULAR PROTEINS	|	-2.41	 |	0	|
|	30	|	TRANSPORT OF THE SLBP DEPENDANT MATURE MRNA	|	-2.39	 |	0	|
|	31	|	MITOTIC METAPHASE AND ANAPHASE	|	-2.37	 |	0	|
|	32	|	TRNA PROCESSING	|	-2.36	 |	0	|
|	33	|	RESOLUTION OF D LOOP STRUCTURES THROUGH SYNTHESIS DEPENDENT STRAND ANNEALING SDSA	|	-2.36	 |	0	|
|	34	|	REGULATION OF CHOLESTEROL BIOSYNTHESIS BY SREBP SREBF	|	-2.36	 |	0	|
|	35	|	AMINO ACIDS REGULATE MTORC1	|	-2.35	 |	2	|
|	36	|	HOMOLOGY DIRECTED REPAIR	|	-2.35	 |	2	|
|	37	|	NUCLEAR ENVELOPE BREAKDOWN	|	-2.35	 |	2	|
|	38	|	MITOTIC PROPHASE	|	-2.35	 |	2	|
|	39	|	CARGO TRAFFICKING TO THE PERICILIARY MEMBRANE	|	-2.34	 |	2	|
|	40	|	SNRNP ASSEMBLY	|	-2.34	 |	2	|
|	41	|	ANTIVIRAL MECHANISM BY IFN STIMULATED GENES	|	-2.34	 |	2	|
|	42	|	MITOCHONDRIAL BIOGENESIS	|	-2.34	 |	2	|
|	43	|	INTRA GOLGI AND RETROGRADE GOLGI TO ER TRAFFIC	|	-2.33	 |	2	|
|	44	|	SEPARATION OF SISTER CHROMATIDS	|	-2.32	 |	2	|
|	45	|	SUMOYLATION OF CHROMATIN ORGANIZATION PROTEINS	|	-2.32	 |	2	|
|	46	|	SUMOYLATION OF SUMOYLATION PROTEINS	|	-2.29	 |	4	|
|	47	|	DNA DAMAGE BYPASS	|	-2.27	 |	4	|
|	48	|	ER TO GOLGI ANTEROGRADE TRANSPORT	|	-2.26	 |	4	|
|	49	|	THE CITRIC ACID TCA CYCLE AND RESPIRATORY ELECTRON TRANSPORT	|	-2.26	 |	4	|
|	50	|	REGULATION OF TP53 ACTIVITY THROUGH PHOSPHORYLATION	|	-2.25	 |	4	|
|	51	|	PROCESSING OF CAPPED INTRON CONTAINING PRE MRNA	|	-2.25	 |	4	|
|	52	|	ISG15 ANTIVIRAL MECHANISM	|	-2.25	 |	4	|
|	53	|	VIRAL MESSENGER RNA SYNTHESIS	|	-2.25	 |	4	|
|	54	|	SUMOYLATION OF UBIQUITINYLATION PROTEINS	|	-2.25	 |	4	|
|	55	|	SYNTHESIS OF PIPS AT THE EARLY ENDOSOME MEMBRANE	|	-2.25	 |	4	|
|	56	|	POSTMITOTIC NUCLEAR PORE COMPLEX NPC REFORMATION	|	-2.24	 |	4	|
|	57	|	AUTOPHAGY	|	-2.24	 |	4	|
|	58	|	RNA POLYMERASE III TRANSCRIPTION INITIATION FROM TYPE 3 PROMOTER	|	-2.24	 |	4	|
|	59	|	S PHASE	|	-2.24	 |	5	|
|	60	|	CYTOSOLIC TRNA AMINOACYLATION	|	-2.22	 |	6	|
|	61	|	FORMATION OF INCISION COMPLEX IN GG NER	|	-2.22	 |	6	|
|	62	|	ACTIVATION OF GENE EXPRESSION BY SREBF SREBP	|	-2.22	 |	6	|
|	63	|	DUAL INCISION IN GG NER	|	-2.22	 |	6	|
|	64	|	G2 M DNA DAMAGE CHECKPOINT	|	-2.21	 |	6	|
|	65	|	RAB GEFS EXCHANGE GTP FOR GDP ON RABS	|	-2.21	 |	6	|
|	66	|	ANCHORING OF THE BASAL BODY TO THE PLASMA MEMBRANE	|	-2.21	 |	6	|
|	67	|	ACTIVATION OF ATR IN RESPONSE TO REPLICATION STRESS	|	-2.21	 |	6	|
|	68	|	G2 M CHECKPOINTS	|	-2.20	 |	7	|
|	69	|	GOLGI TO ER RETROGRADE TRANSPORT	|	-2.20	 |	7	|
|	70	|	MTOR SIGNALLING	|	-2.20	 |	8	|
|	71	|	GLOBAL GENOME NUCLEOTIDE EXCISION REPAIR GG NER	|	-2.19	 |	9	|
</div> 

[Heatmap SNCG (C)] (https://drive.google.com/file/d/1NEVCIKJ0bPD2fdzNX2kR8koU92Cw0OaS/view?usp=drive_link)  
