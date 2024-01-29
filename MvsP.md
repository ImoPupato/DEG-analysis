#### En este apartado nos propusimos realizar un análisis de la expresión de ciertos genes de nuestro interés en personas con Melanoma Cutánemo utilizando la información disponible en el proyecto TCGA-SKCM, diferenciando si provenían de muestras de metástasis o del tumor primario.  
#### Vale aclarar que en todos los casos se trabajó con RNASeq y las personas todas fueron diagnosticadas con melanoma cutáneo como tumor primario.  

#### Para comparar las expresiones (RNASeq) de las familias Vav y Sinucleina en muestras provenientes de metástasis y tumor primario, asignamos a cada muestra su proveniencia, siendo 101 correspondendientes a tumores primarios y 359 a metastásicos. Las expresiones fueron convertidas a CPM (*count per million*) con el paquete edgeR.  

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
|VAV1| 7.94 - 13.27 |17.41 | 6.80 |4.339 x $10^{-14}$|
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

#### Se llevaron adelante dos tipos de análisis: ORA y GSEA. 
#### ORA
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

#### GSEA
Se generaron las siguientes matrices:
- [Matriz de expresión análisis Completo o MvsP](https://drive.google.com/file/d/1TE0wsNRsJk0JiLgzSH_G1ka6et60kKB_/view?usp=drive_link)
- [Matriz de expresión análisis expresión de genes de interes en muestras provenientes de tejido tumoral metastásico](https://drive.google.com/file/d/1_ZK_z0lGFaCkXZgYe6DH6b2WYyKztE6t/view?usp=drive_link)
- [Matriz de expresión análisis expresión de genes de interes en muestras provenientes de tejido tumoral primario](https://drive.google.com/file/d/1_H80pHoUc4c9PdDOTFMcMm9LGmWRuTSx/view?usp=drive_link)

Los fenotipos de contraste son:
- [Metástasis vs Primario](https://drive.google.com/file/d/1TE0wsNRsJk0JiLgzSH_G1ka6et60kKB_/view?usp=drive_link)
- Expresión de Vav1: [Completo](https://drive.google.com/file/d/1_igTwl-hgYfYzTiNlcIN8p61Q8QhQTfb/view?usp=drive_link), [Metastasis](https://drive.google.com/file/d/1_aT2Mj34gWdk4tJut97PoR_OxVgspjVb/view?usp=drive_link), [Primario](https://drive.google.com/file/d/1_L1zmKlKLINVPx1EX9WGgWd9SrbAryE9/view?usp=drive_link)
- Expresión de Vav2: [Completo](https://drive.google.com/file/d/1_jUzpTvgPC3ESDEGV5UL9PONVpAINtHj/view?usp=drive_link), [Metastasis](https://drive.google.com/file/d/1_c7ZFpLVxqsdr-JKmUoc9DOObfUDCu5B/view?usp=drive_link), [Primario](https://drive.google.com/file/d/1_X2TNGNE109nRLOFoTqRwTXDqmmKguGo/view?usp=drive_link)
- Expresión de Vav3: [Completo](https://drive.google.com/file/d/1_mRCUAPjYrvzWccvJwwD5NcpUy87jFYk/view?usp=drive_link), [Metastasis](https://drive.google.com/file/d/1_daAS897b7kYBZXkUn19ZaovEgxsrDNL/view?usp=drive_link), [Primario](https://drive.google.com/file/d/1_LGZdctT86jekzC5B--kXM5zMr1qKZ1g/view?usp=drive_link)
- Expresión de SNCA: [Completo](https://drive.google.com/file/d/1_oCaT0Bgr3j06DU3LTppYxCFuMscZH0K/view?usp=drive_link), [Metastasis](https://drive.google.com/file/d/1_dzjdtkNzaw7lMbOY4kSsadKhFKiMCyE/view?usp=drive_link), [Primario](https://drive.google.com/file/d/1_LneeL2rS02Uuz1MEQu-RYlIM-Ht_N6z/view?usp=drive_link)
- Expresión de SNCB: [Completo](https://drive.google.com/file/d/1_qlTwn2nLClX35w3e_v-xKND5V6-sKvz/view?usp=drive_link), [Metastasis](https://drive.google.com/file/d/1_endX2CGNLx0A4nKGFP_eWP9Op4H243s/view?usp=drive_link), [Primario](https://drive.google.com/file/d/1_UAkkxjCGXzOxitAagC12vPVeNrwKPGv/view?usp=drive_link)
- Expresión de SNCG: [Completo](https://drive.google.com/file/d/1_qsL2sMC2HDmmERVtuaeC29tMoOxDu_B/view?usp=drive_link), [Metastasis](https://drive.google.com/file/d/1_ezJTaxVUrLI-UcsrhxcqXNVEF5dgHI3/view?usp=drive_link), [Primario](https://drive.google.com/file/d/1_WFZ50xb8bYjEMWG_xPfoSjmIGn0S5z5/view?usp=drive_link)
  
Las salidas de GSEA de este análisis están disponibles en:
- [Metástasis vs Primario]()
- Expresión de Vav1: [Completo](https://drive.google.com/file/d/1_xmQSazX3In7gTSk-GiiCSQovgSzbeX6/view?usp=sharing), [Metastasis](https://drive.google.com/file/d/192uFvj194FxrqMJVtehVyiK9BKR9Fl0h/view?usp=sharing), [Primario](https://drive.google.com/file/d/196rG9jW5C2dLNsPQaKZ0lHXTw0Ryqfte/view?usp=sharing)
- Expresión de Vav2: [Completo](https://drive.google.com/file/d/1a-Z7Q9TbFefFoFn-WMkO1Dtub_-B_1BN/view?usp=sharing), [Metastasis](https://drive.google.com/file/d/193WkZ8RZMG_E06MKXOIji40OYFCZJz8O/view?usp=sharing), [Primario](https://drive.google.com/file/d/197N6hMANlz0bSG6SmsfFwvDqRhvZ0tQO/view?usp=drive_link)
- Expresión de Vav3: [Completo](https://drive.google.com/file/d/1_ymrPtvOjzhfz_u5hYsAL0mM9fB6RF7z/view?usp=drive_link), [Metastasis](https://drive.google.com/file/d/193Xno3LGx83ax70YP2L_j3S-px7pvqIR/view?usp=sharing), [Primario](https://drive.google.com/file/d/1940Y27GFJ8c0nYGSy90ceW-cANNd213F/view?usp=drive_link)
- Expresión de SNCA: [Completo](https://drive.google.com/file/d/1a5Zd49GiC-7Tmx6PbFLWr7wJJyCxXvMC/view?usp=sharing), [Metastasis](https://drive.google.com/file/d/191_a_KLDzVxIuZo_U4LEviVpC90601Fu/view?usp=sharing), [Primario](https://drive.google.com/file/d/19631KSQgKSQWgcFkfZjcMHY_sucblhEU/view?usp=drive_link)
- Expresión de SNCB: [Completo](https://drive.google.com/file/d/1a022aZrwjb9qjcsStb7-v91JzaaYKsuw/view?usp=sharing), [Metastasis](https://drive.google.com/file/d/191MqzDOEF_72bEun1-A8HTcl6V5Glc9G/view?usp=sharing), [Primario](https://drive.google.com/file/d/1967KRZD2NwVzci0as58eiMPe15RtlUqI/view?usp=drive_link)
- Expresión de SNCG: [Completo](https://drive.google.com/file/d/1a62NceKzJ3jay8NCknTzt8Q0RpfkQiua/view?usp=sharing), [Metastasis](https://drive.google.com/file/d/190D_S5d3R1T0R7XMHVci7CnXQVhB5QnI/view?usp=sharing), [Primario](https://drive.google.com/file/d/196T3oDDSPgtTzFV1zwC4EcIPh26j0sQR/view?usp=drive_link)  

A continuación se muestra el resumen de las vías significaticas de cada contraste:  
#### Vav1 Completo:
<div align="center">
	
|Orden |Via Reactome |NOM p-val |FWER p-val |
|:-----------:|:--------:|:---------:|:---------:|								
|	1	|	IMMUNOREGULATORY INTERACTIONS BETWEEN A LYMPHOID AND A NON LYMPHOID CELL	|	0	|	0	|
|	2	|	COMPLEMENT CASCADE	|	0	|	0	|
|	3	|	INITIAL TRIGGERING OF COMPLEMENT	|	0	|	0	|
|	4	|	CREATION OF C4 AND C2 ACTIVATORS	|	0	|	0	|
|	5	|	FCGR ACTIVATION	|	0	|	0	|
|	6	|	INTERFERON GAMMA SIGNALING	|	0	|	0	|
|	7	|	ANTIGEN ACTIVATES B CELL RECEPTOR BCR LEADING TO GENERATION OF SECOND MESSENGERS	|	0	|	0	|
|	8	|	FCGR3A MEDIATED IL10 SYNTHESIS	|	0	|	0	|
|	9	|	CHEMOKINE RECEPTORS BIND CHEMOKINES	|	0	|	0	|
|	10	|	PARASITE INFECTION	|	0	|	0	|
|	11	|	FCERI MEDIATED CA 2 MOBILIZATION	|	0	|	0	|
|	12	|	SCAVENGING OF HEME FROM PLASMA	|	0	|	0	|
|	13	|	FCERI MEDIATED MAPK ACTIVATION	|	0	|	0	|
|	14	|	CD22 MEDIATED BCR REGULATION	|	0	|	0	|
|	15	|	BINDING AND UPTAKE OF LIGANDS BY SCAVENGER RECEPTORS	|	0	|	0	|
|	16	|	CELL SURFACE INTERACTIONS AT THE VASCULAR WALL	|	0	|	0	|
|	17	|	ANTI INFLAMMATORY RESPONSE FAVOURING LEISHMANIA PARASITE INFECTION	|	0	|	0	|
|	18	|	INTERFERON ALPHA BETA SIGNALING	|	0	|	0	|
|	19	|	GENERATION OF SECOND MESSENGER MOLECULES	|	0	|	0	|
|	20	|	ROLE OF PHOSPHOLIPIDS IN PHAGOCYTOSIS	|	0	|	0	|
|	21	|	ROLE OF LAT2 NTAL LAB ON CALCIUM MOBILIZATION	|	0	|	0	|
|	22	|	FCGAMMA RECEPTOR FCGR DEPENDENT PHAGOCYTOSIS	|	0	|	0	|
|	23	|	SIGNALING BY THE B CELL RECEPTOR BCR	|	0	|	0	|
|	24	|	TCR SIGNALING	|	0	|	0	|
|	25	|	LEISHMANIA INFECTION	|	0	|	0	|
|	26	|	PD 1 SIGNALING	|	0	|	0	|
|	27	|	INTERLEUKIN 10 SIGNALING	|	0	|	0	|
|	28	|	INTERLEUKIN 2 FAMILY SIGNALING	|	0	|	0	|
|	29	|	COSTIMULATION BY THE CD28 FAMILY	|	0	|	0	|
|	30	|	SIGNALING BY INTERLEUKINS	|	0	|	0	|
|	31	|	POTENTIAL THERAPEUTICS FOR SARS	|	0	|	0	|
|	32	|	CLASS A 1 RHODOPSIN LIKE RECEPTORS	|	0	|	0	|
|	33	|	INTERLEUKIN 4 AND INTERLEUKIN 13 SIGNALING	|	0	|	0	|
|	34	|	ANTIGEN PROCESSING CROSS PRESENTATION	|	0	|	0	|
|	35	|	FCERI MEDIATED NF KB ACTIVATION	|	0	|	0	|
|	36	|	FC EPSILON RECEPTOR FCERI SIGNALING	|	0	|	0	|
|	37	|	PEPTIDE LIGAND BINDING RECEPTORS	|	0	|	0	|
|	38	|	DAP12 INTERACTIONS	|	0	|	0	|
|	39	|	INTEGRIN CELL SURFACE INTERACTIONS	|	0	|	0	|
|	40	|	NEUTROPHIL DEGRANULATION	|	0	|	0	|
|	41	|	INTERLEUKIN 3 INTERLEUKIN 5 AND GM CSF SIGNALING	|	0	|	0	|
|	42	|	INTERFERON SIGNALING	|	0	|	0	|
|	43	|	TOLL LIKE RECEPTOR CASCADES	|	0	|	0	|
|	44	|	TNFR2 NON CANONICAL NF KB PATHWAY	|	0	|	0	|
|	45	|	GPCR LIGAND BINDING	|	0	|	0	|
|	46	|	GPVI MEDIATED ACTIVATION CASCADE	|	0	|	0	|
|	47	|	INTERLEUKIN 12 FAMILY SIGNALING	|	0	|	0	|
|	48	|	THE ROLE OF NEF IN HIV 1 REPLICATION AND DISEASE PATHOGENESIS	|	0	|	0.001	|
|	49	|	TNFS BIND THEIR PHYSIOLOGICAL RECEPTORS	|	0	|	0.003	|
|	50	|	DAP12 SIGNALING	|	0	|	0.003	|
|	51	|	INTERLEUKIN RECEPTOR SHC SIGNALING	|	0	|	0.003	|
|	52	|	SIGNAL REGULATORY PROTEIN FAMILY INTERACTIONS	|	0	|	0.004	|
|	53	|	DECTIN 2 FAMILY	|	0	|	0.004	|
|	54	|	G ALPHA I SIGNALLING EVENTS	|	0	|	0.005	|
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
	
|Orden |Via Reactome |NOM p-val |FWER p-val |
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
