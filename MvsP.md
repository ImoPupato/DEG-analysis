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
