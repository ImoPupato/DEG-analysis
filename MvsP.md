Nos propusimos realizar un análisis de la expresión de ciertos genes de nuestro interés en personas con Melanoma Cutánemo. Utilizando la información disponible en el proyecto TCGA-SKCM, realizamos los pasos que a continuación se detallan para realizar las comparaciones considerando separa a las muestras según sean provenientes Vale aclarar que en todos los casos se trabajó con RNASeq y las personas todas fueron diagnosticadas con melanoma cutáneo como tumor primario.  
Primero decidimos comparar las expresiones (RNASeq) de las familias Vav y Sinucleina en muestras provenientes de metástasis y tumor primario.  
De las 460 muestras 101 corresponden a tumores primarios y 359 a metastásicos.  

Ejemplo del código utilizado para la comparación
```R
t.test(fenotipos$VAV1.cpm[
  fenotipos$Type =="Metastasis"],
  fenotipos$VAV1.cpm[
    fenotipos$Type=="Primario"]
  ,paired = FALSE)
```
Ejemplo de la salida
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
**Tabla de resultados**  
Los valores están en *cuentas por millón* (cpm)

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

Luego, calculamos la mediana (valor de corte) para generar cada grupo de contraste. 

<div align="center">
	
|Grupo|VAV1| VAV2|VAV3| SNCA|SNCB|SNCG|
|:-----------:|:--------:|:---------:|:-----------:|:----------:|:---------:|:----------:|
|Metástasis | 10.19  | 28.18 | 10.34 | 66.52 | 0.47 | 1.62 | 
|Primario | 4.24 | 22.09 | 7.17 | 92.93 | 0.57 | 1.50 | 

</div>  
Tomando en cuenta los valores de corte, los grupos quedaron de la siguiente manera:
