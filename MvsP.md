Nos propusimos realizar un análisis de la expresión de ciertos genes de nuestro interés en personas con Melanoma Cutánemo. Utilizando la información disponible en el proyecto TCGA-SKCM, realizamos los pasos que a continuación se detallan para realizar las comparaciones considerando separa a las muestras según sean provenientes
Ejemplo del código utilizado para el análisis de metástasis o tumos primario. Vale aclarar que en todos los casos se trabajó con RNASeq y las personas todas fueron diagnosticadas con melanoma cutáneo como tumor primario.  
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
Tabla de resultados

<div align="center">
|Gen| IC 95% | $\overline{x}$ metastasis | $\overline{x}$ primario | p-value |
|:---------------------------------:|:---------------------------------:|:---------------------------------:|:---------------------------------:|:---------------------------------:|
|VAV1| 7.94 - 13.27 |17.41 | 6.80 |4.339x10^{-14}|
|VAV2| 2.49 - 10.84 | 33.15 | 26.49 |0.0019 |
|VAV3| 4.66 - 18.64 | 27.59 | 15.94 | 0.0012 |
|SNCA| 4.66 - 18.64 | 27.59 | 15.94 | 0.0012 |
|SNCB| 4.66 - 18.64 | 27.59 | 15.94 | 0.0012 |
|SNCG| 4.66 - 18.64 | 27.59 | 15.94 | 0.0012 |
</div>  
