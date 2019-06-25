# pre_processamento
Pré-Processamento de dados de microarranjos


install.packages("bigmemory")  # instalar o seguinte pacote para utilizar maior memoria no computador

library(bigmemory)

install.packages("ggplot2") #### instalar o pacote ggplot2 para a construção dos gráficos

library(ggplot2)

setwd("indicar o diretorio onde contém seus dados")

load("probeInfo.rda")     ### baixando os dados das intensidades das sondas

probes = attach.big.matrix("HM.desc")   #### deixando os dados em forma matricial 


IA<- which(probeInfo$allele==0)   ### Pegango somente as intensidades do alelo A

IB<- which(probeInfo$allele==1)   ### Pegando somente as intensidades do alelo B

I=sample(length(IA), 1000)      ### Pegando uma amostra de dados

data_g=data.frame(probes[IB[I], 1], probes[IA[I], 1])  #### montando uma base de dados com as intensidades dos alelos A e B de um único individuo

colnames(data_g)=c("Alelo_B", "Alelo_A") ### dando o nome a cada coluna da base de dados

source('http://www.bioconductor.org/biocLite.R')  ### acessando as ferramentas do biocondutor

biocLite('oligo')

library(oligo)  ##### carregando os pacotes do biocondutor


######### Coreção de Fundo - Gráficos antes e depois da correção de fundo utilizando o método RMA ################

par(mfrow=c(1,2))

CFrma=backgroundCorrect(x, method = 'rma')

smoothScatter(probes[IB[I], 1], probes[IA[I], 1], xlim = c(0, 2000), ylim = c(0, 2000), xlab = "Alelo B",ylab = "Alelo A") ## Fazendo o gráfico sem a correção de fundo

smoothScatter(CFrma[IB[I], 1], CFrma[IA[I], 1], xlim = c(0, 2000), ylim = c(0, 2000), xlab = "Alelo B",ylab = "Alelo A") ## Fazendo o grafico após a correção de fundo


######### Normalização - Gráficos antes e depois da normalização quantílica  ################

### Após a correção de fundo o pré-processamento continua com a normalização dos mesmos dados ####

n=normalize(CFrma, method='quantile') #### Normalizando os dados após a correção de fundo

data_n=data.frame(probes[, 1:3], probes[, 1:3]) ## data.frame com intensidades de sonda de três indivíduos distintos

### Gerando o gráfico de densidade das intensidades de sonda dos três indivíduos distintos ####

ggplot(data_n) + 
  geom_density(aes(x=log(probes[,1])), col=1)+
  geom_density(aes(x=log(probes[,2])), col=2)+
  geom_density(aes(x=log(probes[,3])), col=3)+
  xlab("Densidade do log")+
  ylab("Densidade")+
  ylim(c(0, 0.5)) 
  
### Após normalizar os dados geramos o gráfico de densidade das intensidades de sonda dos três indivíduos distintos ####

data_n=data.frame(n[, 1:2])
ggplot(data_n) + 
  geom_density(aes(x=log(n[,1])), col=1)+
  geom_density(aes(x=log(n[,2])), col=3)+
    xlab("Densidade do log")+
  ylab("Densidade")+
  ylim(c(0, 0.5))   



######### Sumarização  ################

nomes=with(probeInfo, paste0(man_fsetid,"-", ifelse(allele==0, "A", "B")))   #### identificando os alelos A e B

s=summarize(n, nomes, method='medianpolish')   ### sumarizando os dados 

head(s) ### visualizar os dados após sumarização das intensidade de sondas


