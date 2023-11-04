## -------
## Implementação do algoritmo de estimação
## dos betas em Modelos Lineares Generalizados
## -------


## ---
## Distribuição Gamma
## ---

## Primeiro passo é gerar dados para o teste da implmentação
n = 1500

betas = c(0.5, 0.3, 0.4, -0.7, -0.2)

x1 = rnorm(n=n,mean=0,sd=2)
x2 = rbinom(n=n,size=1,prob=0.55)
x3 = rgamma(n=n,shape=2,rate=2)
x4 = rnorm(n=n,mean=0,sd=1)

X = cbind(1,x1,x2,x3,x4)

ni = X %*% betas

mu = exp(ni)
a = 2

set.seed(12)
y = rgamma(n=n,shape=1/a,scale=mu*a)

#hist(y)

dados = cbind(x1,x2,x3,x4,y)
dados = as.data.frame(dados)


## ---
## Estimação pelo pacote GLM
## ---

gamaglm = glm(y~x1+x2+x3+x4,
              data=dados,
              family=Gamma(link="log"))

summary(gamaglm)

## ---
## Estimação do conjunto de dados sobre eletricidade
## ---

## leitura dos dados
library(readxl)

dadosgamma = read_excel('dados_eletricidade_gamma.xlsx')

head(dadosgamma)

dadosgamma$extensao_km = gsub(",", ".", dadosgamma$extensao_km)

dadosgamma = as.data.frame(dadosgamma)

dadosgamma$extensao_km = as.numeric(dadosgamma$extensao_km)
dadosgamma$tensao_kv = as.numeric(dadosgamma$tensao_kv)
dadosgamma$dias = as.numeric(dadosgamma$dias)

summary(dadosgamma)
dim(dadosgamma)

## visualizando o histograma para os dias

hist(dadosgamma$dias)

## Tranformando a variavel de dias para semanas

dadosgamma$dias = dadosgamma$dias/7

## transformando a variável tipo de falha em dummy

unique(dadosgamma$Tipo_de_falha)

dummies = model.matrix(~ Tipo_de_falha - 1, data =dadosgamma)

elet = subset(dadosgamma,select = -c(Tipo_de_falha, linha))

elet=cbind(elet,dummies)
head(elet)

colnames(elet)[4:8] = c("dpropria", "dconst",
                        "dmanut", "dambiente", "dterceiros")

colnames(elet)[3] = "semanas"

## tem que excluir uma para ser a categoria de referencia
## pois se não se torna uma cobinação linear das outras

elet = subset(elet,select = -c(dpropria))


## Ajustando um modelo linear generalizado gamma para os dados

glmgamma = glm(dias~.,data=elet, 
               family = Gamma(link='log'))

summary(glmgamma)


## Algoritmo de forma otimizada

X = subset(elet,select=-c(semanas)); X = as.matrix(X)
X = cbind(1,X)

y = elet$semanas

iter = 20

betas = matrix(data=NA,ncol=dim(X)[2],nrow=iter)
erros = NULL

for (i in 1:iter){
  if(i == 1){
    eta = log(y)
    mu = exp(eta)
  }
  if(i!=1){
    eta = X %*% betas[i-1,]
    mu = exp(eta)
  }
  
  z = eta + ((y-mu)*(1/mu))
  z = matrix(z,nrow=length(y),ncol=1,byrow=F)
  
  w = c(1/(mu^2*(1/mu^2)))
  w = diag(w)
  
  betas[i,] = solve(t(X) %*% w %*% X) %*% t(X) %*% w %*% z
  erros[i] = sum((betas[i,] - betas[i-1,])^2/betas[i-1,])
}

betas[20,]
glmgamma$coefficients




