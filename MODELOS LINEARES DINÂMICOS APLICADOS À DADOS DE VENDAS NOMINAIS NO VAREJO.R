## Bibliotecas
library(ggplot2)#Gráficos
library(dlm)#Modelo e previsão
library(lubridate)#Datas dos dados

# Leitura dos dados
dados <- read.csv("dadosmercado.csv", sep = ";", dec = ",", as.is = TRUE)
dados <- dados[, -3]; names(dados) <- c("tempo", "venda")
dados$tempo <- as.Date(paste(dados$tempo, ".01",sep=""), "%Y.%m.%d")
dados2 <- dados
obs <- dados[dados$tempo > '2019-09-01',] #Últimas 13 observações
dados <- dados[dados$tempo <= '2019-09-01',]#Dados sem as últimas 13 observações

# Gráfico da série temporal
ggplot(dados2, aes(x = tempo, y = venda)) +
  geom_line(color = '#251e3e') +
  ggtitle("Vendas Nominais - Varejo de Hipermercados e Supermercados") +
  labs(x = "Mês/Ano", y = "") +
  scale_x_date(NULL, date_labels = "%m/%y", date_breaks = "36 month") + 
  theme(plot.title = element_text(hjust = 0.5))

#Gráfico de Sazonalidade
avarege <- ts(dados2$venda, start = 2000, frequency=12)
ggseasonplot(avarege, col=rainbow(12), year.labels=TRUE, year.labels.left=TRUE)+
  ylab("Salário Mínimo Real") +
  ggtitle("Gráfico de sazonalidade")+
  theme(plot.title = element_text(hjust = 0.5))

#Gráfico Correlograma
ggAcf(dados2$venda)+
  theme(legend.position = "none") +
  ggtitle("Correlograma") +
  labs(x="Defasagem", y="Aurocorrelação")+
  theme(plot.title = element_text(hjust = 0.5))

# Especificação do modelo DLM
model <- function(p){
  return(
    dlmModPoly(2, dV = p[1], dW = p[2:3]) + #Ordem do modelo (2ª ordem).
      dlmModSeas(12, dV = p[4]) #Parte sazonal.
  )
}
(mle <- dlmMLE(dados$venda,parm = c(0.1, 0.001, 1, 1),
               build = model) #Estimação dos parâmetros via Máxima Verossimilhança)
  
  #Ajuste do modelo
  modelfit = model(mle$par)
  
  # kalman filter
  modelfilter <- dlmFilter(dados$venda, modelfit)
  
  # kalman smoothed
  modelsmoothed <- dlmSmooth(dados$venda, modelfit)
  
  # Número de Previsões passos a frente.
  n <- 13
  fore <- dlmForecast(modelfilter, nAhead = n,sampleNew=100) #Realizando a previsão
  x <- dados$tempo
  xf <- seq(max(x) + 1, ymd(as.Date(max(x))) %m+% months(n), "month")
  df <- rbind(
    data.frame(x=x , y=as.numeric(dados$venda), series = 'Original'),
    data.frame(x=x , y=apply(modelfilter$m[-1,1:2],1,sum), series = 'Filtrada'),
    data.frame(x=x , y=apply(modelsmoothed$s[-1,1:2],1,sum), series = 'Suavisada'),
    data.frame(x=xf , y=fore$f, series = 'Previsão'))
  
  #Gráficos das séries
  (dlm <- ggplot(df, aes(x = x, y=y)) + 
      geom_line(aes(col = series),size = 0.8) +
      scale_x_date(NULL, date_labels = "%m/%y", date_breaks = "36 month") +
      xlab("Valor") + ylab("Vendas") + ggtitle('Séries Temporais')
    
    #Limites de predição.
    LI_prev <- (outer(sapply(fore$Q, FUN = function(x) sqrt(diag(x))), qnorm(0.025, lower = FALSE)) +as.vector(t(fore$f)))
    LS_prev <- (outer(sapply(fore$Q, FUN = function(x) sqrt(diag(x))), qnorm(0.975, lower = FALSE)) +as.vector(t(fore$f)))
    
    #Gráfico de Intervalos de Predição
    (ggplot()+
        geom_line(aes(x=1:13, y=LI_prev), color="red",size = 1) +
        geom_point(aes(x=1:13, y=LI_prev, colour="blue")) +
        geom_line(aes(x=1:13, y=LS_prev, colour="blue"),size = 1) +
        geom_point(aes(x=1:13, y=LS_prev, colour="blue")) +
        geom_line(aes(x=1:13, y=obs[,2],  colour="red"),size = 1) +
        geom_point(aes(x=1:13, y=obs[,2], colour="red")) +
        geom_line(aes(x=1:13, y=fore$f, colour="green"),size = 1) +
        geom_point(aes(x=1:13, y=fore$f, colour="green"))+
        ggtitle("Intervalos de predição") +
        scale_color_discrete(name="Valores", labels=c("Limites do intervalo", "Preditos", "Observados"))+
        ylab("Venda")+xlab("Número de previsões") + theme_classic() 
      
      # Avaliação de parâmetros populacionais de interesse e o tempo que demorou para rodar
      ini <- Sys.time()
      set.seed(2020)
      delta = replicate(1000, max(abs(diff(dlmBSample(modelfilter)))))
      zeta = replicate(1000, min(abs(diff(dlmBSample(modelfilter)))))
      fim <- Sys.time()