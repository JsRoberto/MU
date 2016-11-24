#--------------------------------------------------------------------------
#Projeto Time Series - Análise pelo método dos mínimos quadrados
#--------------------------------------------------------------------------

#Definindo a biblioteca e os pacotes não padrões a serem utilizados.
.libPaths(
      "C:/Users/JoséRoberto/AppData/Roaming/SPB_16.6/R/win-library/3.2")
library(ggplot2)

#Caso necessário, o arquivo "dataMU.csv" é baixado.
urlFile <-
      c("https://raw.githubusercontent.com/JsRoberto/MU/master/dataMU.csv")
localFile <- c("./dataMU2.csv")

if (!exists(localFile)) {
      download.file(urlFile, localFile)
}

#A série temporal obtida é armazenada na variável "Z"
dataMU <- read.csv("./dataMU.csv", header = FALSE)
Z <<- ts(dataMU, start = c(1948,1), frequency = 12)

#A função "MMQz()" tem o objetivo de gerar o data.frame "dataMU", que apre-
#senta dados de modelagem da série "Z" obtidos mediante o método dos míni-
#mos quadrados. Entre os quais, temos: (1) série "Zest" que se aproxima de
#"Z", (2) tendência polinomial "Tt" de ordem "pol.order", (3) sazonalidade
#determinística "St", (4) resíduo "at" entre "Zest" e "Z". 
#Além disso, outros parâmetros também são gerados, ligados diretamente aos
#"anos.prev", que indica quantos anos finais da série a modelagem tentará
#prever.
MMQz <- function(Z, pol.order, anos.prev) {
      #Abaixo, são calculados os parâmetros que irão indicar qual subcon-
      #junto de "Z" será utilizado na modelagem e qual seré utilizdo para
      #previsão.
      idx1 <<- 1
      idx2 <<- length(Z)-12*anos.prev
      idx3 <<- length(Z)
      ap <<- anos.prev
      Zoriginal1 <- Z[idx1:idx2]
      Zoriginal2 <- Z[(idx2+1):idx3]
      N <- length(Zoriginal1) 
      anos <- N/12
      
      #Obtenção das matrizes "mT1" e "mT2".
      matrixT <- function(inc, fim, pol.order) {
            expp <<- NULL
            #A função "lst2vct()" transforma uma lista de vetores em um
            #único vetor
            lst2vct <- function(lst) {
                  vct <- vector()
                  for (k in 1:length(lst)) {
                        vct <- c(vct, lst[[k]])
                  }
                  vct
            }
            listT <- tapply(rep(inc:fim, pol.order+1),
                            gl(pol.order+1, fim-inc+1),
                            function(x) {
                                  if (is.null(expp)) {
                                        expp <<- 0
                                  } else expp <<- expp + 1
                                  x^expp
                            })
            vectT <- lst2vct(listT)
            mT <- matrix(vectT, fim-inc+1, pol.order+1)
            mT
      }
      mT1 <- matrixT(idx1, idx2, pol.order)
      mT2 <- matrixT(idx2+1, idx3, pol.order)
      
      #Obtenção das matrizes "mD1" e "mD2".
      mD1 <- mD2 <- mDaux <- rbind(diag(rep(1,11)), rep(-1,11))
      for(i in 1:(anos-1)) mD1 <- rbind(mD1,mDaux)
      for(i in 1:(anos.prev-1)) mD2 <- rbind(mD2,mDaux)
      
      #Obtenção da matriz "X" e da estimador de mínimos quadrados "gamma".
      X <- cbind(mT1, mD1)
      gamma <- chol2inv(chol(t(X)%*%X))%*%t(X)%*%Zoriginal1
      beta <- gamma[1:(pol.order+1)]
      alpha <- gamma[-(1:(pol.order+1))]
      
      #Os dados que serão armazenados em "dataMU" possuem tamanhos diferen-
      #tes, por isso precisam ser preenchidos por missing values ou NAs.
      Zprevisto <- mT2%*%beta + mD2%*%alpha
      Zprevisto <- c(rep(NA, N), Zprevisto)
      Zoriginal1 <- c(Zoriginal1, rep(NA, anos.prev*12))
      Zoriginal2 <- c(rep(NA, N), Zoriginal2)

      #Obtenção da série estimada "Zest", da tendência polinomial "Tt", da
      #sazonalidade determinística "St" e do resíduo "at".
      Zestimado <- c(X%*%gamma, rep(NA, anos.prev*12))
      Tt <- c(mT1%*%beta, rep(NA, anos.prev*12))
      St <- c(mD1%*%alpha, rep(NA, anos.prev*12))
      at <- Zoriginal1 - Zestimado
      
      #Finalmente, os dados são armazenados no data frame "dataMU". 
      dataMU <<- data.frame(Zorg1 = Zoriginal1, Zorg2 = Zoriginal2,
                            Zest = Zestimado, Zprev = Zprevisto,
                            Tt = Tt, St = St, at = at, t = 1:length(Z))
}

#A função "plotz()" produz quatro objetos gráficos, que plotam as caracte-
#rísticas do data.frame "dataMU".
plotz <- function() {
      #A variável abaixo modifica configuração de protagem padrão da bi-
      #blioteca "ggplot2".
      mytheme <- theme(plot.title=element_text(face="bold", size="14",
                                               color="brown"),
                       axis.title=element_text(face="bold", size=10,
                                               color="brown"),
                       axis.text=element_text(face="bold", size=9,
                                              color="darkblue"),
                       panel.background=element_rect(fill="white",
                                                     color="darkblue"),
                       panel.grid.major.x=element_line(color="grey",
                                                       linetype=1),
                       panel.grid.minor.x=element_line(color="grey",
                                                       linetype=2),
                       panel.grid.minor.y=element_blank())
      
      #O objeto "p1" plota a série original "Zorg1" e a série estimada
      #"Zest".
      p1 <<- ggplot(dataMU[idx1:idx2,], aes(x = t, y = Zorg1)) + 
             geom_line(color = "blue3", size = 1.1) + 
             geom_line(mapping = aes(y = Zest), color = "red3",
                       size = 1.1) +
             scale_x_continuous(breaks = seq(12, idx2, by = 24),
                                labels = seq(start(Z)[1], end(Z)[1]-ap,
                                             by = 2)) +
             labs(title = paste("American male unemployment (16-19 years)",
                                paste0("1948-", end(Z)[1] - ap)),
                  x = "Months (December)", y = "Unemployers (1000s)") +
             mytheme
      
      #O objeto "p2" plota a série original "Zorg1", tendência "Tt" e sazo-
      #nalidade "St".
      p2 <<- ggplot(dataMU[idx1:idx2,], aes(x = t, y = Zorg1)) +
             geom_line(color = "blue3", size = 1.1) + 
             geom_line(mapping = aes(y = St), color = "red3", size = 1.1) +
             geom_line(mapping = aes(y = Tt), color = "green3",
                       size = 1.1) +
             scale_x_continuous(breaks = seq(12, idx2, by = 24),
                                labels = seq(start(Z)[1], end(Z)[1]-ap,
                                             by = 2)) +
             labs(title = paste("American male unemployment (16-19 years)",
                                paste0("1948-", end(Z)[1] - ap)),
                  x = "Months (December)", y = "Unemployers (1000s)") +
             mytheme
      
      #O objeto "p3" plota o resíduo "at" e sua média, sendo obtido pela
      #diferença entre a série original "Zorg1" e a série estimada "Zest".
      p3 <<- ggplot(dataMU[idx1:idx2,], aes(x = t, y = at)) +
             geom_line(color = "blue3", size = 1.1) +
             geom_line(mapping = aes(y = mean(at)),
                       color = "red3", lty = "dashed", size = 1.1) +
             scale_x_continuous(breaks = seq(12, idx2, by = 24),
                                labels = seq(start(Z)[1], end(Z)[1]-ap,
                                             by = 2)) +
             labs(title = paste("American male unemployment (16-19 years)",
                                paste0("1948-", end(Z)[1] - ap)),
                  x = "Months (December)", y = "Residue a(t)") + 
             mytheme
      
      #O objeto "p4" plota a séria original "Zorg2" e a série de previsão
      #"Zprev".
      p4 <<- ggplot(dataMU[(idx2+1):idx3,], aes(x = t, y = Zorg2)) +
             geom_line(color = "blue3", size = 1.1) +
             geom_line(mapping = aes(y = Zprev), color = "green3",
                       size = 1.1) +
             scale_x_continuous(breaks = seq(idx2 + 1, idx3, by = 2),
                                labels = rep(c("Jan","Mar","May","Jul",
                                               "Sep","Nov"), ap)) +
             labs(title = paste("Forcasting U.S. male unemployment",
                                paste0("(16-19 years) ", end(Z)[1]-ap+1,
                                       "-1981")),
                  x = "Months (All)", y = "Unemployers (1000s)") +
             mytheme
}

#A função "EMQ()" aplica o erro quadrático médio ao erro normalizado entre
#a série original "Zorg2" e a série de previsão "Zprev". 
EQM <- function() {
      erro <- dataMU$Zorg2 - dataMU$Zprev
      erro.normalizado <- erro/max(sqrt(erro^2), na.rm = TRUE)
      mse <- mean(erro.normalizado^2, na.rm = TRUE)
      mse
}

##########################################################################
#Aplicação de uma tendência linear (pol.order = 1) e previsão para 2 anos
#(anos.prev = 2).
MMQz(Z, 1, 2)
plotz()
EQM()
p1; p2; p3; p4

#Aplicação de uma tendência quadrática (pol.order = 2) e previsão para 3 
#anos (anos.prev = 3).
MMQz(Z, 2, 3)
plotz()
EQM()
p1; p2; p3; p4

#Aplicação de uma tendência de polinômio do 5º grau (pol.order = 5) e pre-
#visão para 3 anos (anos.prev = 3).
MMQz(Z, 5, 3)
plotz()
EQM()
p1; p2; p3; p4