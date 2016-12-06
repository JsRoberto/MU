#--------------------------------------------------------------------------
#Projeto de Filtros Digitais - Elliptic Filter Design 
#--------------------------------------------------------------------------
.libPaths(
      "C:/Users/JoséRoberto/AppData/Roaming/SPB_16.6/R/win-library/3.2")
library(signal)
library(elliptic)
library(ggplot2)
library(gridExtra)

#Cracterísticas dos filtro elíptico
Amax <- 0.5 #dB
Amin <- 60 #dB
fp <- 30 #Hz
fs <- 60 #Hz

# Cálculo das frequências em rad/s
Wp <- 2*pi*fp
Ws <- 2*pi*fs

# Calculo da frequência de corte normalizada
Wsn <- Ws/Wp # pode ser feito fs/fp também

# Cálculo do modulus k e complementary modulus k'
k <- 1/Wsn
m <- k^2
ml <- 1-m

CEI <- function(m) {
      integrand <- function(theta) {1/sqrt(1-m*sin(theta)^2)}
      integrate(integrand, lower = 0, upper = pi/2)
}

K <- CEI(m)
Kl <- CEI(ml)

eps <<- sqrt(10^{0.1*Amax}-1)

k1 <- eps/(sqrt(10^{0.1*Amin}-1))
kl1 <- sqrt(1-k1^2)

K1 <- CEI(k1^2)
Kl1 <- CEI(kl1^2)

N <- K$value*Kl1$value/(Kl$value*K1$value)
N <- round(N)

arcsc <- function(x, m) {
      integrand <- function(t) {1/sqrt((1+t^2)*(1+(1-m)*t^2))}
      integrate(integrand, lower = 0, upper = x)
}

vo <<- (K$value*arcsc(1/eps, kl1)$value)/(N*K1$value)

# Função que calcula os coeficientes do polinômio que define os zeros da
# função de transferência
Zloc <- function(N, K, m) {
      if (exists("Fs.N")) rm(list = "Fs.N", envir = .GlobalEnv)
      # Cálculo dos zeros da função
      i <- seq(N-1, 0, by = -2)
      wz <- 1/(sqrt(m)*sn(i*K$value/N,m)) # os zeros são também os valores
                                    # negativos dos elementos de wz
      wz <- wz[is.finite(wz)]*1i
      Zeros <- lapply(wz, function(zero) c(1, 0, -zero^2)) #lista que armazenará
                                                      #os fatores (s - si)
      lapply(Zeros, function(x) {
            if (!exists("Fs.N")) Fs.N <<- 1
            Fs.N <<- conv(Fs.N, x)
      })
      Fs.N
}

Ploc <- function(N, K, m) {
      if (exists("Fs.D")) rm(list = "Fs.D", envir = .GlobalEnv)
      # Cálculo dos polos da função
      i <- seq(N-1, 0, by = -2)
      spi <- 1i*sn(K$value*i/N + 1i*vo, m)
      Poles <- lapply(spi,
                      function(pole) {
                            if (Im(pole) != 0) {
                                  list(c(1, -pole), c(1, -Conj(pole)))
                            } else c(1, -pole)
      }) #lista que armazenará
      lapply(Poles, function(x) {
            if (!exists("Fs.D")) Fs.D <<- 1
            if (is.list(x)) {
                  for (idx in 1:2) Fs.D <<- conv(Fs.D, x[[idx]])
            } else {
                  Fs.D <<- conv(Fs.D, x)
            }
      })
      Fs.D
}

Zloc(N,K,m)
Ploc(N,K,m)
Gain <- Fs.D[length(Fs.D)]/Fs.N[length(Fs.N)]
TEST <- freqs(Gain*Fs.N, Fs.D, W = seq(0, 10, length.out = 500))
df.test <- data.frame(w_rads = TEST$W,
                      mag_db= 20*log10(Mod(TEST$H)))
p1 <- ggplot(data = df.test, aes(x = w_rads, y = mag_db)) +
      geom_line() +
      scale_x_continuous(breaks = -1:1,
                         labels = 10^{-1:1}) +
      labs(title = "Frequency Response - Elliptic Analog Low Pass Prototype",
           x = "Frequency normalized (rad/s)", y = "Magnitude (dB)")
p1
ttest <- df.test$mag_db + 60
ttest <- ttest^2
10^{df.test$log10w_rads[ttest == min(ttest)]}

###############PAREI AQUI!
freqs.plot <- function(filter.freqs) {
      filter.df <- data.frame(w = filter.freqs$W,
                              mag = abs(filter.freqs$H),
                              mag_dB = 20*log10(abs(filter.freqs$H)),
                              phase_degrees = (180/pi)*Arg(filter.freqs$H))
      filter.df <- na.omit(filter.df)
      cutoff.freq <- function(filter.df, interval = 1:length(filter.df$mag)) {
            vec_aux <- filter.df$mag_dB[interval]
            vec_aux <- vec_aux + 3
            vec_aux <- vec_aux^2
            wc_sample <- (1:length(vec_aux))[vec_aux==min(vec_aux)]
            wc_sample
      }
      if (round(filter.df$mag[length(filter.df$mag)]) == 1
          & round(filter.df$mag[1]) == 1) {#rejeita-faixa
            filter.type <- "rejeita-faixa"
            central.freq <- (1:length(filter.df$mag))[min(filter.df$mag)] 
            wc_sample1 <- cutoff.freq(filter.df, 1:central.freq)
            wc_sample2 <- cutoff.freq(filter.df, central.freq:length(filter.df$mag))
            interval <- wc_sample1:wc_sample2
      } else {
            if (round(filter.df$mag[length(filter.df$mag)]) == 1) {#passa-alta
                  filter.type <- "passa-alta"
                  wc_sample <- cutoff.freq(filter.df)
                  interval <- 1:wc_sample
            }
            if (round(df$mag[1]) == 1) {#passa-baixa
                  filter.type <- "passa-baixa"
                  wc_sample <- cutoff.freq(filter.df)
                  interval <- wc_sample:length(filter.df$mag)
            }
            if (round(filter.df$mag[length(filter.df$mag)]) == 0
                & round(filter.df$mag[1]) == 0) {#passa-faixa
                  filter.type <- "passa-faixa"
                  central.freq <- (1:length(filter.df$mag))[max(filter.df$mag)] 
                  wc_sample1 <- cutoff.freq(filter.df, 1:central.freq)
                  wc_sample2 <- cutoff.freq(filter.df, central.freq:length(filter.df$mag))
                  interval <- c(1:wc_sample1,
                                rep(NA, length(filter.df) + wc_sample1 - wc_sample2 + 1),
                                wc_sample2:length(filter.df$mag))
            }
      }
      
      p1 <- ggplot(data = filter.df, aes(x = w, y = mag)) +
            geom_line(color = "blue3", size = 1) +
            geom_line(data = filter.df[interval,], color = "red", size = 1.2) +
            labs(x = "", y = "Magnitude")
      p2 <- ggplot(data = filter.df, aes(x = w, y = mag_dB)) +
            geom_line(color = "blue3", size = 1) +
            geom_line(data = filter.df[interval,], color = "red", size = 1.2) +
            labs(x = "", y = "Magnitude (dB)")
      p3 <- ggplot(data = filter.df, aes(x = w, y = phase_degrees)) +
            geom_line(color = "blue3", size = 1) +
            labs(x = "", y = "Fase (graus)")
      
      labs.title <- labs(title = paste("Filtro", filter.type,
                                       "- Resposta em Frequência"))
      labs.x <- labs(x = expression(paste("Frequência Normalizada (x ", pi,
                                          " rad/amostra)")))
      
      grid.arrange(p1 + labs.title, p2, p3 + labs.x, nrow = 3)
}

