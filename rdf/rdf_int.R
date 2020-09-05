library(readr)
library(ggplot2)
greigite_3360 <- read_table2("greigite_3360_0.1dr_maxR5",col_names = c('Radius','gofr','intR'))
#View(greigite_3360)

fdn_392 <- read_table2("1fdn_392_0.1dr_maxR5",col_names = c('Radius','gofr','intR'))
#View(greigite_3360)

ggplot() + 
  geom_line(data=greigite_3360,aes(x=Radius,y = intR, colour = "black")) + 
  geom_line(data=fdn_392,aes(x=Radius,y=intR, colour = "red"))

