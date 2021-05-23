library(readxl)
library(reshape2)
library(dplyr)

a1a <- read_excel('./DONOTSYNC/20210311_overview.xlsx', sheet='OK-2')
a1a$study= 'OK2'

a1b <- read_excel('./DONOTSYNC/20210311_overview.xlsx', sheet='OK-3')
a1b$study= 'OK3'

a1c <- read_excel('./DONOTSYNC/20210311_overview.xlsx', sheet='OK-4')
a1c$study= 'OK4'

a1 <- bind_rows(a1a,a1b,a1c)

a1$shared_carriage <-  NULL
a1$shared_serotype <- NULL
a1$serotype <- NULL
a1$Household <- paste(1:nrow(a1), a1$study,sep='_')
a1$study <- NULL
a1.m <- melt(a1, id.vars = 'Household')
a1.m$parent <- 0
a1.m$parent[substr(a1.m$variable,1,1) == 'p'] <- 1
a1.m$variable <-  sub("_[^_]+$", "", a1.m$variable)

a1.m <- a1.m[!(a1.m$variable %in% c('kCE_CS','pCE_CS')),]

a1.m$st <- substring(a1.m$variable,4)
a1.m <- a1.m[,c('Household','parent','st','value')]

a1.m$st[substr(a1.m$st,1,1)=='6'] <- '6'

a1.c <- acast(a1.m, Household ~ parent ~ st, fun.aggregate = sum, na.rm=T)
a1.c[a1.c>1] <- 1 #if a parent is positive for both OP and saliva, just count 1

prev <- t(apply(a1.c,c(2,3) , mean, na.rm=T))
ColonizedN <- t(apply(a1.c,c(2,3) , sum, na.rm=T))


