library(dplyr)

brainfinal <- (data.frame(read.csv(file = "../outputs/brain_ppi_mapping.csv")))[,3:5]
immunofinal <- (data.frame(read.csv(file = "../outputs/immuno_ppi_mapping.csv")))[,3:5]

colnames(brainfinal)[3] <- "braincounts"
colnames(immunofinal)[3] <- "immunocounts"

brainfinal$comb <- paste(brainfinal$GeneA, brainfinal$GeneB)
immunofinal$comb <- paste(immunofinal$GeneA, immunofinal$GeneB)

merge <- inner_join(x = brainfinal, y = immunofinal[3:4], by = "comb")
mapping <- distinct(.data = merge)
mapping <- mapping[,c(4,1,2,3,5)]

library(ggplot2)
library(rbokeh)
mapping$logbrain <- log2(mapping$braincounts + 1)
mapping$logimmuno <- log2(mapping$immunocounts + 1)

logimmunocutoff <- log2(immunocutoff+1)
logbraincutoff <- log2(braincutoff+1)

ggplot(mapping, aes(y=braincounts, x=immunocounts)) + geom_point()
ggplot(mapping, aes(y=log2(braincounts+1), x=log2(immunocounts+1))) + geom_point() +
  geom_vline(xintercept = log2(immunocutoff+1),color="red", linetype="dashed")+
  geom_hline(yintercept = log2(braincutoff+1),color="red", linetype="dashed")

#Figure for raw counts
figure(height = 750, width = 1500, xlab="Log-normlized counts of immune cells", ylab="Log-normlized counts of brain") %>%
  ly_points(x=immunocounts, y=braincounts, data = mapping,
            hover = list(comb ,braincounts, immunocounts),
            legend = FALSE, alpha = 0.6)

#Figure for Log-normalized counts
figure(height = 750, width = 1500, xlab="Log-normlized counts of immune cells", ylab="Log-normlized counts of brain") %>%
    ly_points(x=logimmuno, y=logbrain, data = mapping ,
            hover = list(comb ,logbrain, logimmuno), 
            legend = FALSE, alpha = 0.6) %>%
  ly_abline(v = logimmunocutoff) %>%
  ly_abline(h =logbraincutoff)

mapping_cutoff <- mapping %>%
  filter(logimmuno >= logimmunocutoff & logbrain >= logbraincutoff)