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

ggplot(mapping, aes(y=braincounts, x=immunocounts)) + geom_point()
ggplot(mapping, aes(y=log2(braincounts+1), x=log2(immunocounts+1))) + geom_point()

# p <- figure() %>%
#   ly_points(mapping$braincounts, mapping$mmunocount,
#             hover = list(mapping$comb))

figure(height = 750, width = 1500) %>%
    ly_points(x=logimmuno, y=logbrain, data = mapping , color = mapping$GeneA,
            hover = list(comb ,logbrain, logimmuno), size = 4, legend = FALSE, alpha = 0.6)