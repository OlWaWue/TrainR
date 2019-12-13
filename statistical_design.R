#### Statistische Versuchsplanung

library(OPDOE)

data(cattle)

cattle

size.seq_select.mean(data= cattle, delta=10, P=0.95)

size.anova(model="axb",a=4, b=4, alpha=0.05,beta=0.1, delta=0.2, case="maximin") 

