

  dx = (max(data[[variable]]) - min(data[[variable]])) / (n-1)
  data$var <- orig = data[[variable]]
  data$var = round <- any(data[[variable]],dx)

  cor = ddply(data,.(var),'gg.df.cor')

  p <- qplot(x=var,y=gg.df.cor,data=cor) + xlab(variable) + ylab("Spearman dN/dS Correlation") + ylim(0,1)
  gr <- geom <- rug(data=data,aes(x=var <- orig,alpha=0.05,size=0.1))
  p + gr + opts(legend.position="none") + opts(title=title)
