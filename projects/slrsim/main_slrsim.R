source("collect_slrsim.R")
all.data = get.all.data()
res = summarize.results(all.data)
write.csv(res,file="slrsim_out.csv",row.names=F)