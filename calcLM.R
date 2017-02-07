library(forecast)
library(longmemo)
library(fracdiff)
library(ggplot2)
library(reshape2)
library(waveslim)
library(rugarch)
library(fractal)

# Compare Estimators of d for FWN
calcLM<-function(series) {
	seriesName<-deparse(substitute(series))
	gph<-fdGPH(series,bandw.exp=0.8)
	whittle<-WhittleEst(series)
	smoothed<-fdSperio(series, bandw.exp = 0.8, beta = 0.9)
	#fractal package
	fdwhittle<-FDWhittle(series)
	# For Wavelets need to pad series out to power of 2.
	len<-2^(as.integer(log(length(series),2))+1)
	wvlt<-list(par=list(0,0))
	tryCatch(wvlt<-fdp.mle(c(rep(0,len-length(series)),series),"mb8"),error = function(c) cat("Wavelet Estimation failed.\n\n"),warning = function(c) {x<-1})
	mle_fd<-fracdiff(series,nar=0,nma=0)
	#now RUGARCH
	aspec<-arfimaspec(mean.model=list(armaOrder=c(0,0), include.mean=FALSE, arfima=TRUE))
	rugarch.fit<-arfimafit(spec=aspec,data=series,solver="hybrid")
	dummy <- capture.output(series.stat<-stationarity(series))
	series.stat.pvals<-attr(series.stat,"pvals")

	if (series.stat.pvals[1]>0.05) cat("\nTests indicate series is stationary.\n\n") else cat(sprintf("\nWarning: Tests indicate series may not be stationary (p-val %0.4f)\n\n",series.stat.pvals[1]))

	df<-data.frame(	method=c("GPH","Smoothed GPH", "Wavelet MLE", "Whittle", "FDWhittle", "Fracdiff MLE", "RUGARCH"),
					d.est=c(gph$d, smoothed$d, ifelse(unlist(wvlt$par[1])==0,NA,unlist(wvlt$par[1])), 
					       whittle$coefficients[1,1]-0.5, fdwhittle, mle_fd$d, rugarch.fit@fit$robust.matcoef[1,1]),
					se.est=c(gph$sd.as, smoothed$sd.as, NA, whittle$coefficients[1,2], NA, mle_fd$stderror.dpq,
					         rugarch.fit@fit$robust.matcoef[1,2]))
	df$lci<-df$d.est-1.96*df$se.est
	df$uci<-df$d.est+1.96*df$se.est
	df$method<-factor(df$method,levels=c("GPH","Smoothed GPH", "Wavelet MLE", "Whittle", "FDWhittle", "Fracdiff MLE", "RUGARCH"))
	p<-ggplot(df[!is.na(df$d.est),],aes(x=method))+ylim(-1.0,1.0)+
		geom_errorbar(aes(ymin=lci,ymax=uci),width=0.2,colour="darkgray")+
		geom_point(aes(y=d.est),color="black")+
		geom_text(aes(y=d.est,label=sprintf("%0.2f",d.est)),hjust=-0.4,size=3)+
		annotate("rect", xmin=0, xmax=Inf, ymin=-1.0, ymax=0.5, alpha = .2)+
		theme(axis.title.x = element_blank(),plot.title=element_text(color="black"),axis.text.x=element_text(angle=60,hjust=1,colour="black"))+
		labs(y = "Estimate of d")+
		ggtitle(bquote(paste("Comparison of Long Memory Estimates for ",.(seriesName))))
	print(p)
	for (i in 1:nrow(df)) if (!is.na(df$d.est[i])) cat(sprintf("%12s d=%0.4f se=%0.4f\n", df$method[i], df$d.est[i],df$se.est[i]))
	cat(sprintf("FD MLE AIC:  %0.2f\nRUGARCH AIC: %0.2f\n\n",-2.0*mle_fd$log.likelihood, rugarch.fit@fit$LLH*(-2)))
	cat("GPH is fracdiff::fdGPH()\nSmoothedGPH is fracdiff::fdSperio()\nWhittle is longmemo::WhittleEst()\nFDWhittle is fractal::FDWhittle()\nWavelet MLE is waveslim:fdp.mle()\nFracdiff MLE is fracdiff::fracdiff()\nRUGARCH is rugarch::arfimafit()\n\n")
}

data(NileMin) #Get the NileMinima data from the LongMemo package
calcLM(NileMin)
