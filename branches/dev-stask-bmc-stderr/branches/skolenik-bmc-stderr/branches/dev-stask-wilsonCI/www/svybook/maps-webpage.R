library(maptools)
library(ggplot2)
library(survey)

states<-readShapePoly("GIS/brfss_state_2007_download")
cities<-readShapePoints("GIS/BRFSS_MMSA_2007")


nstate<-length(states@polygons)
statecoords<-lapply(1:nstate,
	 function(j){	
	 	state<-states@polygons[[j]]
	 	m<-length(state@Polygons)
	 	do.call(rbind,lapply(1:m, function(i){
	 	    	data.frame(state@Polygons[[i]]@coords,poly=i,fips=states$ST_FIPS[j],st=states$STATE[j])
	 	    	}))	 	
	 	}
	 )

statedata<-do.call(rbind,statecoords) 
names(statedata)[1:2]<-c("x","y")

citydata<-cbind(cities@data,cities@coords)
names(citydata)[74:75]<-c("x","y")
citydata<-citydata[-6,]  ## no data for Columbus, OH


statefull<-merge(statedata,states@data,by.x="fips",by.y="ST_FIPS")

 qplot(x,y,data=statefull,geom="polygon",group=fips*100+poly,fill=X4416_1 ,colour=I("black"))+geom_point(aes(x,y,fill=NULL,colour=X4416_1,group=NULL),data=citydata)
 

colscale<-scale_fill_brewer(name="Insurance",palette="YlGnBu",limits=c(70,100),levels=7) 

qplot(x,y,data=statefull,geom="polygon",group=fips*100+poly,fill=X4416_1 ,colour=I("grey50"))+geom_point(aes(x,y,fill=X4416_1,group=NULL),shape=21,size=3,data=citydata)


nobreaks <- list(
  scale_x_continuous(breaks = NA),
  scale_y_continuous(breaks = NA),
  opts(axis.ticks.length = unit(0, "cm"))
)


qplot(x,y,data=statefull,geom="polygon",group=fips*100+poly,fill=X868_1 ,colour=I("white"))+geom_point(aes(x,y,fill=X868_1,group=NULL),shape=21,size=3,data=citydata)+theme_bw()+scale_fill_gradient(name="Insurance",low="black",high="grey80")+nobreaks

qplot(x,y,data=statefull,geom="polygon",group=fips*100+poly,fill=cut(X868_1,c(0,82,85,87,89,100),right=FALSE) ,colour=I("black"))+geom_point(aes(x,y,fill=cut(X868_1,c(0,82,85,87,89,100),right=FALSE),group=NULL),shape=21,color=I("grey90"),data=citydata)+scale_fill_grey(name="Insurance")+nobreaks+theme_bw()
