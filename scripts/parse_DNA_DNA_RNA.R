library(ggplot2)

setwd("/home/projects/risearch/subprojects/CRISPR_off_targets/DNA_DNA_RNA_energy_parameters/")

r_nns = c("AA","AC","AG","AU","AN","A-","CA","CC","CG","CU","CN","C-",
          "GA","GC","GG","GU","GN","G-","UA","UC","UG","UU","UN","U-",
          "NA","NC","NG","NU","NN","N-","-A","-C","-G","-U","-N","--")

d_nns =c("AA","AC","AG","AT","AN","A-","CA","CC","CG","CT","CN","C-",
         "GA","GC","GG","GT","GN","G-","TA","TC","TG","TT","TN","T-",
         "NA","NC","NG","NT","NN","N-","-A","-C","-G","-T","-N","--")


r_nns_code=c()
for(nt1 in c('-','A','G','C','U','N')){
  for(nt2 in c('-','A','G','C','U','N')){
    r_nns_code=c(r_nns_code,paste(nt1,nt2,sep=""))
  }
}
d_nns_code=c()
for(nt1 in c('-','A','G','C','T','N')){
  for(nt2 in c('-','A','G','C','T','N')){
    d_nns_code=c(d_nns_code,paste(nt1,nt2,sep=""))
  }
}

r_matches=c("AU","UA","CG","GC","UG","GU")
d_matches=c("AT","TA","CG","GC")
r_nts = c("A","C","G","U","N","-")
d_nts = c("A","C","G","T","N","-")

ris2_minus=559


#####################################





#####################################
######## RIsearch2 parameters #######
RI_RR_all_data=c(-Inf,NA,NA,NA,NA,NA,NA,-40,-40,-40,-40,-40,NA,-40,-40,-40,-40,-40,NA,-40,-40,-40,-40,-40,NA,-40,-40,-40,-40,-40,NA,-40,-40,-40,-40,-40,
                 NA,0,0,0,105,0,NA,-22,-22,-22,-45,-22,NA,-22,-22,-22,-45,-22,NA,-22,-22,-22,-45,-22,NA,-22,-22,-22,-45,-22,NA,-22,-22,-22,-45,-22,NA,0,
                 0,150,105,0,NA,-22,-22,0,-45,-22,NA,-22,-22,0,-45,-22,NA,-22,-22,0,-45,-22,NA,-22,-22,0,-45,-22,NA,-22,-22,0,-45,-22,NA,0,150,0,0,0,NA,
                 -22,0,-22,-22,-22,NA,-22,0,-22,-22,-22,NA,-22,0,-22,-22,-22,NA,-22,0,-22,-22,-22,NA,-22,0,-22,-22,-22,NA,105,105,0,0,0,NA,-45,-45,-22,
                 -22,-22,NA,-45,-45,-22,-22,-22,NA,-45,-45,-22,-22,-22,NA,-45,-45,-22,-22,-22,NA,-45,-45,-22,-22,-22,NA,0,0,0,0,0,NA,-22,-22,-22,-22,-22,
                 NA,-22,-22,-22,-22,-22,NA,-22,-22,-22,-22,-22,NA,-22,-22,-22,-22,-22,NA,-22,-22,-22,-22,-22,NA,NA,NA,NA,NA,NA,-Inf,-71,-71,-71,-71,-71,
                 -Inf,-71,-71,-71,-71,-71,-Inf,-71,-71,-71,-71,-71,-45,-285,-285,-285,-285,-285,-Inf,-71,-71,-71,-71,-71,-40,-22,-22,-22,-45,-22,-71,-22,
                 -22,-22,-70,-22,-71,-22,-22,-22,30,-22,-71,-22,-22,-22,-70,-22,-285,-230,-150,-230,90,-230,-71,-22,-22,-22,-70,-22,-40,-22,-22,0,-45,-22,
                 -71,-22,-22,0,-70,-22,-71,-22,-22,100,30,-22,-71,-22,-22,0,-70,-22,-285,-130,-110,210,60,-230,-71,-22,-22,0,-70,-22,-40,-22,0,-22,-22,-22,
                 -71,-22,0,-22,-22,-22,-71,-22,100,-22,-22,-22,-71,-22,0,-22,-22,-22,-285,-230,220,-230,-230,-230,-71,-22,0,-22,-22,-22,-40,-45,-45,-22,-22,
                 -22,-71,-70,-70,-22,-22,-22,-71,30,30,-22,-22,-22,-71,-70,-70,-22,-22,-22,-285,110,140,-230,-160,-230,-71,-70,-70,-22,-22,-22,-40,-22,-22,-22,
                 -22,-22,-71,-22,-22,-22,-22,-22,-71,-22,-22,-22,-22,-22,-71,-22,-22,-22,-22,-22,-285,-230,-230,-230,-230,-230,-71,-22,-22,-22,-22,-22,NA,NA,
                 NA,NA,NA,NA,-Inf,-71,-71,-71,-71,-71,-Inf,-71,-71,-71,-71,-71,0,-240,-240,-240,-240,-240,-45,-285,-285,-285,-285,-285,-Inf,-71,-71,-71,-71,
                 -71,-40,-22,-22,-22,-45,-22,-71,-22,-22,-22,10,-22,-71,-22,-22,-22,50,-22,-240,-160,-80,-160,240,-160,-285,-230,-150,-230,130,-230,-71,-22,
                 -22,-22,-70,-22,-40,-22,-22,0,-45,-22,-71,-22,-22,80,10,-22,-71,-22,-22,120,50,-22,-240,-60,-40,330,150,-160,-285,-130,-110,210,50,-230,-71,
                 -22,-22,0,-70,-22,-40,-22,0,-22,-22,-22,-71,-22,80,-22,-22,-22,-71,-22,120,-22,-22,-22,-240,-160,340,-160,-160,-160,-285,-230,250,-230,-230,
                 -230,-71,-22,0,-22,-22,-22,-40,-45,-45,-22,-22,-22,-71,10,10,-22,-22,-22,-71,50,50,-22,-22,-22,-240,220,250,-160,-90,-160,-285,140,-130,-230,
                 -160,-230,-71,-70,-70,-22,-22,-22,-40,-22,-22,-22,-22,-22,-71,-22,-22,-22,-22,-22,-71,-22,-22,-22,-22,-22,-240,-160,-160,-160,-160,-160,-285,
                 -230,-230,-230,-230,-230,-71,-22,-22,-22,-22,-22,NA,NA,NA,NA,NA,NA,-Inf,-71,-71,-71,-71,-71,0,-240,-240,-240,-240,-240,-Inf,-71,-71,-71,-71,-71,
                 -Inf,-71,-71,-71,-71,-71,-Inf,-71,-71,-71,-71,-71,-40,-22,-22,-22,-45,-22,-71,-22,-22,-22,-70,-22,-240,-160,-80,-160,210,-160,-71,-22,-22,-22,
                 -70,-22,-71,-22,-22,-22,-70,-22,-71,-22,-22,-22,-70,-22,-40,-22,-22,0,-45,-22,-71,-22,-22,0,-70,-22,-240,-60,-40,240,140,-160,-71,-22,-22,0,-70,
                 -22,-71,-22,-22,0,-70,-22,-71,-22,-22,0,-70,-22,-40,-22,0,-22,-22,-22,-71,-22,0,-22,-22,-22,-240,-160,330,-160,-160,-160,-71,-22,0,-22,-22,-22,
                 -71,-22,0,-22,-22,-22,-71,-22,0,-22,-22,-22,-40,-45,-45,-22,-22,-22,-71,-70,-70,-22,-22,-22,-240,210,210,-160,-90,-160,-71,-70,-70,-22,-22,-22,
                 -71,-70,-70,-22,-22,-22,-71,-70,-70,-22,-22,-22,-40,-22,-22,-22,-22,-22,-71,-22,-22,-22,-22,-22,-240,-160,-160,-160,-160,-160,-71,-22,-22,-22,
                 -22,-22,-71,-22,-22,-22,-22,-22,-71,-22,-22,-22,-22,-22,NA,NA,NA,NA,NA,NA,-45,-285,-285,-285,-285,-285,-45,-285,-285,-285,-285,-285,-Inf,-71,
                 -71,-71,-71,-71,-Inf,-71,-71,-71,-71,-71,-Inf,-71,-71,-71,-71,-71,-40,-22,-22,-22,-45,-22,-285,-230,-150,-230,130,-230,-285,-230,-150,-230,100,
                 -230,-71,-22,-22,-22,-70,-22,-71,-22,-22,-22,0,-22,-71,-22,-22,-22,-70,-22,-40,-22,-22,0,-45,-22,-285,-130,-110,210,100,-230,-285,-130,-110,140,
                 -30,-230,-71,-22,-22,0,-70,-22,-71,-22,-22,70,0,-22,-71,-22,-22,0,-70,-22,-40,-22,0,-22,-22,-22,-285,-230,240,-230,-230,-230,-285,-230,150,-230,
                 -230,-230,-71,-22,0,-22,-22,-22,-71,-22,70,-22,-22,-22,-71,-22,0,-22,-22,-22,-40,-45,-45,-22,-22,-22,-285,90,130,-230,-160,-230,-285,60,50,-230,
                 -160,-230,-71,-70,-70,-22,-22,-22,-71,0,0,-22,-22,-22,-71,-70,-70,-22,-22,-22,-40,-22,-22,-22,-22,-22,-285,-230,-230,-230,-230,-230,-285,-230,
                 -230,-230,-230,-230,-71,-22,-22,-22,-22,-22,-71,-22,-22,-22,-22,-22,-71,-22,-22,-22,-22,-22,NA,NA,NA,NA,NA,NA,-Inf,-71,-71,-71,-71,-71,-Inf,
                 -71,-71,-71,-71,-71,-Inf,-71,-71,-71,-71,-71,-Inf,-71,-71,-71,-71,-71,-Inf,-71,-71,-71,-71,-71,-40,-22,-22,-22,-45,-22,-71,-22,-22,-22,-70,
                 -22,-71,-22,-22,-22,-70,-22,-71,-22,-22,-22,-70,-22,-71,-22,-22,-22,-70,-22,-71,-22,-22,-22,-70,-22,-40,-22,-22,0,-45,-22,-71,-22,-22,0,-70,
                 -22,-71,-22,-22,0,-70,-22,-71,-22,-22,0,-70,-22,-71,-22,-22,0,-70,-22,-71,-22,-22,0,-70,-22,-40,-22,0,-22,-22,-22,-71,-22,0,-22,-22,-22,-71,
                 -22,0,-22,-22,-22,-71,-22,0,-22,-22,-22,-71,-22,0,-22,-22,-22,-71,-22,0,-22,-22,-22,-40,-45,-45,-22,-22,-22,-71,-70,-70,-22,-22,-22,-71,-70,
                 -70,-22,-22,-22,-71,-70,-70,-22,-22,-22,-71,-70,-70,-22,-22,-22,-71,-70,-70,-22,-22,-22,-40,-22,-22,-22,-22,-22,-71,-22,-22,-22,-22,-22,-71,
                 -22,-22,-22,-22,-22,-71,-22,-22,-22,-22,-22,-71,-22,-22,-22,-22,-22,-71,-22,-22,-22,-22,-22)


RI_RR_all=NULL
for (i in 1:36){
  for (j in 1:36){
    RI_RR_all=rbind(RI_RR_all,data.frame(x=r_nns_code[i],y=r_nns_code[j],num=RI_RR_all_data[((i-1)*36)+j]))
  }
}
RI_RR_all$x=as.character(RI_RR_all$x)
RI_RR_all$y=as.character(RI_RR_all$y)

for (i in 0:length(r_nns)){
  if (i==0) {cat("    ",paste(r_nns,collapse = "   "))}
  for (j in 0:length(r_nns)){
    if (j==0) {cat(r_nns[i])}
    else{
      val = as.character(RI_RR_all$num[RI_RR_all$x==r_nns[i] & RI_RR_all$y==r_nns[j]])
      if (length(val)==0) {val=""}
      else if (is.na(val) | val=="NA") {val="*"}
      else if (val==-Inf) {val="-inf"}
      valc=nchar(val)
      cat(paste(rep(" ",(4-valc)), collapse=""),val)
    }
  }
  cat("\n")
}


RI_RR_internal_loop_data = data.frame(size=c(4:30),eng=c(110,200,200,210,230,240,250,260,270,280,290,290,300,310,310,320,330,330,340,340,350,350,350,360,360,370,370))
for (i in 3:length(RI_RR_internal_loop_data$size)){
  test=lm(RI_RR_internal_loop_data$eng[1:i] ~ RI_RR_internal_loop_data$size[1:i])
  print(ggplot(RI_RR_internal_loop_data,aes(x=size,y=eng))+geom_point()+geom_abline(intercept=round(test$coefficients[1]),slope = round(test$coefficients[2]))+
          ggtitle(paste("Internal Loop",as.character(i)))+ylim(0,700))
  print(paste(as.character(i),as.character(round(test$coefficients[1])),(as.character(round(test$coefficients[2])))))
}
ggplot(RI_RR_internal_loop_data,aes(x=size,y=eng))+geom_point()+
  geom_abline(intercept=138,slope = 11)+xlim(0,30)

RI_RR_ter_AU_penalty = -45
RI_RR_init_eng = -409
RI_RR_bulge_initiation = -200
RI_RR_bulge_extension = -40
RI_RR_internal_loop_asymetry = -60
RI_RR_internal_loop_initiation = -138
RI_RR_internal_loop_symetric_extension = -22
RI_RR_internal_loop_asymetric_extension = -11


################################################
########################## DNA - DNA matrix ###############################
######### WARNING all energies are already multiplied with -100 
dd_commons = NULL
for (d_nt in d_nts){
  dd_commons = rbind(dd_commons, data.frame(x=paste("-",d_nt,sep=""), y=paste(d_nts,"-",sep=""), num=NA))
  dd_commons = rbind(dd_commons, data.frame(x=paste(d_nt,"-",sep=""), y=paste("-",d_nts,sep=""), num=NA))
}
dd_commons$num[dd_commons$x=="--"&dd_commons$y=="--"]=-Inf
dd_commons = rbind(dd_commons, data.frame(x="N-", y="N-", num=-Inf))
dd_commons=unique(dd_commons)

SL_H_DD_stack = NULL
SL_H_DD_stack = rbind(SL_H_DD_stack, data.frame(x="AA",y="TT",num=100))
SL_H_DD_stack = rbind(SL_H_DD_stack, data.frame(x="AC",y="TG",num=144)) #
SL_H_DD_stack = rbind(SL_H_DD_stack, data.frame(x="AG",y="TC",num=128)) #
SL_H_DD_stack = rbind(SL_H_DD_stack, data.frame(x="AT",y="TA",num=88)) 
SL_H_DD_stack = rbind(SL_H_DD_stack, data.frame(x="CA",y="GT",num=145)) 
SL_H_DD_stack = rbind(SL_H_DD_stack, data.frame(x="CC",y="GG",num=184)) #
SL_H_DD_stack = rbind(SL_H_DD_stack, data.frame(x="CG",y="GC",num=217))
SL_H_DD_stack = rbind(SL_H_DD_stack, data.frame(x="CT",y="GA",num=128)) 
SL_H_DD_stack = rbind(SL_H_DD_stack, data.frame(x="GA",y="CT",num=130))
SL_H_DD_stack = rbind(SL_H_DD_stack, data.frame(x="GC",y="CG",num=224))
SL_H_DD_stack = rbind(SL_H_DD_stack, data.frame(x="GG",y="CC",num=184))
SL_H_DD_stack = rbind(SL_H_DD_stack, data.frame(x="GT",y="CA",num=144))
SL_H_DD_stack = rbind(SL_H_DD_stack, data.frame(x="TA",y="AT",num=58))
SL_H_DD_stack = rbind(SL_H_DD_stack, data.frame(x="TC",y="AG",num=130)) # 
SL_H_DD_stack = rbind(SL_H_DD_stack, data.frame(x="TG",y="AC",num=145)) #
SL_H_DD_stack = rbind(SL_H_DD_stack, data.frame(x="TT",y="AA",num=100)) #

SL_H_DD_ter_AT_penalty=-5
SL_H_DD_ter_AT = NULL
SL_H_DD_ter_AT = rbind(SL_H_DD_ter_AT, data.frame(x=paste(c("A","T"),"-",sep=""),y=paste(c("T","A"),"-",sep=""),num=SL_H_DD_ter_AT_penalty))
SL_H_DD_ter_AT = rbind(SL_H_DD_ter_AT, data.frame(x=paste(c("G","C"),"-",sep=""),y=paste(c("C","G"),"-",sep=""),num=0))
SL_H_DD_ter_AT = rbind(SL_H_DD_ter_AT, data.frame(x=paste(c("A","A","A","C","C","C","G","G","G","T","T","T","A","C","G","T","N","N","N","N"),"-",sep=""),
                                                  y=paste(c("A","C","G","A","C","T","A","G","T","C","G","T","N","N","N","N","A","C","G","T"),"-",sep=""),
                                                  num=-Inf))
SL_H_DD_init_eng = -196
SL_H_DD_init = NULL
SL_H_DD_init = rbind(SL_H_DD_init, data.frame(x=paste("-",c("C","G"),sep=""),y=paste("-",c("G","C"),sep=""),num=(ris2_minus+SL_H_DD_init_eng )))
SL_H_DD_init = rbind(SL_H_DD_init, data.frame(x=paste("-",c("A","T"),sep=""),y=paste("-",c("T","A"),sep=""),num=(ris2_minus+SL_H_DD_init_eng+SL_H_DD_ter_AT_penalty)))
SL_H_DD_init = rbind(SL_H_DD_init, data.frame(x=paste("-",c("A","A","A","C","C","C","G","G","G","T","T","T","A","C","G","T","N","N","N","N","N"),sep=""),
                                              y=paste("-",c("A","C","G","A","C","T","A","G","T","C","G","T","N","N","N","N","A","C","G","T","N"),sep=""),num=0))


SL_H_DD_bulge_data = data.frame(size=c(2,3,4,5,6,7,8,9,10,12,14,16,18,20,25,30),eng=c(290,310,320,330,350,370,390,410,430,450,480,500,520,530,560,590))
for (i in 3:20){
  test=lm(SL_H_DD_bulge_data$eng[1:i] ~ SL_H_DD_bulge_data$size[1:i])
  print(ggplot(SL_H_DD_bulge_data,aes(x=size,y=eng))+geom_point()+geom_abline(intercept=round(test$coefficients[1]),slope = round(test$coefficients[2]))+
          ggtitle(paste("Bulge",as.character(i)))+ylim(0,700))
  print(paste(as.character(i),as.character(round(test$coefficients[1])),(as.character(round(test$coefficients[2])))))
}
ggplot(SL_H_DD_bulge_data,aes(x=size,y=eng))+geom_point()+
  geom_abline(intercept=261,slope = 16)+ggtitle(i)+ylim(0,700)+xlim(0,30)
# bulges intercept=261 slope=16 selected 
# therefore -261 will be added to bulge initiations
# HOWEVER A-T closure will also be added on top of it. -261-5=-266 
# bulge extensions will be -16
# Single bulges cannot be treated specially
SL_H_DD_bulge_initiation = -261
SL_H_DD_bulge_extension = -16
SL_H_DD_bulges = NULL
SL_H_DD_bulges = rbind(SL_H_DD_bulges, data.frame(x="A-",y=c(paste("T",d_nts[d_nts!='-'],sep = "")),num=SL_H_DD_bulge_initiation+SL_H_DD_bulge_extension+SL_H_DD_ter_AT_penalty))
SL_H_DD_bulges = rbind(SL_H_DD_bulges, data.frame(x="T-",y=c(paste("A",d_nts[d_nts!='-'],sep = "")),num=SL_H_DD_bulge_initiation+SL_H_DD_bulge_extension+SL_H_DD_ter_AT_penalty))
SL_H_DD_bulges = rbind(SL_H_DD_bulges, data.frame(x="-A",y=c(paste(d_nts[d_nts!='-'],"T",sep = "")),num=SL_H_DD_ter_AT_penalty))
SL_H_DD_bulges = rbind(SL_H_DD_bulges, data.frame(x="-T",y=c(paste(d_nts[d_nts!='-'],"A",sep = "")),num=SL_H_DD_ter_AT_penalty))
SL_H_DD_bulges = rbind(SL_H_DD_bulges, data.frame(x="C-",y=c(paste("G",d_nts[d_nts!='-'],sep = "")),num=SL_H_DD_bulge_initiation+SL_H_DD_bulge_extension))
SL_H_DD_bulges = rbind(SL_H_DD_bulges, data.frame(x="G-",y=c(paste("C",d_nts[d_nts!='-'],sep = "")),num=SL_H_DD_bulge_initiation+SL_H_DD_bulge_extension))
SL_H_DD_bulges = rbind(SL_H_DD_bulges, data.frame(x="-C",y=c(paste(d_nts[d_nts!='-'],"G",sep = "")),num=0))
SL_H_DD_bulges = rbind(SL_H_DD_bulges, data.frame(x="-G",y=c(paste(d_nts[d_nts!='-'],"C",sep = "")),num=0))
SL_H_DD_bulges = rbind(SL_H_DD_bulges, data.frame(y="A-",x=c(paste("T",d_nts[d_nts!='-'],sep = "")),num=SL_H_DD_bulge_initiation+SL_H_DD_bulge_extension+SL_H_DD_ter_AT_penalty))
SL_H_DD_bulges = rbind(SL_H_DD_bulges, data.frame(y="T-",x=c(paste("A",d_nts[d_nts!='-'],sep = "")),num=SL_H_DD_bulge_initiation+SL_H_DD_bulge_extension+SL_H_DD_ter_AT_penalty))
SL_H_DD_bulges = rbind(SL_H_DD_bulges, data.frame(y="-A",x=c(paste(d_nts[d_nts!='-'],"T",sep = "")),num=SL_H_DD_ter_AT_penalty))
SL_H_DD_bulges = rbind(SL_H_DD_bulges, data.frame(y="-T",x=c(paste(d_nts[d_nts!='-'],"A",sep = "")),num=SL_H_DD_ter_AT_penalty))
SL_H_DD_bulges = rbind(SL_H_DD_bulges, data.frame(y="C-",x=c(paste("G",d_nts[d_nts!='-'],sep = "")),num=SL_H_DD_bulge_initiation+SL_H_DD_bulge_extension))
SL_H_DD_bulges = rbind(SL_H_DD_bulges, data.frame(y="G-",x=c(paste("C",d_nts[d_nts!='-'],sep = "")),num=SL_H_DD_bulge_initiation+SL_H_DD_bulge_extension))
SL_H_DD_bulges = rbind(SL_H_DD_bulges, data.frame(y="-C",x=c(paste(d_nts[d_nts!='-'],"G",sep = "")),num=0))
SL_H_DD_bulges = rbind(SL_H_DD_bulges, data.frame(y="-G",x=c(paste(d_nts[d_nts!='-'],"C",sep = "")),num=0))
for (d_nt in d_nts[d_nts!='-']){
  SL_H_DD_bulges = rbind(SL_H_DD_bulges, data.frame(x=paste(d_nts[d_nts!='-'],d_nt,sep = ""),y="--",num=SL_H_DD_bulge_extension))
  SL_H_DD_bulges = rbind(SL_H_DD_bulges, data.frame(x="--",y=paste(d_nts[d_nts!='-'],d_nt,sep = ""),num=SL_H_DD_bulge_extension))
}



SL_H_DD_internal_loop_data = data.frame(size=c(3,4,5,6,7,8,9,10,12,14,16,18,20,25,30),eng=c(320,360,400,440,460,480,490,490,520,540,560,580,590,630,660))
for (i in 3:length(SL_H_DD_internal_loop_data$size)){
  test=lm(SL_H_DD_internal_loop_data$eng[1:i] ~ SL_H_DD_internal_loop_data$size[1:i])
  print(ggplot(SL_H_DD_internal_loop_data,aes(x=size,y=eng))+geom_point()+geom_abline(intercept=round(test$coefficients[1]),slope = round(test$coefficients[2]))+
          ggtitle(paste("Internal Loop",as.character(i)))+ylim(0,700))
  print(paste(as.character(i),as.character(round(test$coefficients[1])),(as.character(round(test$coefficients[2])))))
}
ggplot(SL_H_DD_internal_loop_data,aes(x=size,y=eng))+geom_point()+
  geom_abline(intercept=293,slope = 22)+ylim(0,700)

ggplot(SL_H_DD_internal_loop_data,aes(x=size,y=eng))+geom_point()+
  geom_abline(intercept=231,slope = 33)+ylim(0,700)+xlim(0,32)
# internal loops intercept=293 slope=22 selected 
# therefore -293 will be added to loop initiations 
# symetric loop extensions will be -22*2
# for asymetric loop extensions bulge extension will be applied except the first one
SL_H_DD_internal_loop_asymetry = -30
SL_H_DD_internal_loop_initiation = -231
SL_H_DD_internal_loop_symetric_extension = -66
SL_H_DD_internal_loop_asymetric_extension = -33
SL_H_DD_loops = NULL
for (d_nt in d_nts){
  for (d_nt2 in d_nts){
    mm = paste(d_nt,d_nt2,sep = "")
    if (!(mm%in%c(d_matches,"--"))) {
      for (d_nt3 in d_nts[d_nts!='-']){
        for (d_nt4 in d_nts[d_nts!='-']){
          mm = paste(d_nt3,d_nt4,sep = "")
          if (!(mm%in%c(d_matches,"--"))) {
            SL_H_DD_loops = rbind(SL_H_DD_loops, data.frame(x=paste(d_nt,d_nt3,sep=""),y=paste(d_nt2,d_nt4,sep=""),num=SL_H_DD_internal_loop_symetric_extension))
          }
        }
      }
    }
  }
}

for (d_nt in d_nts[d_nts!='-']){
  for (d_nt2 in d_nts[d_nts!='-']){
    mm = paste(d_nt,d_nt2,sep = "")
    if (!(mm%in%c(d_matches,"--"))) {
      for (d_nt3 in d_nts[d_nts!='-']){
        SL_H_DD_loops = rbind(SL_H_DD_loops, data.frame(x=paste(d_nt,'-',sep=""),y=paste(d_nt2,d_nt3,sep=""),
                                                        num=SL_H_DD_internal_loop_asymetric_extension+SL_H_DD_internal_loop_asymetry))
        SL_H_DD_loops = rbind(SL_H_DD_loops, data.frame(y=paste(d_nt,'-',sep=""),x=paste(d_nt2,d_nt3,sep=""),
                                                        num=SL_H_DD_internal_loop_asymetric_extension+SL_H_DD_internal_loop_asymetry))
      }
    }
  }
}

SL_H_DD_mm = NULL
SL_H_DD_mm = rbind(SL_H_DD_mm, data.frame(x=c("GA","AC"),y=c("CA","AG"),num=c(-17+(SL_H_DD_internal_loop_initiation+SL_H_DD_internal_loop_symetric_extension),-17)))
SL_H_DD_mm = rbind(SL_H_DD_mm, data.frame(x=c("GA","CC"),y=c("CC","AG"),num=c(-81+(SL_H_DD_internal_loop_initiation+SL_H_DD_internal_loop_symetric_extension),-81)))
SL_H_DD_mm = rbind(SL_H_DD_mm, data.frame(x=c("GA","GC"),y=c("CG","AG"),num=c(25+(SL_H_DD_internal_loop_initiation+SL_H_DD_internal_loop_symetric_extension),25)))
SL_H_DD_mm = rbind(SL_H_DD_mm, data.frame(x=c("GC","AC"),y=c("CA","CG"),num=c(-47+(SL_H_DD_internal_loop_initiation+SL_H_DD_internal_loop_symetric_extension),-47)))
SL_H_DD_mm = rbind(SL_H_DD_mm, data.frame(x=c("GC","CC"),y=c("CC","CG"),num=c(-79+(SL_H_DD_internal_loop_initiation+SL_H_DD_internal_loop_symetric_extension),-79)))
SL_H_DD_mm = rbind(SL_H_DD_mm, data.frame(x=c("GC","TC"),y=c("CT","CG"),num=c(-62+(SL_H_DD_internal_loop_initiation+SL_H_DD_internal_loop_symetric_extension),-62)))
SL_H_DD_mm = rbind(SL_H_DD_mm, data.frame(x=c("GG","AC"),y=c("CA","GG"),num=c(52+(SL_H_DD_internal_loop_initiation+SL_H_DD_internal_loop_symetric_extension),52)))
SL_H_DD_mm = rbind(SL_H_DD_mm, data.frame(x=c("GG","GC"),y=c("CG","GG"),num=c(111+(SL_H_DD_internal_loop_initiation+SL_H_DD_internal_loop_symetric_extension),111)))
SL_H_DD_mm = rbind(SL_H_DD_mm, data.frame(x=c("GG","TC"),y=c("CT","GG"),num=c(-8+(SL_H_DD_internal_loop_initiation+SL_H_DD_internal_loop_symetric_extension),-8)))
SL_H_DD_mm = rbind(SL_H_DD_mm, data.frame(x=c("GT","CC"),y=c("CC","TG"),num=c(-98+(SL_H_DD_internal_loop_initiation+SL_H_DD_internal_loop_symetric_extension),-98)))
SL_H_DD_mm = rbind(SL_H_DD_mm, data.frame(x=c("GT","GC"),y=c("CG","TG"),num=c(59+(SL_H_DD_internal_loop_initiation+SL_H_DD_internal_loop_symetric_extension),59)))
SL_H_DD_mm = rbind(SL_H_DD_mm, data.frame(x=c("GT","TC"),y=c("CT","TG"),num=c(-45+(SL_H_DD_internal_loop_initiation+SL_H_DD_internal_loop_symetric_extension),-45)))

SL_H_DD_mm = rbind(SL_H_DD_mm, data.frame(x=c("CA","AG"),y=c("GA","AC"),num=c(-43+(SL_H_DD_internal_loop_initiation+SL_H_DD_internal_loop_symetric_extension),-43)))
SL_H_DD_mm = rbind(SL_H_DD_mm, data.frame(x=c("CA","CG"),y=c("GC","AC"),num=c(-75+(SL_H_DD_internal_loop_initiation+SL_H_DD_internal_loop_symetric_extension),-75)))
SL_H_DD_mm = rbind(SL_H_DD_mm, data.frame(x=c("CA","GG"),y=c("GG","AC"),num=c(-3+(SL_H_DD_internal_loop_initiation+SL_H_DD_internal_loop_symetric_extension),-3)))
SL_H_DD_mm = rbind(SL_H_DD_mm, data.frame(x=c("CC","AG"),y=c("GA","CC"),num=c(-79+(SL_H_DD_internal_loop_initiation+SL_H_DD_internal_loop_symetric_extension),-79)))
SL_H_DD_mm = rbind(SL_H_DD_mm, data.frame(x=c("CC","CG"),y=c("GC","CC"),num=c(-70+(SL_H_DD_internal_loop_initiation+SL_H_DD_internal_loop_symetric_extension),-70)))
SL_H_DD_mm = rbind(SL_H_DD_mm, data.frame(x=c("CC","TG"),y=c("GT","CC"),num=c(-62+(SL_H_DD_internal_loop_initiation+SL_H_DD_internal_loop_symetric_extension),-62)))
SL_H_DD_mm = rbind(SL_H_DD_mm, data.frame(x=c("CG","AG"),y=c("GA","GC"),num=c(-11+(SL_H_DD_internal_loop_initiation+SL_H_DD_internal_loop_symetric_extension),-11)))
SL_H_DD_mm = rbind(SL_H_DD_mm, data.frame(x=c("CG","GG"),y=c("GG","GC"),num=c(11+(SL_H_DD_internal_loop_initiation+SL_H_DD_internal_loop_symetric_extension),11)))
SL_H_DD_mm = rbind(SL_H_DD_mm, data.frame(x=c("CG","TG"),y=c("GT","GC"),num=c(47+(SL_H_DD_internal_loop_initiation+SL_H_DD_internal_loop_symetric_extension),47)))
SL_H_DD_mm = rbind(SL_H_DD_mm, data.frame(x=c("CT","CG"),y=c("GC","TC"),num=c(-40+(SL_H_DD_internal_loop_initiation+SL_H_DD_internal_loop_symetric_extension),-40)))
SL_H_DD_mm = rbind(SL_H_DD_mm, data.frame(x=c("CT","GG"),y=c("GG","TC"),num=c(32+(SL_H_DD_internal_loop_initiation+SL_H_DD_internal_loop_symetric_extension),32)))
SL_H_DD_mm = rbind(SL_H_DD_mm, data.frame(x=c("CT","TG"),y=c("GT","TC"),num=c(12+(SL_H_DD_internal_loop_initiation+SL_H_DD_internal_loop_symetric_extension),12)))

SL_H_DD_mm = rbind(SL_H_DD_mm, data.frame(x=c("AA","AT"),y=c("TA","AA"),num=c(-61+(SL_H_DD_internal_loop_initiation+SL_H_DD_internal_loop_symetric_extension),-61)+SL_H_DD_ter_AT_penalty))
SL_H_DD_mm = rbind(SL_H_DD_mm, data.frame(x=c("AA","CT"),y=c("TC","AA"),num=c(-88+(SL_H_DD_internal_loop_initiation+SL_H_DD_internal_loop_symetric_extension),-88)+SL_H_DD_ter_AT_penalty))
SL_H_DD_mm = rbind(SL_H_DD_mm, data.frame(x=c("AA","GT"),y=c("TG","AA"),num=c(-14+(SL_H_DD_internal_loop_initiation+SL_H_DD_internal_loop_symetric_extension),-14)+SL_H_DD_ter_AT_penalty))
SL_H_DD_mm = rbind(SL_H_DD_mm, data.frame(x=c("AC","AT"),y=c("TA","CA"),num=c(-77+(SL_H_DD_internal_loop_initiation+SL_H_DD_internal_loop_symetric_extension),-77)+SL_H_DD_ter_AT_penalty))
SL_H_DD_mm = rbind(SL_H_DD_mm, data.frame(x=c("AC","CT"),y=c("TC","CA"),num=c(-133+(SL_H_DD_internal_loop_initiation+SL_H_DD_internal_loop_symetric_extension),-133)+SL_H_DD_ter_AT_penalty))
SL_H_DD_mm = rbind(SL_H_DD_mm, data.frame(x=c("AC","TT"),y=c("TT","CA"),num=c(-64+(SL_H_DD_internal_loop_initiation+SL_H_DD_internal_loop_symetric_extension),-64)+SL_H_DD_ter_AT_penalty))
SL_H_DD_mm = rbind(SL_H_DD_mm, data.frame(x=c("AG","AT"),y=c("TA","GA"),num=c(-2+(SL_H_DD_internal_loop_initiation+SL_H_DD_internal_loop_symetric_extension),-2)+SL_H_DD_ter_AT_penalty))
SL_H_DD_mm = rbind(SL_H_DD_mm, data.frame(x=c("AG","GT"),y=c("TG","GA"),num=c(13+(SL_H_DD_internal_loop_initiation+SL_H_DD_internal_loop_symetric_extension),13)+SL_H_DD_ter_AT_penalty))
SL_H_DD_mm = rbind(SL_H_DD_mm, data.frame(x=c("AG","TT"),y=c("TT","GA"),num=c(-71+(SL_H_DD_internal_loop_initiation+SL_H_DD_internal_loop_symetric_extension),-71)+SL_H_DD_ter_AT_penalty))
SL_H_DD_mm = rbind(SL_H_DD_mm, data.frame(x=c("AT","CT"),y=c("TC","TA"),num=c(-73+(SL_H_DD_internal_loop_initiation+SL_H_DD_internal_loop_symetric_extension),-73)+SL_H_DD_ter_AT_penalty))
SL_H_DD_mm = rbind(SL_H_DD_mm, data.frame(x=c("AT","GT"),y=c("TG","TA"),num=c(-7+(SL_H_DD_internal_loop_initiation+SL_H_DD_internal_loop_symetric_extension),-7)+SL_H_DD_ter_AT_penalty))
SL_H_DD_mm = rbind(SL_H_DD_mm, data.frame(x=c("AT","TT"),y=c("TT","TA"),num=c(-71+(SL_H_DD_internal_loop_initiation+SL_H_DD_internal_loop_symetric_extension),-71)+SL_H_DD_ter_AT_penalty))

SL_H_DD_mm = rbind(SL_H_DD_mm, data.frame(x=c("TA","AA"),y=c("AA","AT"),num=c(-71+(SL_H_DD_internal_loop_initiation+SL_H_DD_internal_loop_symetric_extension),-71)+SL_H_DD_ter_AT_penalty))
SL_H_DD_mm = rbind(SL_H_DD_mm, data.frame(x=c("TA","CA"),y=c("AC","AT"),num=c(-92+(SL_H_DD_internal_loop_initiation+SL_H_DD_internal_loop_symetric_extension),-92)+SL_H_DD_ter_AT_penalty))
SL_H_DD_mm = rbind(SL_H_DD_mm, data.frame(x=c("TA","GA"),y=c("AG","AT"),num=c(-42+(SL_H_DD_internal_loop_initiation+SL_H_DD_internal_loop_symetric_extension),-42)+SL_H_DD_ter_AT_penalty))
SL_H_DD_mm = rbind(SL_H_DD_mm, data.frame(x=c("TC","AA"),y=c("AA","CT"),num=c(-133+(SL_H_DD_internal_loop_initiation+SL_H_DD_internal_loop_symetric_extension),-133)+SL_H_DD_ter_AT_penalty))
SL_H_DD_mm = rbind(SL_H_DD_mm, data.frame(x=c("TC","CA"),y=c("AC","CT"),num=c(-105+(SL_H_DD_internal_loop_initiation+SL_H_DD_internal_loop_symetric_extension),-105)+SL_H_DD_ter_AT_penalty))
SL_H_DD_mm = rbind(SL_H_DD_mm, data.frame(x=c("TC","TA"),y=c("AT","CT"),num=c(-97+(SL_H_DD_internal_loop_initiation+SL_H_DD_internal_loop_symetric_extension),-97)+SL_H_DD_ter_AT_penalty))
SL_H_DD_mm = rbind(SL_H_DD_mm, data.frame(x=c("TG","AA"),y=c("AA","GT"),num=c(-74+(SL_H_DD_internal_loop_initiation+SL_H_DD_internal_loop_symetric_extension),-74)+SL_H_DD_ter_AT_penalty))
SL_H_DD_mm = rbind(SL_H_DD_mm, data.frame(x=c("TG","GA"),y=c("AG","GT"),num=c(-44+(SL_H_DD_internal_loop_initiation+SL_H_DD_internal_loop_symetric_extension),-44)+SL_H_DD_ter_AT_penalty))
SL_H_DD_mm = rbind(SL_H_DD_mm, data.frame(x=c("TG","TA"),y=c("AT","GT"),num=c(-43+(SL_H_DD_internal_loop_initiation+SL_H_DD_internal_loop_symetric_extension),-43)+SL_H_DD_ter_AT_penalty))
SL_H_DD_mm = rbind(SL_H_DD_mm, data.frame(x=c("TT","CA"),y=c("AC","TT"),num=c(-75+(SL_H_DD_internal_loop_initiation+SL_H_DD_internal_loop_symetric_extension),-75)+SL_H_DD_ter_AT_penalty))
SL_H_DD_mm = rbind(SL_H_DD_mm, data.frame(x=c("TT","GA"),y=c("AG","TT"),num=c(-34+(SL_H_DD_internal_loop_initiation+SL_H_DD_internal_loop_symetric_extension),-34)+SL_H_DD_ter_AT_penalty))
SL_H_DD_mm = rbind(SL_H_DD_mm, data.frame(x=c("TT","TA"),y=c("AT","TT"),num=c(-68+(SL_H_DD_internal_loop_initiation+SL_H_DD_internal_loop_symetric_extension),-68)+SL_H_DD_ter_AT_penalty))

SL_H_DD_mm = rbind(SL_H_DD_mm, data.frame(x="AN",y=paste("T",d_nts[d_nts!='-'],sep=""),num=(-133)+SL_H_DD_internal_loop_initiation+SL_H_DD_internal_loop_symetric_extension+SL_H_DD_ter_AT_penalty))
SL_H_DD_mm = rbind(SL_H_DD_mm, data.frame(x="CN",y=paste("G",d_nts[d_nts!='-'],sep=""),num=(-98)+SL_H_DD_internal_loop_initiation+SL_H_DD_internal_loop_symetric_extension))
SL_H_DD_mm = rbind(SL_H_DD_mm, data.frame(x="GN",y=paste("C",d_nts[d_nts!='-'],sep=""),num=(-98)+SL_H_DD_internal_loop_initiation+SL_H_DD_internal_loop_symetric_extension))
SL_H_DD_mm = rbind(SL_H_DD_mm, data.frame(x="TN",y=paste("A",d_nts[d_nts!='-'],sep=""),num=(-133)+SL_H_DD_internal_loop_initiation+SL_H_DD_internal_loop_symetric_extension+SL_H_DD_ter_AT_penalty))
SL_H_DD_mm = rbind(SL_H_DD_mm, data.frame(y="AN",x=paste("T",d_nts[!d_nts%in%c('-','N')],sep=""),num=(-133)+SL_H_DD_internal_loop_initiation+SL_H_DD_internal_loop_symetric_extension+SL_H_DD_ter_AT_penalty))
SL_H_DD_mm = rbind(SL_H_DD_mm, data.frame(y="CN",x=paste("G",d_nts[!d_nts%in%c('-','N')],sep=""),num=(-98)+SL_H_DD_internal_loop_initiation+SL_H_DD_internal_loop_symetric_extension))
SL_H_DD_mm = rbind(SL_H_DD_mm, data.frame(y="GN",x=paste("C",d_nts[!d_nts%in%c('-','N')],sep=""),num=(-98)+SL_H_DD_internal_loop_initiation+SL_H_DD_internal_loop_symetric_extension))
SL_H_DD_mm = rbind(SL_H_DD_mm, data.frame(y="TN",x=paste("A",d_nts[!d_nts%in%c('-','N')],sep=""),num=(-133)+SL_H_DD_internal_loop_initiation+SL_H_DD_internal_loop_symetric_extension+SL_H_DD_ter_AT_penalty))

SL_H_DD_mm = rbind(SL_H_DD_mm, data.frame(x="NA",y=paste(d_nts[d_nts!='-'],"T",sep=""),num=SL_H_DD_ter_AT_penalty))
SL_H_DD_mm = rbind(SL_H_DD_mm, data.frame(x="NC",y=paste(d_nts[d_nts!='-'],"G",sep=""),num=0))
SL_H_DD_mm = rbind(SL_H_DD_mm, data.frame(x="NG",y=paste(d_nts[d_nts!='-'],"C",sep=""),num=0))
SL_H_DD_mm = rbind(SL_H_DD_mm, data.frame(x="NT",y=paste(d_nts[d_nts!='-'],"A",sep=""),num=SL_H_DD_ter_AT_penalty))
SL_H_DD_mm = rbind(SL_H_DD_mm, data.frame(y="NA",x=paste(d_nts[!d_nts%in%c('-','N')],"T",sep=""),num=SL_H_DD_ter_AT_penalty))
SL_H_DD_mm = rbind(SL_H_DD_mm, data.frame(y="NC",x=paste(d_nts[!d_nts%in%c('-','N')],"G",sep=""),num=0))
SL_H_DD_mm = rbind(SL_H_DD_mm, data.frame(y="NG",x=paste(d_nts[!d_nts%in%c('-','N')],"C",sep=""),num=0))
SL_H_DD_mm = rbind(SL_H_DD_mm, data.frame(y="NT",x=paste(d_nts[!d_nts%in%c('-','N')],"A",sep=""),num=SL_H_DD_ter_AT_penalty))


SL_H_DD_all = rbind(SL_H_DD_stack, SL_H_DD_init, SL_H_DD_ter_AT, SL_H_DD_mm, dd_commons, SL_H_DD_bulges, SL_H_DD_loops)
SL_H_DD_all$x=as.character(SL_H_DD_all$x)
SL_H_DD_all$y=as.character(SL_H_DD_all$y)

for (i in 0:length(d_nns)){
  if (i==0) {cat("    ",paste(d_nns,collapse = "   "))}
  for (j in 0:length(d_nns)){
    if (j==0) {cat(d_nns[i])}
    else{
      val = as.character(SL_H_DD_all$num[SL_H_DD_all$x==d_nns[i] & SL_H_DD_all$y==d_nns[j]])
      if (length(val)==0) {val=""}
      else if (is.na(val) | val=="NA") {val="*"}
      else if (val==-Inf) {val="-inf"}
      valc=nchar(val)
      cat(paste(rep(" ",(4-valc)), collapse=""),val)
    }
  }
  cat("\n")
}



################################################
########################## DNA - RNA matrix ######################
######### WARNING all energies are already multiplied with -100
# mean_RD_all=NULL
# for (i in 1:36){
#   for (j in 1:36){
#     rr_num = RI_RR_all$num[RI_RR_all$x==r_nns[i] & RI_RR_all$y==r_nns[j]]
#     dd_num = SL_H_DD_all$num[SL_H_DD_all$x==d_nns[i] & SL_H_DD_all$y==d_nns[j]]
#     if (is.infinite(rr_num))
#       mean_RD_all=rbind(mean_RD_all,data.frame(x=r_nns[i],y=d_nns[j],num=-Inf))
#     else if (is.na(rr_num))
#       mean_RD_all=rbind(mean_RD_all,data.frame(x=r_nns[i],y=d_nns[j],num=NA))
#     else
#       mean_RD_all=rbind(mean_RD_all,data.frame(x=r_nns[i],y=d_nns[j],num=round((rr_num+dd_num)/2)))
#   }
# }
# 
# for (i in 0:length(r_nns)){
#   if (i==0) {cat("    ",paste(d_nns,collapse = "   "))}
#   for (j in 0:length(d_nns)){
#     if (j==0) {cat(r_nns[i])}
#     else{
#       val = as.character(mean_RD_all$num[mean_RD_all$x==r_nns[i] & mean_RD_all$y==d_nns[j]])
#       if (length(val)==0) {val=""}
#       else if (is.na(val) | val=="NA") {val="*"}
#       else if (val==-Inf) {val="-inf"}
#       valc=nchar(val)
#       cat(paste(rep(" ",(4-valc)), collapse=""),val)
#     }
#   }
#   cat("\n")
# }
###########################################################################################
rd_commons = NULL
for (i in 1:6){
  rd_commons = rbind(rd_commons, data.frame(x=paste("-",r_nts[i],sep=""), y=paste(d_nts,"-",sep=""), num=NA))
  rd_commons = rbind(rd_commons, data.frame(x=paste(r_nts[i],"-",sep=""), y=paste("-",d_nts,sep=""), num=NA))
}
rd_commons$num[rd_commons$x=="--"&rd_commons$y=="--"]=-Inf
rd_commons = rbind(rd_commons, data.frame(x="N-", y="N-", num=-Inf))
rd_commons=unique(rd_commons)

############################################################################################
SU_RD_stack = NULL
SU_RD_stack = rbind(SU_RD_stack, data.frame(x="AA",y="TT",num=100))
SU_RD_stack = rbind(SU_RD_stack, data.frame(x="AC",y="TG",num=210)) #
SU_RD_stack = rbind(SU_RD_stack, data.frame(x="AG",y="TC",num=180)) #
SU_RD_stack = rbind(SU_RD_stack, data.frame(x="AU",y="TA",num=90)) 
SU_RD_stack = rbind(SU_RD_stack, data.frame(x="CA",y="GT",num=90)) 
SU_RD_stack = rbind(SU_RD_stack, data.frame(x="CC",y="GG",num=210)) #
SU_RD_stack = rbind(SU_RD_stack, data.frame(x="CG",y="GC",num=170))
SU_RD_stack = rbind(SU_RD_stack, data.frame(x="CU",y="GA",num=90)) 
SU_RD_stack = rbind(SU_RD_stack, data.frame(x="GA",y="CT",num=130))
SU_RD_stack = rbind(SU_RD_stack, data.frame(x="GC",y="CG",num=270))
SU_RD_stack = rbind(SU_RD_stack, data.frame(x="GG",y="CC",num=290))
SU_RD_stack = rbind(SU_RD_stack, data.frame(x="GU",y="CA",num=110))
SU_RD_stack = rbind(SU_RD_stack, data.frame(x="UA",y="AT",num=60))
SU_RD_stack = rbind(SU_RD_stack, data.frame(x="UC",y="AG",num=150)) # 
SU_RD_stack = rbind(SU_RD_stack, data.frame(x="UG",y="AC",num=160)) #
SU_RD_stack = rbind(SU_RD_stack, data.frame(x="UU",y="AA",num=20)) #

SU_RD_ter_AT_penalty = round((SL_H_DD_ter_AT_penalty+RI_RR_ter_AU_penalty)/2)
SU_RD_ter_AU_penalty = round((SL_H_DD_ter_AT_penalty+RI_RR_ter_AU_penalty)/2)

SU_RD_ter_AT = NULL
SU_RD_ter_AT = rbind(SU_RD_ter_AT, data.frame(x="A-",y="T-",num=SU_RD_ter_AT_penalty))
SU_RD_ter_AT = rbind(SU_RD_ter_AT, data.frame(x="U-",y="A-",num=SU_RD_ter_AU_penalty))
#SU_RD_ter_AT = rbind(SU_RD_ter_AT, data.frame(x="U-",y="G-",num=SU_RD_ter_GU_penalty))
SU_RD_ter_AT = rbind(SU_RD_ter_AT, data.frame(x=paste(c("G","C"),"-",sep=""),y=paste(c("C","G"),"-",sep=""),num=0))
SU_RD_ter_AT = rbind(SU_RD_ter_AT, data.frame(x=paste(c("A","A","A","C","C","C","G","G","G","U","U","U","A","C","G","U","N","N","N","N"),"-",sep=""),
                                              y=paste(c("A","C","G","A","C","T","A","G","T","C","G","T","N","N","N","N","A","C","G","T"),"-",sep=""),
                                              num=-Inf))

SU_RD_init_eng = -310
SU_RD_init = NULL
SU_RD_init = rbind(SU_RD_init, data.frame(x=paste("-",c("C","G"),sep=""),y=paste("-",c("G","C"),sep=""),num=(ris2_minus+SU_RD_init_eng )))
SU_RD_init = rbind(SU_RD_init, data.frame(x="-A",y="-T",num=(ris2_minus+SU_RD_init_eng+SU_RD_ter_AT_penalty)))
SU_RD_init = rbind(SU_RD_init, data.frame(x="-U",y="-A",num=(ris2_minus+SU_RD_init_eng+SU_RD_ter_AU_penalty)))
#SU_RD_init = rbind(SU_RD_init, data.frame(x="-U",y="-G",num=(ris2_minus+SU_RD_init_eng+SU_RD_ter_GU_penalty)))
SU_RD_init = rbind(SU_RD_init, data.frame(x=paste("-",c("A","A","A","C","C","C","G","G","G","U","U","U","A","C","G","U","N","N","N","N","N"),sep=""),
                                          y=paste("-",c("A","C","G","A","C","T","A","G","T","C","G","T","N","N","N","N","A","C","G","T","N"),sep=""),num=0))

SU_RD_bulge_initiation = round((RI_RR_bulge_initiation+SL_H_DD_bulge_initiation)/2)
SU_RD_bulge_extension = round((RI_RR_bulge_extension+SL_H_DD_bulge_extension)/2)
SU_RD_bulges = NULL
SU_RD_bulges = rbind(SU_RD_bulges, data.frame(x="A-",y=c(paste("T",d_nts[d_nts!='-'],sep = "")),num=SU_RD_bulge_initiation+SU_RD_bulge_extension+SU_RD_ter_AT_penalty))
SU_RD_bulges = rbind(SU_RD_bulges, data.frame(x="U-",y=c(paste("A",d_nts[d_nts!='-'],sep = "")),num=SU_RD_bulge_initiation+SU_RD_bulge_extension+SU_RD_ter_AT_penalty))
SU_RD_bulges = rbind(SU_RD_bulges, data.frame(x="-A",y=c(paste(d_nts[d_nts!='-'],"T",sep = "")),num=SU_RD_ter_AT_penalty))
SU_RD_bulges = rbind(SU_RD_bulges, data.frame(x="-U",y=c(paste(d_nts[d_nts!='-'],"A",sep = "")),num=SU_RD_ter_AT_penalty))
SU_RD_bulges = rbind(SU_RD_bulges, data.frame(x="C-",y=c(paste("G",d_nts[d_nts!='-'],sep = "")),num=SU_RD_bulge_initiation+SU_RD_bulge_extension))
SU_RD_bulges = rbind(SU_RD_bulges, data.frame(x="G-",y=c(paste("C",d_nts[d_nts!='-'],sep = "")),num=SU_RD_bulge_initiation+SU_RD_bulge_extension))
SU_RD_bulges = rbind(SU_RD_bulges, data.frame(x="-C",y=c(paste(d_nts[d_nts!='-'],"G",sep = "")),num=0))
SU_RD_bulges = rbind(SU_RD_bulges, data.frame(x="-G",y=c(paste(d_nts[d_nts!='-'],"C",sep = "")),num=0))
SU_RD_bulges = rbind(SU_RD_bulges, data.frame(y="A-",x=c(paste("U",r_nts[r_nts!='-'],sep = "")),num=SU_RD_bulge_initiation+SU_RD_bulge_extension+SU_RD_ter_AT_penalty))
SU_RD_bulges = rbind(SU_RD_bulges, data.frame(y="T-",x=c(paste("A",r_nts[r_nts!='-'],sep = "")),num=SU_RD_bulge_initiation+SU_RD_bulge_extension+SU_RD_ter_AT_penalty))
SU_RD_bulges = rbind(SU_RD_bulges, data.frame(y="-A",x=c(paste(r_nts[r_nts!='-'],"U",sep = "")),num=SU_RD_ter_AT_penalty))
SU_RD_bulges = rbind(SU_RD_bulges, data.frame(y="-T",x=c(paste(r_nts[r_nts!='-'],"A",sep = "")),num=SU_RD_ter_AT_penalty))
SU_RD_bulges = rbind(SU_RD_bulges, data.frame(y="C-",x=c(paste("G",r_nts[r_nts!='-'],sep = "")),num=SU_RD_bulge_initiation+SU_RD_bulge_extension))
SU_RD_bulges = rbind(SU_RD_bulges, data.frame(y="G-",x=c(paste("C",r_nts[r_nts!='-'],sep = "")),num=SU_RD_bulge_initiation+SU_RD_bulge_extension))
SU_RD_bulges = rbind(SU_RD_bulges, data.frame(y="-C",x=c(paste(r_nts[r_nts!='-'],"G",sep = "")),num=0))
SU_RD_bulges = rbind(SU_RD_bulges, data.frame(y="-G",x=c(paste(r_nts[r_nts!='-'],"C",sep = "")),num=0))
for (i in 1:5){
  SU_RD_bulges = rbind(SU_RD_bulges, data.frame(x=paste(r_nts[r_nts!='-'],r_nts[i],sep = ""),y="--",num=SU_RD_bulge_extension))
  SU_RD_bulges = rbind(SU_RD_bulges, data.frame(x="--",y=paste(d_nts[d_nts!='-'],d_nts[i],sep = ""),num=SU_RD_bulge_extension))
}

SU_RD_internal_loop_asymetry = round((RI_RR_internal_loop_asymetry+SL_H_DD_internal_loop_asymetry)/2)
SU_RD_internal_loop_initiation = round((RI_RR_internal_loop_initiation+SL_H_DD_internal_loop_initiation)/2)
SU_RD_internal_loop_symetric_extension = round((RI_RR_internal_loop_symetric_extension+SL_H_DD_internal_loop_symetric_extension)/2)
SU_RD_internal_loop_asymetric_extension = round((RI_RR_internal_loop_asymetric_extension+SL_H_DD_internal_loop_asymetric_extension)/2)
SU_RD_loops = NULL
for (i in 1:6){
  for (j in 1:6){
    mm = paste(d_nts[i],d_nts[j],sep = "")
    if (!(mm%in%c(d_matches,"--"))) {
      for (k in 1:5){
        for (l in 1:5){
          mm = paste(d_nts[k],d_nts[l],sep = "")
          if (!(mm%in%c(d_matches,"--"))) {
            SU_RD_loops = rbind(SU_RD_loops, data.frame(x=paste(r_nts[i],r_nts[k],sep=""),y=paste(d_nts[j],d_nts[l],sep=""),num=SU_RD_internal_loop_symetric_extension))
          }
        }
      }
    }
  }
}

for (i in 1:5){
  for (j in 1:5){
    mm = paste(d_nts[i],d_nts[j],sep = "")
    if (!(mm%in%c(d_matches,"--"))) {
      for (k in 1:5){
        SU_RD_loops = rbind(SU_RD_loops, data.frame(x=paste(r_nts[i],'-',sep=""),y=paste(d_nts[j],d_nts[k],sep=""),
                                                    num=SU_RD_internal_loop_asymetric_extension+SU_RD_internal_loop_asymetry))
        SU_RD_loops = rbind(SU_RD_loops, data.frame(y=paste(d_nts[i],'-',sep=""),x=paste(r_nts[j],r_nts[k],sep=""),
                                                    num=SU_RD_internal_loop_asymetric_extension+SU_RD_internal_loop_asymetry))
      }
    }
  }
}


##########################################################################
## RIS2 RNA RNA mismatches unknown, only AG GA GG UU => 80 100 120 70 rest accepted as 0, and GU UG cases from *Turner (8 different energies)
half <- function(a,b){return(round((a+b)/2))}
SU_RD_mm = NULL
SU_RD_mm = rbind(SU_RD_mm, data.frame(x=c("GA","AC"),y=c("CA","AG"),num=c(half(0,-17)+(SU_RD_internal_loop_initiation+SU_RD_internal_loop_symetric_extension),half(0,-17))))
SU_RD_mm = rbind(SU_RD_mm, data.frame(x=c("CA","AG"),y=c("GA","AC"),num=c(half(0,-43)+(SU_RD_internal_loop_initiation+SU_RD_internal_loop_symetric_extension),half(0,-43))))
SU_RD_mm = rbind(SU_RD_mm, data.frame(x=c("AA","AU"),y=c("TA","AA"),num=c(half(0,-61)+(SU_RD_internal_loop_initiation+SU_RD_internal_loop_symetric_extension),half(0,-61))+SU_RD_ter_AU_penalty))
SU_RD_mm = rbind(SU_RD_mm, data.frame(x=c("UA","AA"),y=c("AA","AT"),num=c(half(0,-71)+(SU_RD_internal_loop_initiation+SU_RD_internal_loop_symetric_extension),half(0,-71))+SU_RD_ter_AU_penalty))

SU_RD_mm = rbind(SU_RD_mm, data.frame(x=c("GA","CC"),y=c("CC","AG"),num=c(half(0,-81)+(SU_RD_internal_loop_initiation+SU_RD_internal_loop_symetric_extension),half(0,-81))))
SU_RD_mm = rbind(SU_RD_mm, data.frame(x=c("CA","CG"),y=c("GC","AC"),num=c(half(0,-75)+(SU_RD_internal_loop_initiation+SU_RD_internal_loop_symetric_extension),half(0,-75))))
SU_RD_mm = rbind(SU_RD_mm, data.frame(x=c("AA","CU"),y=c("TC","AA"),num=c(half(0,-88)+(SU_RD_internal_loop_initiation+SU_RD_internal_loop_symetric_extension),half(0,-88))+SU_RD_ter_AU_penalty))
SU_RD_mm = rbind(SU_RD_mm, data.frame(x=c("UA","CA"),y=c("AC","AT"),num=c(half(0,-92)+(SU_RD_internal_loop_initiation+SU_RD_internal_loop_symetric_extension),half(0,-92))+SU_RD_ter_AU_penalty))

#AG terminal mismatch specialty from R-R
SU_RD_mm = rbind(SU_RD_mm, data.frame(x=c("GA","GC"),y=c("CG","AG"),num=c(half(80, 25)+(SU_RD_internal_loop_initiation+SU_RD_internal_loop_symetric_extension),half(80, 25))))
SU_RD_mm = rbind(SU_RD_mm, data.frame(x=c("CA","GG"),y=c("GG","AC"),num=c(half(80, -3)+(SU_RD_internal_loop_initiation+SU_RD_internal_loop_symetric_extension),half(80, -3))))
SU_RD_mm = rbind(SU_RD_mm, data.frame(x=c("AA","GU"),y=c("TG","AA"),num=c(half(80,-14)+(SU_RD_internal_loop_initiation+SU_RD_internal_loop_symetric_extension),half(80,-14))+SU_RD_ter_AU_penalty))
SU_RD_mm = rbind(SU_RD_mm, data.frame(x=c("UA","GA"),y=c("AG","AT"),num=c(half(80,-42)+(SU_RD_internal_loop_initiation+SU_RD_internal_loop_symetric_extension),half(80,-42))+SU_RD_ter_AU_penalty))

SU_RD_mm = rbind(SU_RD_mm, data.frame(x=c("GC","AC"),y=c("CA","CG"),num=c(half(0,-47)+(SU_RD_internal_loop_initiation+SU_RD_internal_loop_symetric_extension),half(0,-47))))
SU_RD_mm = rbind(SU_RD_mm, data.frame(x=c("CC","AG"),y=c("GA","CC"),num=c(half(0,-79)+(SU_RD_internal_loop_initiation+SU_RD_internal_loop_symetric_extension),half(0,-79))))
SU_RD_mm = rbind(SU_RD_mm, data.frame(x=c("AC","AU"),y=c("TA","CA"),num=c(half(0,-77)+(SU_RD_internal_loop_initiation+SU_RD_internal_loop_symetric_extension),half(0,-77))+SU_RD_ter_AU_penalty))
SU_RD_mm = rbind(SU_RD_mm, data.frame(x=c("UC","AA"),y=c("AA","CT"),num=c(half(0,-133)+(SU_RD_internal_loop_initiation+SU_RD_internal_loop_symetric_extension),half(0,-133))+SU_RD_ter_AU_penalty))

SU_RD_mm = rbind(SU_RD_mm, data.frame(x=c("GC","CC"),y=c("CC","CG"),num=c(half(0, -79)+(SU_RD_internal_loop_initiation+SU_RD_internal_loop_symetric_extension),half(0, -79))))
SU_RD_mm = rbind(SU_RD_mm, data.frame(x=c("CC","CG"),y=c("GC","CC"),num=c(half(0, -70)+(SU_RD_internal_loop_initiation+SU_RD_internal_loop_symetric_extension),half(0, -70))))
SU_RD_mm = rbind(SU_RD_mm, data.frame(x=c("AC","CU"),y=c("TC","CA"),num=c(half(0,-133)+(SU_RD_internal_loop_initiation+SU_RD_internal_loop_symetric_extension),half(0,-133))+SU_RD_ter_AU_penalty))
SU_RD_mm = rbind(SU_RD_mm, data.frame(x=c("UC","CA"),y=c("AC","CT"),num=c(half(0,-105)+(SU_RD_internal_loop_initiation+SU_RD_internal_loop_symetric_extension),half(0,-105))+SU_RD_ter_AU_penalty))

SU_RD_mm = rbind(SU_RD_mm, data.frame(x=c("GC","UC"),y=c("CT","CG"),num=c(half(0,-62)+(SU_RD_internal_loop_initiation+SU_RD_internal_loop_symetric_extension),half(0,-62))))
SU_RD_mm = rbind(SU_RD_mm, data.frame(x=c("CC","UG"),y=c("GT","CC"),num=c(half(0,-62)+(SU_RD_internal_loop_initiation+SU_RD_internal_loop_symetric_extension),half(0,-62))))
SU_RD_mm = rbind(SU_RD_mm, data.frame(x=c("AC","UU"),y=c("TT","CA"),num=c(half(0,-64)+(SU_RD_internal_loop_initiation+SU_RD_internal_loop_symetric_extension),half(0,-64))+SU_RD_ter_AU_penalty))
SU_RD_mm = rbind(SU_RD_mm, data.frame(x=c("UC","UA"),y=c("AT","CT"),num=c(half(0,-97)+(SU_RD_internal_loop_initiation+SU_RD_internal_loop_symetric_extension),half(0,-97))+SU_RD_ter_AU_penalty))

#GA terminal mismatch specialty from R-R
SU_RD_mm = rbind(SU_RD_mm, data.frame(x=c("GG","AC"),y=c("CA","GG"),num=c(half(100, 52)+(SU_RD_internal_loop_initiation+SU_RD_internal_loop_symetric_extension),half(100,52)))) 
SU_RD_mm = rbind(SU_RD_mm, data.frame(x=c("CG","AG"),y=c("GA","GC"),num=c(half(100,-11)+(SU_RD_internal_loop_initiation+SU_RD_internal_loop_symetric_extension),half(100,-11))))
SU_RD_mm = rbind(SU_RD_mm, data.frame(x=c("AG","AU"),y=c("TA","GA"),num=c(half(100, -2)+(SU_RD_internal_loop_initiation+SU_RD_internal_loop_symetric_extension),half(100,-2))+SU_RD_ter_AU_penalty))
SU_RD_mm = rbind(SU_RD_mm, data.frame(x=c("UG","AA"),y=c("AA","GT"),num=c(half(100,-74)+(SU_RD_internal_loop_initiation+SU_RD_internal_loop_symetric_extension),half(100,-74))+SU_RD_ter_AU_penalty))

#GG terminal mismatch specialty from R-R
SU_RD_mm = rbind(SU_RD_mm, data.frame(x=c("GG","GC"),y=c("CG","GG"),num=c(half(120,111)+(SU_RD_internal_loop_initiation+SU_RD_internal_loop_symetric_extension),half(120,111)))) 
SU_RD_mm = rbind(SU_RD_mm, data.frame(x=c("CG","GG"),y=c("GG","GC"),num=c(half(120,11)+(SU_RD_internal_loop_initiation+SU_RD_internal_loop_symetric_extension),half(120,11))))
SU_RD_mm = rbind(SU_RD_mm, data.frame(x=c("AG","GU"),y=c("TG","GA"),num=c(half(120,13)+(SU_RD_internal_loop_initiation+SU_RD_internal_loop_symetric_extension),half(120,13))+SU_RD_ter_AU_penalty))
SU_RD_mm = rbind(SU_RD_mm, data.frame(x=c("UG","GA"),y=c("AG","GT"),num=c(half(120,-44)+(SU_RD_internal_loop_initiation+SU_RD_internal_loop_symetric_extension),half(120,-44))+SU_RD_ter_AU_penalty))

SU_RD_mm = rbind(SU_RD_mm, data.frame(x=c("GU","CC"),y=c("CC","TG"),num=c(half(0,-98)+(SU_RD_internal_loop_initiation+SU_RD_internal_loop_symetric_extension),half(0,-98))))
SU_RD_mm = rbind(SU_RD_mm, data.frame(x=c("CU","CG"),y=c("GC","TC"),num=c(half(0,-40)+(SU_RD_internal_loop_initiation+SU_RD_internal_loop_symetric_extension),half(0,-40))))
SU_RD_mm = rbind(SU_RD_mm, data.frame(x=c("AU","CU"),y=c("TC","TA"),num=c(half(0,-73)+(SU_RD_internal_loop_initiation+SU_RD_internal_loop_symetric_extension),half(0,-73))+SU_RD_ter_AU_penalty))
SU_RD_mm = rbind(SU_RD_mm, data.frame(x=c("UU","CA"),y=c("AC","TT"),num=c(half(0,-75)+(SU_RD_internal_loop_initiation+SU_RD_internal_loop_symetric_extension),half(0,-75))+SU_RD_ter_AU_penalty))

SU_RD_mm = rbind(SU_RD_mm, data.frame(x=c("GU","UC"),y=c("CT","TG"),num=c(half(70,-45)+(SU_RD_internal_loop_initiation+SU_RD_internal_loop_symetric_extension),half(70,-45))))
SU_RD_mm = rbind(SU_RD_mm, data.frame(x=c("CU","UG"),y=c("GT","TC"),num=c(half(70,12)+(SU_RD_internal_loop_initiation+SU_RD_internal_loop_symetric_extension),half(70,12))))
SU_RD_mm = rbind(SU_RD_mm, data.frame(x=c("AU","UU"),y=c("TT","TA"),num=c(half(70,-71)+(SU_RD_internal_loop_initiation+SU_RD_internal_loop_symetric_extension),half(70,-71))+SU_RD_ter_AU_penalty))
SU_RD_mm = rbind(SU_RD_mm, data.frame(x=c("UU","UA"),y=c("AT","TT"),num=c(half(70,-68)+(SU_RD_internal_loop_initiation+SU_RD_internal_loop_symetric_extension),half(70,-68))+SU_RD_ter_AU_penalty))

SU_RD_mm = rbind(SU_RD_mm, data.frame(x=c("GG","UC"),y=c("CT","GG"),num=c(half(153,-8)+(SU_RD_internal_loop_initiation+SU_RD_internal_loop_symetric_extension),half(153,-8))))
SU_RD_mm = rbind(SU_RD_mm, data.frame(x=c("CG","UG"),y=c("GT","GC"),num=c(half(141,47)+(SU_RD_internal_loop_initiation+SU_RD_internal_loop_symetric_extension),half(141,47))))
SU_RD_mm = rbind(SU_RD_mm, data.frame(x=c("AG","UU"),y=c("TT","GA"),num=c(half(55,-71)+(SU_RD_internal_loop_initiation+SU_RD_internal_loop_symetric_extension),half(55,-71))+SU_RD_ter_AU_penalty))
SU_RD_mm = rbind(SU_RD_mm, data.frame(x=c("UG","UA"),y=c("AT","GT"),num=c(half(100,-43)+(SU_RD_internal_loop_initiation+SU_RD_internal_loop_symetric_extension),half(100,-43))+SU_RD_ter_AU_penalty))

SU_RD_mm = rbind(SU_RD_mm, data.frame(x=c("GU","GC"),y=c("CG","TG"),num=c(half(251,59)+(SU_RD_internal_loop_initiation+SU_RD_internal_loop_symetric_extension),half(251,59))))
SU_RD_mm = rbind(SU_RD_mm, data.frame(x=c("CU","GG"),y=c("GG","TC"),num=c(half(211,32)+(SU_RD_internal_loop_initiation+SU_RD_internal_loop_symetric_extension),half(211,32))))
SU_RD_mm = rbind(SU_RD_mm, data.frame(x=c("AU","GU"),y=c("TG","TA"),num=c(half(136,-7)+(SU_RD_internal_loop_initiation+SU_RD_internal_loop_symetric_extension),half(136,-7))+SU_RD_ter_AU_penalty))
SU_RD_mm = rbind(SU_RD_mm, data.frame(x=c("UU","GA"),y=c("AG","TT"),num=c(half(127,-34)+(SU_RD_internal_loop_initiation+SU_RD_internal_loop_symetric_extension),half(127,-34))+SU_RD_ter_AU_penalty))

SU_RD_mm = rbind(SU_RD_mm, data.frame(x="AN",y=paste("T",d_nts[d_nts!='-'],sep=""),num=min(SU_RD_mm$num[SU_RD_mm$x%in%paste("A",r_nts,sep="")&SU_RD_mm$y%in%paste("T",d_nts,sep="")])))
SU_RD_mm = rbind(SU_RD_mm, data.frame(x="CN",y=paste("G",d_nts[d_nts!='-'],sep=""),num=min(SU_RD_mm$num[SU_RD_mm$x%in%paste("C",r_nts,sep="")&SU_RD_mm$y%in%paste("G",d_nts,sep="")])))
SU_RD_mm = rbind(SU_RD_mm, data.frame(x="GN",y=paste("C",d_nts[d_nts!='-'],sep=""),num=min(SU_RD_mm$num[SU_RD_mm$x%in%paste("G",r_nts,sep="")&SU_RD_mm$y%in%paste("C",d_nts,sep="")])))
SU_RD_mm = rbind(SU_RD_mm, data.frame(x="UN",y=paste("A",d_nts[d_nts!='-'],sep=""),num=min(SU_RD_mm$num[SU_RD_mm$x%in%paste("U",r_nts,sep="")&SU_RD_mm$y%in%paste("A",d_nts,sep="")])))
SU_RD_mm = rbind(SU_RD_mm, data.frame(y="AN",x=paste("U",r_nts[!r_nts%in%c('-','N')],sep=""),num=min(SU_RD_mm$num[SU_RD_mm$y%in%paste("A",d_nts,sep="")&SU_RD_mm$x%in%paste("U",r_nts,sep="")])))
SU_RD_mm = rbind(SU_RD_mm, data.frame(y="CN",x=paste("G",r_nts[!r_nts%in%c('-','N')],sep=""),num=min(SU_RD_mm$num[SU_RD_mm$y%in%paste("C",d_nts,sep="")&SU_RD_mm$x%in%paste("G",r_nts,sep="")])))
SU_RD_mm = rbind(SU_RD_mm, data.frame(y="GN",x=paste("C",r_nts[!r_nts%in%c('-','N')],sep=""),num=min(SU_RD_mm$num[SU_RD_mm$y%in%paste("G",d_nts,sep="")&SU_RD_mm$x%in%paste("C",r_nts,sep="")])))
SU_RD_mm = rbind(SU_RD_mm, data.frame(y="TN",x=paste("A",r_nts[!r_nts%in%c('-','N')],sep=""),num=min(SU_RD_mm$num[SU_RD_mm$y%in%paste("T",d_nts,sep="")&SU_RD_mm$x%in%paste("A",r_nts,sep="")])))

SU_RD_mm = rbind(SU_RD_mm, data.frame(x="NA",y=paste(d_nts[d_nts!='-'],"T",sep=""),num=SU_RD_ter_AT_penalty))
SU_RD_mm = rbind(SU_RD_mm, data.frame(x="NC",y=paste(d_nts[d_nts!='-'],"G",sep=""),num=0))
SU_RD_mm = rbind(SU_RD_mm, data.frame(x="NG",y=paste(d_nts[d_nts!='-'],"C",sep=""),num=0))
SU_RD_mm = rbind(SU_RD_mm, data.frame(x="NU",y=paste(d_nts[d_nts!='-'],"A",sep=""),num=SU_RD_ter_AT_penalty))
SU_RD_mm = rbind(SU_RD_mm, data.frame(y="NA",x=paste(r_nts[!r_nts%in%c('-','N')],"U",sep=""),num=SU_RD_ter_AT_penalty))
SU_RD_mm = rbind(SU_RD_mm, data.frame(y="NC",x=paste(r_nts[!r_nts%in%c('-','N')],"G",sep=""),num=0))
SU_RD_mm = rbind(SU_RD_mm, data.frame(y="NG",x=paste(r_nts[!r_nts%in%c('-','N')],"C",sep=""),num=0))
SU_RD_mm = rbind(SU_RD_mm, data.frame(y="NT",x=paste(r_nts[!r_nts%in%c('-','N')],"A",sep=""),num=SU_RD_ter_AT_penalty))
##########################################################################


SU_RD_all = rbind(SU_RD_stack, SU_RD_init, SU_RD_ter_AT, rd_commons, SU_RD_bulges, SU_RD_loops, SU_RD_mm)
SU_RD_all$x=as.character(SU_RD_all$x)
SU_RD_all$y=as.character(SU_RD_all$y)

for (i in 0:length(d_nns)){
  if (i==0) {cat("    ",paste(d_nns,collapse = "   "))}
  for (j in 0:length(d_nns)){
    if (j==0) {cat(r_nns[i])}
    else{
      val = as.character(SU_RD_all$num[SU_RD_all$x==r_nns[i] & SU_RD_all$y==d_nns[j]])
      if (length(val)==0) {val=""}
      else if (is.na(val) | val=="NA") {val="*"}
      else if (val==-Inf) {val="-inf"}
      valc=nchar(val)
      cat(paste(rep(" ",(4-valc)), collapse=""),val)
    }
  }
  cat("\n")
}

######################################################################
################ SU RD with GU allowed ###############################
###### IN PROGRESS

SU_RD_wGU_stack = NULL
SU_RD_wGU_stack = rbind(SU_RD_wGU_stack, data.frame(x="AA",y="TT",num=100))
SU_RD_wGU_stack = rbind(SU_RD_wGU_stack, data.frame(x="AC",y="TG",num=210)) #
SU_RD_wGU_stack = rbind(SU_RD_wGU_stack, data.frame(x="AG",y="TC",num=180)) #
SU_RD_wGU_stack = rbind(SU_RD_wGU_stack, data.frame(x="AU",y="TA",num=90)) 
SU_RD_wGU_stack = rbind(SU_RD_wGU_stack, data.frame(x="CA",y="GT",num=90)) 
SU_RD_wGU_stack = rbind(SU_RD_wGU_stack, data.frame(x="CC",y="GG",num=210)) #
SU_RD_wGU_stack = rbind(SU_RD_wGU_stack, data.frame(x="CG",y="GC",num=170))
SU_RD_wGU_stack = rbind(SU_RD_wGU_stack, data.frame(x="CU",y="GA",num=90)) 
SU_RD_wGU_stack = rbind(SU_RD_wGU_stack, data.frame(x="GA",y="CT",num=130))
SU_RD_wGU_stack = rbind(SU_RD_wGU_stack, data.frame(x="GC",y="CG",num=270))
SU_RD_wGU_stack = rbind(SU_RD_wGU_stack, data.frame(x="GG",y="CC",num=290))
SU_RD_wGU_stack = rbind(SU_RD_wGU_stack, data.frame(x="GU",y="CA",num=110))
SU_RD_wGU_stack = rbind(SU_RD_wGU_stack, data.frame(x="UA",y="AT",num=60))
SU_RD_wGU_stack = rbind(SU_RD_wGU_stack, data.frame(x="UC",y="AG",num=150)) # 
SU_RD_wGU_stack = rbind(SU_RD_wGU_stack, data.frame(x="UG",y="AC",num=160)) #
SU_RD_wGU_stack = rbind(SU_RD_wGU_stack, data.frame(x="UU",y="AA",num=20)) #

SU_RD_wGU_stack = rbind(SU_RD_wGU_stack, data.frame(x=c("GG","UC"),y=c("CT","GG"),num=half(153,-8))) #GU
SU_RD_wGU_stack = rbind(SU_RD_wGU_stack, data.frame(x=c("CG","UG"),y=c("GT","GC"),num=half(141,47))) #GU
SU_RD_wGU_stack = rbind(SU_RD_wGU_stack, data.frame(x=c("AG","UU"),y=c("TT","GA"),num=half(55,-71))) #GU
SU_RD_wGU_stack = rbind(SU_RD_wGU_stack, data.frame(x=c("UG","UA"),y=c("AT","GT"),num=half(100,-43))) #GU
SU_RD_wGU_stack = rbind(SU_RD_wGU_stack, data.frame(x=c("GU","GC"),y=c("CG","TG"),num=half(251,59))) #GU
SU_RD_wGU_stack = rbind(SU_RD_wGU_stack, data.frame(x=c("CU","GG"),y=c("GG","TC"),num=half(211,32))) #GU
SU_RD_wGU_stack = rbind(SU_RD_wGU_stack, data.frame(x=c("AU","GU"),y=c("TG","TA"),num=half(136,-7))) #GU
SU_RD_wGU_stack = rbind(SU_RD_wGU_stack, data.frame(x=c("UU","GA"),y=c("AG","TT"),num=half(127,-34))) #GU
SU_RD_wGU_stack = rbind(SU_RD_wGU_stack, data.frame(x=c("GG","UU"),y=c("TT","GG"),num=half(50,-74))) #GU
SU_RD_wGU_stack = rbind(SU_RD_wGU_stack, data.frame(x="GU",y="TG",num=half(-129,-115))) #GU
SU_RD_wGU_stack = rbind(SU_RD_wGU_stack, data.frame(x="UG",y="GT",num=half(-30,-52))) #GU


SU_RD_wGU_ter_AT_penalty = round((SL_H_DD_ter_AT_penalty+RI_RR_ter_AU_penalty)/2)
SU_RD_wGU_ter_AU_penalty = round((SL_H_DD_ter_AT_penalty+RI_RR_ter_AU_penalty)/2)

SU_RD_wGU_ter_AT = NULL
SU_RD_wGU_ter_AT = rbind(SU_RD_wGU_ter_AT, data.frame(x="A-",y="T-",num=SU_RD_wGU_ter_AT_penalty))
SU_RD_wGU_ter_AT = rbind(SU_RD_wGU_ter_AT, data.frame(x="U-",y="A-",num=SU_RD_wGU_ter_AU_penalty))
SU_RD_wGU_ter_AT = rbind(SU_RD_wGU_ter_AT, data.frame(x="G-",y="T-",num=SU_RD_wGU_ter_AU_penalty))
SU_RD_wGU_ter_AT = rbind(SU_RD_wGU_ter_AT, data.frame(x="U-",y="G-",num=SU_RD_wGU_ter_AU_penalty))
#SU_RD_wGU_ter_AT = rbind(SU_RD_wGU_ter_AT, data.frame(x="U-",y="G-",num=SU_RD_wGU_ter_GU_penalty))
SU_RD_wGU_ter_AT = rbind(SU_RD_wGU_ter_AT, data.frame(x=paste(c("G","C"),"-",sep=""),y=paste(c("C","G"),"-",sep=""),num=0))
SU_RD_wGU_ter_AT = rbind(SU_RD_wGU_ter_AT, data.frame(x=paste(c("A","A","A","C","C","C","G","G","U","U","A","C","G","U","N","N","N","N"),"-",sep=""),
                                                      y=paste(c("A","C","G","A","C","T","A","G","C","T","N","N","N","N","A","C","G","T"),"-",sep=""),num=-Inf))

SU_RD_wGU_init_eng = -310
SU_RD_wGU_init = NULL
SU_RD_wGU_init = rbind(SU_RD_wGU_init, data.frame(x=paste("-",c("C","G"),sep=""),y=paste("-",c("G","C"),sep=""),num=(ris2_minus+SU_RD_wGU_init_eng )))
SU_RD_wGU_init = rbind(SU_RD_wGU_init, data.frame(x="-A",y="-T",num=(ris2_minus+SU_RD_wGU_init_eng+SU_RD_wGU_ter_AT_penalty)))
SU_RD_wGU_init = rbind(SU_RD_wGU_init, data.frame(x=c("-U","-G","-U"),y=c("-A","-T","-G"),num=(ris2_minus+SU_RD_wGU_init_eng+SU_RD_wGU_ter_AU_penalty)))
#SU_RD_wGU_init = rbind(SU_RD_wGU_init, data.frame(x="-U",y="-G",num=(ris2_minus+SU_RD_wGU_init_eng+SU_RD_wGU_ter_GU_penalty)))
SU_RD_wGU_init = rbind(SU_RD_wGU_init, data.frame(x=paste("-",c("A","A","A","C","C","C","G","G","U","U","A","C","G","U","N","N","N","N","N"),sep=""),
                                                  y=paste("-",c("A","C","G","A","C","T","A","G","C","T","N","N","N","N","A","C","G","T","N"),sep=""),num=0))

SU_RD_wGU_bulge_initiation = round((RI_RR_bulge_initiation+SL_H_DD_bulge_initiation)/2)
SU_RD_wGU_bulge_extension = round((RI_RR_bulge_extension+SL_H_DD_bulge_extension)/2)
SU_RD_wGU_bulges = NULL
SU_RD_wGU_bulges = rbind(SU_RD_wGU_bulges, data.frame(x="C-",y=c(paste("G",d_nts[d_nts!='-'],sep = "")),num=SU_RD_wGU_bulge_initiation+SU_RD_wGU_bulge_extension))
SU_RD_wGU_bulges = rbind(SU_RD_wGU_bulges, data.frame(x="G-",y=c(paste("C",d_nts[d_nts!='-'],sep = "")),num=SU_RD_wGU_bulge_initiation+SU_RD_wGU_bulge_extension))
SU_RD_wGU_bulges = rbind(SU_RD_wGU_bulges, data.frame(y="C-",x=c(paste("G",r_nts[r_nts!='-'],sep = "")),num=SU_RD_wGU_bulge_initiation+SU_RD_wGU_bulge_extension))
SU_RD_wGU_bulges = rbind(SU_RD_wGU_bulges, data.frame(y="G-",x=c(paste("C",r_nts[r_nts!='-'],sep = "")),num=SU_RD_wGU_bulge_initiation+SU_RD_wGU_bulge_extension))
SU_RD_wGU_bulges = rbind(SU_RD_wGU_bulges, data.frame(x="-C",y=c(paste(d_nts[d_nts!='-'],"G",sep = "")),num=0))
SU_RD_wGU_bulges = rbind(SU_RD_wGU_bulges, data.frame(x="-G",y=c(paste(d_nts[d_nts!='-'],"C",sep = "")),num=0))
SU_RD_wGU_bulges = rbind(SU_RD_wGU_bulges, data.frame(y="-C",x=c(paste(r_nts[r_nts!='-'],"G",sep = "")),num=0))
SU_RD_wGU_bulges = rbind(SU_RD_wGU_bulges, data.frame(y="-G",x=c(paste(r_nts[r_nts!='-'],"C",sep = "")),num=0))
SU_RD_wGU_bulges = rbind(SU_RD_wGU_bulges, data.frame(x="A-",y=c(paste("T",d_nts[d_nts!='-'],sep = "")),num=SU_RD_wGU_bulge_initiation+SU_RD_wGU_bulge_extension+SU_RD_wGU_ter_AT_penalty))
SU_RD_wGU_bulges = rbind(SU_RD_wGU_bulges, data.frame(x="U-",y=c(paste("A",d_nts[d_nts!='-'],sep = "")),num=SU_RD_wGU_bulge_initiation+SU_RD_wGU_bulge_extension+SU_RD_wGU_ter_AT_penalty))
SU_RD_wGU_bulges = rbind(SU_RD_wGU_bulges, data.frame(y="A-",x=c(paste("U",r_nts[r_nts!='-'],sep = "")),num=SU_RD_wGU_bulge_initiation+SU_RD_wGU_bulge_extension+SU_RD_wGU_ter_AT_penalty))
SU_RD_wGU_bulges = rbind(SU_RD_wGU_bulges, data.frame(y="T-",x=c(paste("A",r_nts[r_nts!='-'],sep = "")),num=SU_RD_wGU_bulge_initiation+SU_RD_wGU_bulge_extension+SU_RD_wGU_ter_AT_penalty))
SU_RD_wGU_bulges = rbind(SU_RD_wGU_bulges, data.frame(x="-A",y=c(paste(d_nts[d_nts!='-'],"T",sep = "")),num=SU_RD_wGU_ter_AT_penalty))
SU_RD_wGU_bulges = rbind(SU_RD_wGU_bulges, data.frame(x="-U",y=c(paste(d_nts[d_nts!='-'],"A",sep = "")),num=SU_RD_wGU_ter_AT_penalty))
SU_RD_wGU_bulges = rbind(SU_RD_wGU_bulges, data.frame(y="-A",x=c(paste(r_nts[r_nts!='-'],"U",sep = "")),num=SU_RD_wGU_ter_AT_penalty))
SU_RD_wGU_bulges = rbind(SU_RD_wGU_bulges, data.frame(y="-T",x=c(paste(r_nts[r_nts!='-'],"A",sep = "")),num=SU_RD_wGU_ter_AT_penalty))
# GUs
SU_RD_wGU_bulges = rbind(SU_RD_wGU_bulges, data.frame(x="G-",y=c(paste("T",d_nts[d_nts!='-'],sep = "")),num=SU_RD_wGU_bulge_initiation+SU_RD_wGU_bulge_extension+SU_RD_wGU_ter_AT_penalty))
SU_RD_wGU_bulges = rbind(SU_RD_wGU_bulges, data.frame(x="U-",y=c(paste("G",d_nts[d_nts!='-'],sep = "")),num=SU_RD_wGU_bulge_initiation+SU_RD_wGU_bulge_extension+SU_RD_wGU_ter_AT_penalty))
SU_RD_wGU_bulges = rbind(SU_RD_wGU_bulges, data.frame(y="G-",x=c(paste("U",r_nts[r_nts!='-'],sep = "")),num=SU_RD_wGU_bulge_initiation+SU_RD_wGU_bulge_extension+SU_RD_wGU_ter_AT_penalty))
SU_RD_wGU_bulges = rbind(SU_RD_wGU_bulges, data.frame(y="T-",x=c(paste("G",r_nts[r_nts!='-'],sep = "")),num=SU_RD_wGU_bulge_initiation+SU_RD_wGU_bulge_extension+SU_RD_wGU_ter_AT_penalty))
SU_RD_wGU_bulges = rbind(SU_RD_wGU_bulges, data.frame(x="-G",y=c(paste(d_nts[d_nts!='-'],"T",sep = "")),num=SU_RD_wGU_ter_AT_penalty))
SU_RD_wGU_bulges = rbind(SU_RD_wGU_bulges, data.frame(x="-U",y=c(paste(d_nts[d_nts!='-'],"G",sep = "")),num=SU_RD_wGU_ter_AT_penalty))
SU_RD_wGU_bulges = rbind(SU_RD_wGU_bulges, data.frame(y="-G",x=c(paste(r_nts[r_nts!='-'],"U",sep = "")),num=SU_RD_wGU_ter_AT_penalty))
SU_RD_wGU_bulges = rbind(SU_RD_wGU_bulges, data.frame(y="-T",x=c(paste(r_nts[r_nts!='-'],"G",sep = "")),num=SU_RD_wGU_ter_AT_penalty))
for (i in 1:5){
  SU_RD_wGU_bulges = rbind(SU_RD_wGU_bulges, data.frame(x=paste(r_nts[r_nts!='-'],r_nts[i],sep = ""),y="--",num=SU_RD_wGU_bulge_extension))
  SU_RD_wGU_bulges = rbind(SU_RD_wGU_bulges, data.frame(x="--",y=paste(d_nts[d_nts!='-'],d_nts[i],sep = ""),num=SU_RD_wGU_bulge_extension))
}

SU_RD_wGU_internal_loop_asymetry = round((RI_RR_internal_loop_asymetry+SL_H_DD_internal_loop_asymetry)/2)
SU_RD_wGU_internal_loop_initiation = round((RI_RR_internal_loop_initiation+SL_H_DD_internal_loop_initiation)/2)
SU_RD_wGU_internal_loop_symetric_extension = round((RI_RR_internal_loop_symetric_extension+SL_H_DD_internal_loop_symetric_extension)/2)
SU_RD_wGU_internal_loop_asymetric_extension = round((RI_RR_internal_loop_asymetric_extension+SL_H_DD_internal_loop_asymetric_extension)/2)
SU_RD_wGU_loops = NULL
for (i in 1:6){
  for (j in 1:6){
    mm = paste(d_nts[i],d_nts[j],sep = "")
    if (!(mm%in%c(d_matches,"--","GT","TG"))) {
      for (k in 1:5){
        for (l in 1:5){
          mm = paste(d_nts[k],d_nts[l],sep = "")
          if (!(mm%in%c(d_matches,"--","GT","TG"))) {
            SU_RD_wGU_loops = rbind(SU_RD_wGU_loops, data.frame(x=paste(r_nts[i],r_nts[k],sep=""),y=paste(d_nts[j],d_nts[l],sep=""),num=SU_RD_wGU_internal_loop_symetric_extension))
          }
        }
      }
    }
  }
}

for (i in 1:5){
  for (j in 1:5){
    mm = paste(d_nts[i],d_nts[j],sep = "")
    if (!(mm%in%c(d_matches,"--","GT","TG"))) {
      for (k in 1:5){
        SU_RD_wGU_loops = rbind(SU_RD_wGU_loops, data.frame(x=paste(r_nts[i],'-',sep=""),y=paste(d_nts[j],d_nts[k],sep=""),
                                                            num=SU_RD_wGU_internal_loop_asymetric_extension+SU_RD_wGU_internal_loop_asymetry))
        SU_RD_wGU_loops = rbind(SU_RD_wGU_loops, data.frame(y=paste(d_nts[i],'-',sep=""),x=paste(r_nts[j],r_nts[k],sep=""),
                                                            num=SU_RD_wGU_internal_loop_asymetric_extension+SU_RD_wGU_internal_loop_asymetry))
      }
    }
  }
}


##########################################################################
## RIS2 RNA RNA mismatches unknown, only AG GA GG UU => 80 100 120 70 rest accepted as 0, and GU UG cases from *Turner (8 different energies)
SU_RD_wGU_mm = NULL
SU_RD_wGU_mm = rbind(SU_RD_wGU_mm, data.frame(x=c("GA","AC"),y=c("CA","AG"),num=c(half(0,-17)+(SU_RD_wGU_internal_loop_initiation+SU_RD_wGU_internal_loop_symetric_extension),half(0,-17))))
SU_RD_wGU_mm = rbind(SU_RD_wGU_mm, data.frame(x=c("CA","AG"),y=c("GA","AC"),num=c(half(0,-43)+(SU_RD_wGU_internal_loop_initiation+SU_RD_wGU_internal_loop_symetric_extension),half(0,-43))))
SU_RD_wGU_mm = rbind(SU_RD_wGU_mm, data.frame(x=c("AA","AU"),y=c("TA","AA"),num=c(half(0,-61)+(SU_RD_wGU_internal_loop_initiation+SU_RD_wGU_internal_loop_symetric_extension),half(0,-61))+SU_RD_wGU_ter_AU_penalty))
SU_RD_wGU_mm = rbind(SU_RD_wGU_mm, data.frame(x=c("UA","AA"),y=c("AA","AT"),num=c(half(0,-71)+(SU_RD_wGU_internal_loop_initiation+SU_RD_wGU_internal_loop_symetric_extension),half(0,-71))+SU_RD_wGU_ter_AU_penalty))
SU_RD_wGU_mm = rbind(SU_RD_wGU_mm, data.frame(x=c("GA","AU","UA","AG"),y=c("TA","AG","GA","AT"),
                                              num=c(half(0,SU_RD_wGU_internal_loop_symetric_extension)+(SU_RD_wGU_internal_loop_initiation+SU_RD_wGU_internal_loop_symetric_extension),
                                                    half(0,SU_RD_wGU_internal_loop_symetric_extension))+SU_RD_wGU_ter_AU_penalty))

SU_RD_wGU_mm = rbind(SU_RD_wGU_mm, data.frame(x=c("GA","CC"),y=c("CC","AG"),num=c(half(0,-81)+(SU_RD_wGU_internal_loop_initiation+SU_RD_wGU_internal_loop_symetric_extension),half(0,-81))))
SU_RD_wGU_mm = rbind(SU_RD_wGU_mm, data.frame(x=c("CA","CG"),y=c("GC","AC"),num=c(half(0,-75)+(SU_RD_wGU_internal_loop_initiation+SU_RD_wGU_internal_loop_symetric_extension),half(0,-75))))
SU_RD_wGU_mm = rbind(SU_RD_wGU_mm, data.frame(x=c("AA","CU"),y=c("TC","AA"),num=c(half(0,-88)+(SU_RD_wGU_internal_loop_initiation+SU_RD_wGU_internal_loop_symetric_extension),half(0,-88))+SU_RD_wGU_ter_AU_penalty))
SU_RD_wGU_mm = rbind(SU_RD_wGU_mm, data.frame(x=c("UA","CA"),y=c("AC","AT"),num=c(half(0,-92)+(SU_RD_wGU_internal_loop_initiation+SU_RD_wGU_internal_loop_symetric_extension),half(0,-92))+SU_RD_wGU_ter_AU_penalty))
SU_RD_wGU_mm = rbind(SU_RD_wGU_mm, data.frame(x=c("GA","CU","UA","CG"),y=c("TC","AG","GC","AT"),
                                              num=c(half(0,SU_RD_wGU_internal_loop_symetric_extension)+(SU_RD_wGU_internal_loop_initiation+SU_RD_wGU_internal_loop_symetric_extension),
                                                    half(0,SU_RD_wGU_internal_loop_symetric_extension))+SU_RD_wGU_ter_AU_penalty))

#AG terminal mismatch specialty from R-R
SU_RD_wGU_mm = rbind(SU_RD_wGU_mm, data.frame(x=c("GA","GC"),y=c("CG","AG"),num=c(half(80, 25)+(SU_RD_wGU_internal_loop_initiation+SU_RD_wGU_internal_loop_symetric_extension),half(80, 25))))
SU_RD_wGU_mm = rbind(SU_RD_wGU_mm, data.frame(x=c("CA","GG"),y=c("GG","AC"),num=c(half(80, -3)+(SU_RD_wGU_internal_loop_initiation+SU_RD_wGU_internal_loop_symetric_extension),half(80, -3))))
SU_RD_wGU_mm = rbind(SU_RD_wGU_mm, data.frame(x=c("AA","GU"),y=c("TG","AA"),num=c(half(80,-14)+(SU_RD_wGU_internal_loop_initiation+SU_RD_wGU_internal_loop_symetric_extension),half(80,-14))+SU_RD_wGU_ter_AU_penalty))
SU_RD_wGU_mm = rbind(SU_RD_wGU_mm, data.frame(x=c("UA","GA"),y=c("AG","AT"),num=c(half(80,-42)+(SU_RD_wGU_internal_loop_initiation+SU_RD_wGU_internal_loop_symetric_extension),half(80,-42))+SU_RD_wGU_ter_AU_penalty))
SU_RD_wGU_mm = rbind(SU_RD_wGU_mm,data.frame(x=c("GA","GU","UA","GG"),y=c("TG","AG","GG","AT"),
                                             num=c(half(80,SU_RD_wGU_internal_loop_symetric_extension)+(SU_RD_wGU_internal_loop_initiation+SU_RD_wGU_internal_loop_symetric_extension),
                                                   half(80,SU_RD_wGU_internal_loop_symetric_extension))+SU_RD_wGU_ter_AU_penalty))

SU_RD_wGU_mm = rbind(SU_RD_wGU_mm, data.frame(x=c("GC","AC"),y=c("CA","CG"),num=c(half(0,-47)+(SU_RD_wGU_internal_loop_initiation+SU_RD_wGU_internal_loop_symetric_extension),half(0,-47))))
SU_RD_wGU_mm = rbind(SU_RD_wGU_mm, data.frame(x=c("CC","AG"),y=c("GA","CC"),num=c(half(0,-79)+(SU_RD_wGU_internal_loop_initiation+SU_RD_wGU_internal_loop_symetric_extension),half(0,-79))))
SU_RD_wGU_mm = rbind(SU_RD_wGU_mm, data.frame(x=c("AC","AU"),y=c("TA","CA"),num=c(half(0,-77)+(SU_RD_wGU_internal_loop_initiation+SU_RD_wGU_internal_loop_symetric_extension),half(0,-77))+SU_RD_wGU_ter_AU_penalty))
SU_RD_wGU_mm = rbind(SU_RD_wGU_mm, data.frame(x=c("UC","AA"),y=c("AA","CT"),num=c(half(0,-133)+(SU_RD_wGU_internal_loop_initiation+SU_RD_wGU_internal_loop_symetric_extension),half(0,-133))+SU_RD_wGU_ter_AU_penalty))
SU_RD_wGU_mm = rbind(SU_RD_wGU_mm,data.frame(x=c("GC","AU","UC","AG"),y=c("TA","CG","GA","CT"),
                                             num=c(half(0,SU_RD_wGU_internal_loop_symetric_extension)+(SU_RD_wGU_internal_loop_initiation+SU_RD_wGU_internal_loop_symetric_extension),
                                                   half(0,SU_RD_wGU_internal_loop_symetric_extension))+SU_RD_wGU_ter_AU_penalty))

SU_RD_wGU_mm = rbind(SU_RD_wGU_mm, data.frame(x=c("GC","CC"),y=c("CC","CG"),num=c(half(0, -79)+(SU_RD_wGU_internal_loop_initiation+SU_RD_wGU_internal_loop_symetric_extension),half(0, -79))))
SU_RD_wGU_mm = rbind(SU_RD_wGU_mm, data.frame(x=c("CC","CG"),y=c("GC","CC"),num=c(half(0, -70)+(SU_RD_wGU_internal_loop_initiation+SU_RD_wGU_internal_loop_symetric_extension),half(0, -70))))
SU_RD_wGU_mm = rbind(SU_RD_wGU_mm, data.frame(x=c("AC","CU"),y=c("TC","CA"),num=c(half(0,-133)+(SU_RD_wGU_internal_loop_initiation+SU_RD_wGU_internal_loop_symetric_extension),half(0,-133))+SU_RD_wGU_ter_AU_penalty))
SU_RD_wGU_mm = rbind(SU_RD_wGU_mm, data.frame(x=c("UC","CA"),y=c("AC","CT"),num=c(half(0,-105)+(SU_RD_wGU_internal_loop_initiation+SU_RD_wGU_internal_loop_symetric_extension),half(0,-105))+SU_RD_wGU_ter_AU_penalty))
SU_RD_wGU_mm = rbind(SU_RD_wGU_mm,data.frame(x=c("GC","CU","UC","CG"),y=c("TC","CG","GC","CT"),
                                             num=c(half(0,SU_RD_wGU_internal_loop_symetric_extension)+(SU_RD_wGU_internal_loop_initiation+SU_RD_wGU_internal_loop_symetric_extension),
                                                   half(0,SU_RD_wGU_internal_loop_symetric_extension))+SU_RD_wGU_ter_AU_penalty))

SU_RD_wGU_mm = rbind(SU_RD_wGU_mm, data.frame(x=c("GC","UC"),y=c("CT","CG"),num=c(half(0,-62)+(SU_RD_wGU_internal_loop_initiation+SU_RD_wGU_internal_loop_symetric_extension),half(0,-62))))
SU_RD_wGU_mm = rbind(SU_RD_wGU_mm, data.frame(x=c("CC","UG"),y=c("GT","CC"),num=c(half(0,-62)+(SU_RD_wGU_internal_loop_initiation+SU_RD_wGU_internal_loop_symetric_extension),half(0,-62))))
SU_RD_wGU_mm = rbind(SU_RD_wGU_mm, data.frame(x=c("AC","UU"),y=c("TT","CA"),num=c(half(0,-64)+(SU_RD_wGU_internal_loop_initiation+SU_RD_wGU_internal_loop_symetric_extension),half(0,-64))+SU_RD_wGU_ter_AU_penalty))
SU_RD_wGU_mm = rbind(SU_RD_wGU_mm, data.frame(x=c("UC","UA"),y=c("AT","CT"),num=c(half(0,-97)+(SU_RD_wGU_internal_loop_initiation+SU_RD_wGU_internal_loop_symetric_extension),half(0,-97))+SU_RD_wGU_ter_AU_penalty))
SU_RD_wGU_mm = rbind(SU_RD_wGU_mm,data.frame(x=c("GC","UU","UC","UG"),y=c("TT","CG","GT","CT"),
                                             num=c(half(0,SU_RD_wGU_internal_loop_symetric_extension)+(SU_RD_wGU_internal_loop_initiation+SU_RD_wGU_internal_loop_symetric_extension),
                                                   half(0,SU_RD_wGU_internal_loop_symetric_extension))+SU_RD_wGU_ter_AU_penalty))

#GA terminal mismatch specialty from R-R
SU_RD_wGU_mm = rbind(SU_RD_wGU_mm, data.frame(x=c("GG","AC"),y=c("CA","GG"),num=c(half(100, 52)+(SU_RD_wGU_internal_loop_initiation+SU_RD_wGU_internal_loop_symetric_extension),half(100,52)))) 
SU_RD_wGU_mm = rbind(SU_RD_wGU_mm, data.frame(x=c("CG","AG"),y=c("GA","GC"),num=c(half(100,-11)+(SU_RD_wGU_internal_loop_initiation+SU_RD_wGU_internal_loop_symetric_extension),half(100,-11))))
SU_RD_wGU_mm = rbind(SU_RD_wGU_mm, data.frame(x=c("AG","AU"),y=c("TA","GA"),num=c(half(100, -2)+(SU_RD_wGU_internal_loop_initiation+SU_RD_wGU_internal_loop_symetric_extension),half(100,-2))+SU_RD_wGU_ter_AU_penalty))
SU_RD_wGU_mm = rbind(SU_RD_wGU_mm, data.frame(x=c("UG","AA"),y=c("AA","GT"),num=c(half(100,-74)+(SU_RD_wGU_internal_loop_initiation+SU_RD_wGU_internal_loop_symetric_extension),half(100,-74))+SU_RD_wGU_ter_AU_penalty))
SU_RD_wGU_mm = rbind(SU_RD_wGU_mm,data.frame(x=c("GG","AU","UG","AG"),y=c("TA","GG","GA","GT"),
                                             num=c(half(100,SU_RD_wGU_internal_loop_symetric_extension)+(SU_RD_wGU_internal_loop_initiation+SU_RD_wGU_internal_loop_symetric_extension),
                                                   half(100,SU_RD_wGU_internal_loop_symetric_extension))+SU_RD_wGU_ter_AU_penalty))

#GG terminal mismatch specialty from R-R
SU_RD_wGU_mm = rbind(SU_RD_wGU_mm, data.frame(x=c("GG","GC"),y=c("CG","GG"),num=c(half(120,111)+(SU_RD_wGU_internal_loop_initiation+SU_RD_wGU_internal_loop_symetric_extension),half(120,111)))) 
SU_RD_wGU_mm = rbind(SU_RD_wGU_mm, data.frame(x=c("CG","GG"),y=c("GG","GC"),num=c(half(120,11)+(SU_RD_wGU_internal_loop_initiation+SU_RD_wGU_internal_loop_symetric_extension),half(120,11))))
SU_RD_wGU_mm = rbind(SU_RD_wGU_mm, data.frame(x=c("AG","GU"),y=c("TG","GA"),num=c(half(120,13)+(SU_RD_wGU_internal_loop_initiation+SU_RD_wGU_internal_loop_symetric_extension),half(120,13))+SU_RD_wGU_ter_AU_penalty))
SU_RD_wGU_mm = rbind(SU_RD_wGU_mm, data.frame(x=c("UG","GA"),y=c("AG","GT"),num=c(half(120,-44)+(SU_RD_wGU_internal_loop_initiation+SU_RD_wGU_internal_loop_symetric_extension),half(120,-44))+SU_RD_wGU_ter_AU_penalty))
SU_RD_wGU_mm = rbind(SU_RD_wGU_mm,data.frame(x=c("GG","GU","UG","GG"),y=c("TG","GG","GG","GT"),
                                             num=c(half(120,SU_RD_wGU_internal_loop_symetric_extension)+(SU_RD_wGU_internal_loop_initiation+SU_RD_wGU_internal_loop_symetric_extension),
                                                   half(120,SU_RD_wGU_internal_loop_symetric_extension))+SU_RD_wGU_ter_AU_penalty))

SU_RD_wGU_mm = rbind(SU_RD_wGU_mm, data.frame(x=c("GU","CC"),y=c("CC","TG"),num=c(half(0,-98)+(SU_RD_wGU_internal_loop_initiation+SU_RD_wGU_internal_loop_symetric_extension),half(0,-98))))
SU_RD_wGU_mm = rbind(SU_RD_wGU_mm, data.frame(x=c("CU","CG"),y=c("GC","TC"),num=c(half(0,-40)+(SU_RD_wGU_internal_loop_initiation+SU_RD_wGU_internal_loop_symetric_extension),half(0,-40))))
SU_RD_wGU_mm = rbind(SU_RD_wGU_mm, data.frame(x=c("AU","CU"),y=c("TC","TA"),num=c(half(0,-73)+(SU_RD_wGU_internal_loop_initiation+SU_RD_wGU_internal_loop_symetric_extension),half(0,-73))+SU_RD_wGU_ter_AU_penalty))
SU_RD_wGU_mm = rbind(SU_RD_wGU_mm, data.frame(x=c("UU","CA"),y=c("AC","TT"),num=c(half(0,-75)+(SU_RD_wGU_internal_loop_initiation+SU_RD_wGU_internal_loop_symetric_extension),half(0,-75))+SU_RD_wGU_ter_AU_penalty))
SU_RD_wGU_mm = rbind(SU_RD_wGU_mm,data.frame(x=c("GU","CU","UU","CG"),y=c("TC","TG","GC","TT"),
                                             num=c(half(0,SU_RD_wGU_internal_loop_symetric_extension)+(SU_RD_wGU_internal_loop_initiation+SU_RD_wGU_internal_loop_symetric_extension),
                                                   half(0,SU_RD_wGU_internal_loop_symetric_extension))+SU_RD_wGU_ter_AU_penalty))

SU_RD_wGU_mm = rbind(SU_RD_wGU_mm, data.frame(x=c("GU","UC"),y=c("CT","TG"),num=c(half(70,-45)+(SU_RD_wGU_internal_loop_initiation+SU_RD_wGU_internal_loop_symetric_extension),half(70,-45))))
SU_RD_wGU_mm = rbind(SU_RD_wGU_mm, data.frame(x=c("CU","UG"),y=c("GT","TC"),num=c(half(70,12)+(SU_RD_wGU_internal_loop_initiation+SU_RD_wGU_internal_loop_symetric_extension),half(70,12))))
SU_RD_wGU_mm = rbind(SU_RD_wGU_mm, data.frame(x=c("AU","UU"),y=c("TT","TA"),num=c(half(70,-71)+(SU_RD_wGU_internal_loop_initiation+SU_RD_wGU_internal_loop_symetric_extension),half(70,-71))+SU_RD_wGU_ter_AU_penalty))
SU_RD_wGU_mm = rbind(SU_RD_wGU_mm, data.frame(x=c("UU","UA"),y=c("AT","TT"),num=c(half(70,-68)+(SU_RD_wGU_internal_loop_initiation+SU_RD_wGU_internal_loop_symetric_extension),half(70,-68))+SU_RD_wGU_ter_AU_penalty))
SU_RD_wGU_mm = rbind(SU_RD_wGU_mm,data.frame(x=c("GU","UU","UU","UG"),y=c("TT","TG","GT","TT"),
                                             num=c(half(0,SU_RD_wGU_internal_loop_symetric_extension)+(SU_RD_wGU_internal_loop_initiation+SU_RD_wGU_internal_loop_symetric_extension),
                                                   half(0,SU_RD_wGU_internal_loop_symetric_extension))+SU_RD_wGU_ter_AU_penalty))

SU_RD_wGU_mm = rbind(SU_RD_wGU_mm, data.frame(x="AN",y=paste("T",d_nts[d_nts!='-'],sep=""),num=min(SU_RD_wGU_mm$num[SU_RD_wGU_mm$x%in%paste("A",r_nts,sep="")&SU_RD_wGU_mm$y%in%paste("T",d_nts,sep="")])))
SU_RD_wGU_mm = rbind(SU_RD_wGU_mm, data.frame(x="CN",y=paste("G",d_nts[d_nts!='-'],sep=""),num=min(SU_RD_wGU_mm$num[SU_RD_wGU_mm$x%in%paste("C",r_nts,sep="")&SU_RD_wGU_mm$y%in%paste("G",d_nts,sep="")])))
SU_RD_wGU_mm = rbind(SU_RD_wGU_mm, data.frame(x="GN",y=paste("C",d_nts[d_nts!='-'],sep=""),num=min(SU_RD_wGU_mm$num[SU_RD_wGU_mm$x%in%paste("G",r_nts,sep="")&SU_RD_wGU_mm$y%in%paste("C",d_nts,sep="")])))
SU_RD_wGU_mm = rbind(SU_RD_wGU_mm, data.frame(x="UN",y=paste("A",d_nts[d_nts!='-'],sep=""),num=min(SU_RD_wGU_mm$num[SU_RD_wGU_mm$x%in%paste("U",r_nts,sep="")&SU_RD_wGU_mm$y%in%paste("A",d_nts,sep="")])))
SU_RD_wGU_mm = rbind(SU_RD_wGU_mm, data.frame(x="GN",y=paste("T",d_nts[d_nts!='-'],sep=""),num=min(SU_RD_wGU_mm$num[SU_RD_wGU_mm$x%in%paste("G",r_nts,sep="")&SU_RD_wGU_mm$y%in%paste("C",d_nts,sep="")])))
SU_RD_wGU_mm = rbind(SU_RD_wGU_mm, data.frame(x="UN",y=paste("G",d_nts[d_nts!='-'],sep=""),num=min(SU_RD_wGU_mm$num[SU_RD_wGU_mm$x%in%paste("U",r_nts,sep="")&SU_RD_wGU_mm$y%in%paste("A",d_nts,sep="")])))

SU_RD_wGU_mm = rbind(SU_RD_wGU_mm, data.frame(y="AN",x=paste("U",r_nts[!r_nts%in%c('-','N')],sep=""),num=min(SU_RD_wGU_mm$num[SU_RD_wGU_mm$y%in%paste("A",d_nts,sep="")&SU_RD_wGU_mm$x%in%paste("U",r_nts,sep="")])))
SU_RD_wGU_mm = rbind(SU_RD_wGU_mm, data.frame(y="CN",x=paste("G",r_nts[!r_nts%in%c('-','N')],sep=""),num=min(SU_RD_wGU_mm$num[SU_RD_wGU_mm$y%in%paste("C",d_nts,sep="")&SU_RD_wGU_mm$x%in%paste("G",r_nts,sep="")])))
SU_RD_wGU_mm = rbind(SU_RD_wGU_mm, data.frame(y="GN",x=paste("C",r_nts[!r_nts%in%c('-','N')],sep=""),num=min(SU_RD_wGU_mm$num[SU_RD_wGU_mm$y%in%paste("G",d_nts,sep="")&SU_RD_wGU_mm$x%in%paste("C",r_nts,sep="")])))
SU_RD_wGU_mm = rbind(SU_RD_wGU_mm, data.frame(y="TN",x=paste("A",r_nts[!r_nts%in%c('-','N')],sep=""),num=min(SU_RD_wGU_mm$num[SU_RD_wGU_mm$y%in%paste("T",d_nts,sep="")&SU_RD_wGU_mm$x%in%paste("A",r_nts,sep="")])))
SU_RD_wGU_mm = rbind(SU_RD_wGU_mm, data.frame(y="GN",x=paste("U",r_nts[!r_nts%in%c('-','N')],sep=""),num=min(SU_RD_wGU_mm$num[SU_RD_wGU_mm$y%in%paste("G",d_nts,sep="")&SU_RD_wGU_mm$x%in%paste("C",r_nts,sep="")])))
SU_RD_wGU_mm = rbind(SU_RD_wGU_mm, data.frame(y="TN",x=paste("G",r_nts[!r_nts%in%c('-','N')],sep=""),num=min(SU_RD_wGU_mm$num[SU_RD_wGU_mm$y%in%paste("T",d_nts,sep="")&SU_RD_wGU_mm$x%in%paste("A",r_nts,sep="")])))

SU_RD_wGU_mm = rbind(SU_RD_wGU_mm, data.frame(x="NC",y=paste(d_nts[d_nts!='-'],"G",sep=""),num=0))
SU_RD_wGU_mm = rbind(SU_RD_wGU_mm, data.frame(x="NG",y=paste(d_nts[d_nts!='-'],"C",sep=""),num=0))
SU_RD_wGU_mm = rbind(SU_RD_wGU_mm, data.frame(y="NC",x=paste(r_nts[!r_nts%in%c('-','N')],"G",sep=""),num=0))
SU_RD_wGU_mm = rbind(SU_RD_wGU_mm, data.frame(y="NG",x=paste(r_nts[!r_nts%in%c('-','N')],"C",sep=""),num=0))
SU_RD_wGU_mm = rbind(SU_RD_wGU_mm, data.frame(x="NA",y=paste(d_nts[d_nts!='-'],"T",sep=""),num=SU_RD_wGU_ter_AT_penalty))
SU_RD_wGU_mm = rbind(SU_RD_wGU_mm, data.frame(x="NU",y=paste(d_nts[d_nts!='-'],"A",sep=""),num=SU_RD_wGU_ter_AT_penalty))
SU_RD_wGU_mm = rbind(SU_RD_wGU_mm, data.frame(y="NA",x=paste(r_nts[!r_nts%in%c('-','N')],"U",sep=""),num=SU_RD_wGU_ter_AT_penalty))
SU_RD_wGU_mm = rbind(SU_RD_wGU_mm, data.frame(y="NT",x=paste(r_nts[!r_nts%in%c('-','N')],"A",sep=""),num=SU_RD_wGU_ter_AT_penalty))
SU_RD_wGU_mm = rbind(SU_RD_wGU_mm, data.frame(x="NG",y=paste(d_nts[d_nts!='-'],"T",sep=""),num=SU_RD_wGU_ter_AT_penalty))
SU_RD_wGU_mm = rbind(SU_RD_wGU_mm, data.frame(x="NU",y=paste(d_nts[d_nts!='-'],"G",sep=""),num=SU_RD_wGU_ter_AT_penalty))
SU_RD_wGU_mm = rbind(SU_RD_wGU_mm, data.frame(y="NG",x=paste(r_nts[!r_nts%in%c('-','N')],"U",sep=""),num=SU_RD_wGU_ter_AT_penalty))
SU_RD_wGU_mm = rbind(SU_RD_wGU_mm, data.frame(y="NT",x=paste(r_nts[!r_nts%in%c('-','N')],"G",sep=""),num=SU_RD_wGU_ter_AT_penalty))
##########################################################################


SU_RD_wGU_all = rbind(SU_RD_wGU_stack, SU_RD_wGU_init, SU_RD_wGU_ter_AT, rd_commons, SU_RD_wGU_bulges, SU_RD_wGU_loops, SU_RD_wGU_mm)
SU_RD_wGU_all$x=as.character(SU_RD_wGU_all$x)
SU_RD_wGU_all$y=as.character(SU_RD_wGU_all$y)

for (i in 0:length(d_nns)){
  if (i==0) {cat("    ",paste(d_nns,collapse = "   "))}
  for (j in 0:length(d_nns)){
    if (j==0) {cat(r_nns[i])}
    else{
      val = as.character(SU_RD_wGU_all$num[SU_RD_wGU_all$x==r_nns[i] & SU_RD_wGU_all$y==d_nns[j]])
      if (length(val)==0) {val=""}
      else if (is.na(val) | val=="NA") {val="*"}
      else if (val==-Inf) {val="-inf"}
      valc=nchar(val)
      cat(paste(rep(" ",(4-valc)), collapse=""),val)
    }
  }
  cat("\n")
}


##################################################
############## TRANSFORMATIONS ###################
transform_into_code <-function(eng_df,nns_1,nns_2,instance_name){
  cat(instance_name)
  cat(" = {\\\n")
  for(i in 1:6){
    cat("{\\\n")
    for(j in 1:6){
      cat("{")
      for(k in 1:6){
        if (k!=1)
          cat(", ")
        cat("{")
        for(l in 1:6){
          if (l!=1)
            cat(", ")
          if (is.na(eng_df$num[eng_df$x==nns_1[((i-1)*6)+j] & eng_df$y==nns_2[((k-1)*6)+l]]))
            cat("-123")
          else if (is.infinite(eng_df$num[eng_df$x==nns_1[((i-1)*6)+j] & eng_df$y==nns_2[((k-1)*6)+l]]))
            cat("-2000")
          else
            cat(as.character(eng_df$num[eng_df$x==nns_1[((i-1)*6)+j] & eng_df$y==nns_2[((k-1)*6)+l]]))
        }
        cat("}")
      } 
      if(j==6)
        cat(" }\\\n")
      else
        cat(" },\\\n")
    }
    if(i==6)
      cat("}\\\n};\n")
    else
      cat("},\\\n")
  }
}

transform_into_code(RI_RR_all,r_nns_code,r_nns_code,"dsm_t04_pos")

transform_into_code(SL_H_DD_all,d_nns_code,d_nns_code,"dsm_slh04_woGU_pos")

transform_into_code(SU_RD_all,r_nns_code,d_nns_code,"dsm_su95_rev_woGU_pos")

transform_into_code(SU_RD_wGU_all,r_nns_code,d_nns_code,"dsm_su95_rev_wGU_pos")
