library(openxlsx)
library(stringr)
library(caTools)
library(preprocessCore)

FL = c("Plasma", "Saliva")
NNN = length(FL)

for (nnn in 1:NNN ) {

   print (nnn)
   
   file = FL[nnn]
  
   IFN = paste( file, "_SPIA.txt", sep = "")
   data = read.table(IFN, header = TRUE, sep = "\t")
   
   exp = data[,-1]
   
   SYMBOL = as.vector(data[,1])
   
   SN = colnames(exp)
   
   NS = length(SN)
   
   indMM = c()
   indHL = c()
   
   for (ns in 1:NS) {
   
      record = SN[ns]
      
      rec = str_sub(record, start = 1, end = 2)
      
      if ( rec == "MM" ) indMM = c(indMM, ns)
      if ( rec == "HL" ) indHL = c(indHL, ns)
   
   }
   
   NMM = length(indMM)
   NHL = length(indHL)
   
   ind = c(indMM,indHL)

   exp = exp[,ind]
   
   ScoreMM = array(1,dim = c(NMM))
   ScoreHL = array(0,dim = c(NHL))
   
   Score = c(ScoreMM,ScoreHL)
   
   NG = length(SYMBOL)
   
   SN = colnames(exp)
   
   exp0 = as.matrix(exp)
   
   NRR = nrow(exp0)
   
   NCC = ncol(exp0)
   
   exp = mapply(exp0, FUN = as.numeric)
   
   exp <- matrix(data=exp, ncol=NCC, nrow=NRR)
   
   SN = colnames(exp0)
   
   NS = length(SN)
   
   print (NS)
   
   for (ns in 1:NS) {
   
       print (ns)
       
       e_x_p = exp[,-ns]
       s_c_o_r_e = Score [-ns]
       S_N = SN[-ns]
       
       AUC = array(100500,dim = c (NG))
           
       #print (NG) 
       
       for ( ng in 1:NG ) {
           ggg = as.numeric(as.vector(t(e_x_p[ng,])))
           AUC[ng] = colAUC(ggg,s_c_o_r_e) 
       }
       
       #print ("maxAUC")
       #print (max(AUC))
      
       ORDER = order(AUC)
       REVORDER = ORDER

       for ( ng in 1:NG ) {
           REVORDER[ng] = ORDER[NG+1-ng] 
       }

       AUC_decr = AUC[REVORDER]

       MGI = REVORDER[1:30]

       GENES = SYMBOL[MGI]
     
       if ( ns == 1 ) GENES_all = GENES
       if ( ns > 1 ) GENES_all = intersect(GENES_all,GENES)

       Texp = t(e_x_p[MGI,])
       SN_ = as.vector(S_N)
       Flag = as.vector(t(s_c_o_r_e))
       T = cbind (S_N,Flag,Texp)
     
       CNcur = c("Name", "Score", GENES)
       
       OFN = paste ("spia/", SN[ns], ".csv", sep = "" )
       write.table(T, OFN, row.names = FALSE, col.names = CNcur, sep = ",")
       
       print (length(GENES_all))
       
   }
   
   NCG = length(GENES_all)
   
   indCore = c()
   
   for (ncg in 1:NCG) {
   
       gene = GENES_all[ncg]
       
       indC = which( SYMBOL == gene )
   
       indCore = c(indCore, indC)
   
   }
   
   SYMBOL = SYMBOL[indCore]
   
   exp = exp[indCore,]
   
   SN_ = as.vector(SN)
   Flag = as.vector(t(Score))
   Texp = t(exp)
       
   T = cbind (SN,Flag,Texp)
     
   CNcur = c("Name", "Score", GENES_all)
 
   OFN = paste (file, "_robust_spia.csv", sep = "" )
   write.table(T, OFN, row.names = FALSE, col.names = CNcur, sep = ",")
   
}


