library(dplyr)
library(data.table)
library(stringr)
options(scipen=5)
data <- read.csv("nuclear-wallet-cards.txt", col.names = "strData", header = FALSE, sep = "\n", 
                 stringsAsFactors = FALSE )

decimalnumcount<-function(x){
  if(!grepl("\\d", x) | class(x)!="character")
    return(0)
  exp_num <- 0
  if(grepl("[E|e][+0-9]+", x)){
    exp_num <- (str_extract(x, "[E|e][+0-9]+")%>%str_extract(., "[+0-9]+") )%>%as.numeric()
    x <- sub("[E|e][+0-9]+", "", x)
  }
  x <- as.numeric(x)
  while((x %% 1) != 0){
    exp_num <- exp_num - 1
    x <- x * 10
  }
  return(exp_num)
}


cardName <- list("A" = 2: 4,
                 "Isomeric" = 5,
                 "Z" = 7: 9,
                 "Element" = 11: 12,
                 "Jp" = 17:26,
                 "DecayMode" = 31: 34,
                 "BranchRatio" = 35: 41,
                 "Ex_MeV" = 43: 49,
                 "Q_MeV" = 50: 60,
                 "HalfLife" = 64: 80,
                 "Abundance" = 82: 96,
                 "AtomicMass_MeV" = 98: 105,
                 "MassUncertainty_MeV" = 106: 113,
                 "IfSystematics" = 115,
                 "HalfLife_sec" = 125: 132)

inequal <- c("LT", "LE", "GT", "GE", "AP")

csv_data <- data %>%rowwise()%>%
  mutate("A" = (strsplit(strData, ""))[[1]][cardName$"A"]%>% paste(collapse = ""))%>%
  mutate("Isomeric" = (strsplit(strData, ""))[[1]][cardName$"Isomeric"]%>% paste(collapse = "")) %>%
  mutate("Z" = (strsplit(strData, ""))[[1]][cardName$"Z"]%>% paste(collapse = ""))%>%
  mutate("Element" = (strsplit(strData, ""))[[1]][cardName$"Element"]%>% paste(collapse = ""))%>%
  mutate("Jp" = (strsplit(strData, ""))[[1]][cardName$"Jp"]%>% paste(collapse = ""))%>%
  mutate("DecayMode" = (strsplit(strData, ""))[[1]][cardName$"DecayMode"]%>% paste(collapse = ""))%>%
  mutate("BranchRatio" = (strsplit(strData, ""))[[1]][cardName$"BranchRatio"]%>% paste(collapse = ""))%>%
  mutate("Ex_MeV" = (strsplit(strData, ""))[[1]][cardName$"Ex_MeV"]%>% paste(collapse = ""))%>%
  mutate("Q_MeV" = (strsplit(strData, ""))[[1]][cardName$"Q_MeV"]%>% paste(collapse = ""))%>%
  mutate_at(vars("Q_MeV"), ~ gsub("S", "",.))%>%
  mutate("HalfLife" = (strsplit(strData, ""))[[1]][cardName$"HalfLife"]%>% paste(collapse = ""))%>%
  mutate("Status" = ifelse(grepl("STABLE", HalfLife), "stable", 
                           ifelse(grepl("unbound", HalfLife), "unbound", "unstable")))%>%
  mutate("Hf_inequality" = ifelse(Status=="unstable", ifelse(
    length(which(sapply(inequal, grepl, HalfLife))) == 0, "EQUAL", 
    inequal[which(sapply(inequal, grepl, HalfLife))]
  ), NA))%>%
  mutate("Hf" = ifelse(Status=="unstable",  strsplit(HalfLife, " ")[[1]][1], NA))%>%
  mutate("Hf_err" = ifelse(Status=="unstable" & nchar(Hf) > 0,  sub(Hf, "", HalfLife, fixed = T)%>%gsub("[a-zA-Z ]","", .), NA))%>%
  mutate("Hf_err_u" = 10^(decimalnumcount(Hf))*as.numeric(ifelse(grepl("-",Hf_err),str_extract(Hf_err, "\\+[\\d]+"), Hf_err)))%>%
  mutate("Hf_err_d" = 10^(decimalnumcount(Hf))*abs(as.numeric(ifelse(grepl("-",Hf_err), str_extract(Hf_err, "-[\\d]+"), Hf_err))))%>%
  mutate("Hf_unit" = ifelse(Status=="unstable" &  nchar(HalfLife) > 0,  str_extract(HalfLife, "\ [a-zA-Z]+"),
                            NA))%>%
  select(-c(HalfLife, Hf_err))%>%
  mutate("Abundance" = (strsplit(strData, ""))[[1]][cardName$"Abundance"]%>% paste(collapse = ""))%>%
  mutate("Abundance_per" = strsplit(Abundance, "%")[[1]]%>%nth(1))%>%
  mutate("Abundance_err_per" = (10^(decimalnumcount(Abundance_per))*as.numeric(strsplit(Abundance, "%")[[1]]%>%nth(2))))%>%
  select(-c(Abundance))%>%
  mutate("Mass" = (strsplit(strData, ""))[[1]][cardName$"AtomicMass_MeV"]%>% paste(collapse = ""))%>%
  mutate("Mass_err" = (strsplit(strData, ""))[[1]][cardName$"MassUncertainty_MeV"]%>% paste(collapse = ""))%>%
  mutate("IfSystematics" = (strsplit(strData, ""))[[1]][cardName$"IfSystematics"]%>% paste(collapse = ""))%>%
  mutate("HalfLife_sec" = (strsplit(strData, ""))[[1]][cardName$"HalfLife_sec"]%>% paste(collapse = ""))%>%
  mutate_at(vars("A", "Z"), ~as.integer(.))%>%
  mutate_at(vars("Isomeric", "Element", "Jp", "DecayMode", "BranchRatio","IfSystematics","HalfLife_sec","Hf_unit"), ~gsub(" ","", .))%>%
  mutate_at(vars( "Ex_MeV", "Q_MeV", "Mass", "Mass_err", "Abundance_per", "Abundance_err_per"), 
            ~  as.numeric(.))%>%
  as.data.frame()%>%
  rowwise()%>%mutate_at(vars(Element), ~ str_to_title(.)) %>%
  mutate(Nuclide = paste(A, Element, sep = "")) %>% relocate(Nuclide)%>%
  select(-c(strData, Element, A))

write.csv(csv_data, "./NuData.csv", na = "", row.names = FALSE, quote=FALSE)
