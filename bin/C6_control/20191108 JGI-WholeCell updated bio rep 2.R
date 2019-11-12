# Install packages
pacman::p_load("tidyverse", "readxl")

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/synbio-data-analysis/")

# Read in the data
fils <- list.files("data/C6_control/", pattern = "20191108", full.names = T)
fils
Datab <- read_excel(fils, col_names = T)
Datab
Molarity=(Datab[,c(3:98)] + 0.0217)/(0.58*15546)
Table=(Molarity*0.0002)*1e+09
Time=c(0:60)

# calculate each condition
Blank=rowMeans(Table[,c("E7", "E8")])

D2.1=rowMeans(Table[,c("A1", "A2")]-Blank)
D2.sd=apply((Table[,c("A1", "A2")]-Blank), 1,sd)
D2.sd

F2.1=rowMeans(Table[,c("B1","B2")]-Blank)
F2.sd=apply((Table[,c("B1", "B2")]-Blank), 1,sd)

G2.1=rowMeans(Table[,c("C1", "C2")]-Blank)
G2.sd=apply((Table[,c("C1", "C2")]-Blank), 1,sd)

B3.1=rowMeans(Table[,c("D1", "D2")]-Blank) 
B3.sd=apply((Table[,c("D1", "D2")]-Blank), 1,sd)

D3.1=rowMeans(Table[,c("E1", "E2")]-Blank)
D3.sd=apply((Table[,c("E1", "E2")]-Blank), 1,sd)

E3.1=rowMeans(Table[,c("F1", "F2")]-Blank)
E3.sd=apply((Table[,c("F1", "F2")]-Blank), 1,sd)

F3.1=rowMeans(Table[,c("G1", "G2")]-Blank)
F3.sd=apply((Table[,c("G1", "G2")]-Blank), 1,sd)

G3.1=rowMeans(Table[,c("H1", "H2")]-Blank)
G3.sd=apply((Table[,c("H1", "H2")]-Blank), 1,sd)

H3.1=rowMeans(Table[,c("A4", "A3")]-Blank)
H3.sd=apply((Table[,c("A4", "A3")]-Blank), 1,sd)

A4.1=rowMeans(Table[,c("B4", "B3")]-Blank)
A4.sd=apply((Table[,c("B4", "B3")]-Blank), 1,sd)

B4.1=rowMeans(Table[,c("C4", "C3")]-Blank)
B4.sd=apply((Table[,c("C4", "C3")]-Blank), 1,sd)

C4.1=rowMeans(Table[,c("D4", "D3")]-Blank) 
C4.sd=apply((Table[,c("D4", "D3")]-Blank), 1,sd)

D4.1=rowMeans(Table[,c("E4", "E3")]-Blank)
D4.sd=apply((Table[,c("E4", "E3")]-Blank), 1,sd)

E4.1=rowMeans(Table[,c("F4", "F3")]-Blank)
E4.sd=apply((Table[,c("F4", "F3")]-Blank), 1,sd)

F4.1=rowMeans(Table[,c("G4", "G3")]-Blank)
F4.sd=apply((Table[,c("G4", "G3")]-Blank), 1,sd)

G4.1=rowMeans(Table[,c("H4", "H3")]-Blank)
G4.sd=apply((Table[,c("H4", "H3")]-Blank), 1,sd)

H4.1=rowMeans(Table[,c("A5", "A6")]-Blank)
H4.sd=apply((Table[,c("A5", "A6")]-Blank), 1,sd)

A5.1=rowMeans(Table[,c("B5", "B6")]-Blank)
A5.sd=apply((Table[,c("B5", "B6")]-Blank), 1,sd)

B5.1=rowMeans(Table[,c("C5", "C6")]-Blank)
B5.sd=apply((Table[,c("C5", "C6")]-Blank), 1,sd)

C5.1=rowMeans(Table[,c("D5", "D6")]-Blank) 
C5.sd=apply((Table[,c("D5", "D6")]-Blank), 1,sd)

D5.1=rowMeans(Table[,c("E5","E6")]-Blank)
D5.sd=apply((Table[,c("E5", "E6")]-Blank), 1,sd)

E5.1=rowMeans(Table[,c("F5", "F6")]-Blank)
E5.sd=apply((Table[,c("F5", "F6")]-Blank), 1,sd)

F5.1=rowMeans(Table[,c("G5", "G6")]-Blank)
F5.sd=apply((Table[,c("G5", "G6")]-Blank), 1,sd)

G5.1=rowMeans(Table[,c("H5", "H6")]-Blank)
G5.sd=apply((Table[,c("H5", "H6")]-Blank), 1,sd)

H5.1=rowMeans(Table[,c("A7", "A8")]-Blank)
H5.sd=apply((Table[,c("A7", "A8")]-Blank), 1,sd)

A6.1=rowMeans(Table[,c("B7", "B8")]-Blank)
A6.sd=apply((Table[,c("B7", "B8")]-Blank), 1,sd)

B6.1=rowMeans(Table[,c("C7", "C8")]-Blank)
B6.sd=apply((Table[,c("C7", "C8")]-Blank), 1,sd)

OleA=rowMeans(Table[,c("D7", "D8")]-Blank) 
OleA.sd=apply((Table[,c("D7", "D8")]-Blank), 1,sd)

A2=rowMeans(Table[,c("E7", "E8")]-Blank)
A2.sd=apply((Table[,c("E7", "E8")]-Blank), 1,sd)

B2=rowMeans(Table[,c("F10", "F11", "F12")]-Blank)
B2.sd=apply((Table[,c("F10", "F11", "F12")]-Blank), 1,sd)

E2=rowMeans(Table[,c("G10", "G11", "G12")]-Blank)
E2.sd=apply((Table[,c("G10", "G11", "G12")]-Blank), 1,sd)

#graphs-chain length graph
library(reshape2)
library(ggplot2)
library(ggpubr)
chain_table=data.frame(Time,F2.1,C5.1,A4.1,B4.1,OleA,B5.1,D2.1,G2.1,B3.1,D3.1,E3.1,F3.1,G3.1,H3.1,C4.1,D4.1,E4.1,F4.1,G4.1,H4.1,A5.1,D5.1,E5.1,F5.1,G5.1,H5.1,A6.1,B6.1)
chain.melt=melt(chain_table, id.vars="Time")
write_csv(chain.melt, "data/C6_control/2019-11-08_whole_cell_JGI_1_melted_bio_rep_2.csv")      

chain_SD=data.frame(Time,F2.sd,C5.sd,A4.sd,B4.sd,OleA.sd,B5.sd,D2.sd,G2.sd,B3.sd,D3.sd,E3.sd,F3.sd,G3.sd,H3.sd,C4.sd,D4.sd,E4.sd,F4.sd,G4.sd,H4.sd,A5.sd,D5.sd,E5.sd,F5.sd,G5.sd,H5.sd,A6.sd,B6.sd)



chain_SD.melt=melt(chain_SD, id.vars="Time")
dim(chain_SD.melt)
dim(chain.melt)

Chain.ymin=chain.melt$value - chain_SD.melt$value

Chain.ymax=chain.melt$value + chain_SD.melt$value

cbbPalette <- c("#DCDCDC", "#FF0000", "#00FF00", "#FFFF00", "#FF00FF", "#0072B2", "#00FA9A", "#FF8080", "#0066CC","#008080", "#CC99FF", "#FFCC00", "#666699", "#FFB6C1", "#DB7093", "#333333", "#00FFFF", "#CCCCFF", "#C0C0C0", "#808000", "#FFE4E1", "#008000", "#800080", "#99CC00", "#FF6600","#7FFFD4", "#0000FF", "#A52A2A", "#7FFF00")


cbbPaletteGrey <- c("#99CC00", "#B526AB", "#0066CC", "#C17F0D", "#000000", "#D3D3D3", "#D3D3D3", "#D3D3D3", "#D3D3D3","#D3D3D3", "#D3D3D3", "#D3D3D3", "#D3D3D3", "#D3D3D3", "#D3D3D3", "#D3D3D3", "#D3D3D3", "#D3D3D3", "#D3D3D3", "#D3D3D3", "#D3D3D3", "#D3D3D3", "#D3D3D3", "#D3D3D3", "#D3D3D3","#D3D3D3", "#D3D3D3", "#D3D3D3", "#D3D3D3")

something <-c("circle","circle","circle","circle","square","circle","circle","circle","circle","circle","circle","circle","circle","circle","circle","circle","circle","circle","circle","circle","circle","circle","circle","circle","circle","circle","circle","circle")

y.axis <- expression(paste("nmol ", italic("p"), "-NP"))


ggplot(chain.melt, aes(x=Time, y=value, color=variable)) +
  labs( y=y.axis, x="Time (minutes)")+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=Chain.ymin,ymax=Chain.ymax),size=2)+
  theme_classic()+
  scale_color_manual(labels=c("F2",
                              "C5",
                              "A4",
                              "B4",
                              "OleA",
                              "D2",
                              "G2",
                              "B3",
                              "D3",
                              "E3",
                              "F3",
                              "G3",
                              "H3",
                              "C4",
                              "D4",
                              "E4",
                              "F4",
                              "G4",
                              "H4",
                              "A5",
                              "B5",
                              "D5",
                              "E5",
                              "F5",
                              "G5",
                              "H5",
                              "A6",
                              "B6"),
                     values=cbbPaletteGrey)+
  scale_shape_manual(values=something)+
  coord_cartesian(xlim=c(0,40.5), ylim=(c(-1, 30)))

dev.off()

