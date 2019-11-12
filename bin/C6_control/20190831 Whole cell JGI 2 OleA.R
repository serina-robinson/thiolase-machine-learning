# Install packages
pacman::p_load("tidyverse")

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/synbio-data-analysis/")

# Read in the data
fils <- list.files("data/C6_control/", pattern = "20190909", full.names = T)
fils
Datab <- read_csv(fils, col_names = T)


Molarity=(Datab[,c(3:98)] + 0.0217)/(0.58*15546)
Table=(Molarity*0.0002)*1e+09
Time=c(0:60)

# calculate each condition
Blank=rowMeans(Table[,c("A1","A2")])

A6.2=rowMeans(Table[,c("B1", "B2")]-Blank)
A6.2.sd=apply((Table[,c("B1","B2")]-Blank), 1,sd)

A7.2=rowMeans(Table[,c("C1", "C2")]-Blank)
A7.2.sd=apply((Table[,c("C1", "C2")]-Blank), 1,sd)

A8.2=rowMeans(Table[,c("D1", "D2")]-Blank)
A8.2.sd=apply((Table[,c("D1", "D2")]-Blank), 1,sd)

A9.2=rowMeans(Table[,c("E1", "E2")]-Blank) 
A9.2.sd=apply((Table[,c("E1", "E2")]-Blank), 1,sd)

A10.2=rowMeans(Table[,c("F1", "F2")]-Blank)
A10.2.sd=apply((Table[,c("F1", "F2")]-Blank), 1,sd)

B5.2=rowMeans(Table[,c("G1", "G2")]-Blank)
B5.2.sd=apply((Table[,c("G1", "G2")]-Blank), 1,sd)

B6.2=rowMeans(Table[,c("H1", "H2")]-Blank)
B6.2.sd=apply((Table[,c("H1", "H2")]-Blank), 1,sd)

B7.2=rowMeans(Table[,c("A3", "A4")]-Blank)
B7.2.sd=apply((Table[,c("A3", "A4")]-Blank), 1,sd)

B8.2=rowMeans(Table[,c("B4", "B3")]-Blank)
B8.2.sd=apply((Table[,c("B4", "B3")]-Blank), 1,sd)

B9.2=rowMeans(Table[,c("C3","C4")]-Blank)
B9.2.sd=apply((Table[,c("C3")]-Blank), 1,sd)

B10.2=rowMeans(Table[,c("D3","D4")]-Blank)
B10.2.sd=apply((Table[,c("D3")]-Blank), 1,sd)

C5.2=rowMeans(Table[,c("E3", "E4")]-Blank)
C5.2.sd=apply((Table[,c("E3", "E4")]-Blank), 1,sd)

C6.2=rowMeans(Table[,c("F3", "F4")]-Blank) 
C6.2.sd=apply((Table[,c("F3", "F4")]-Blank), 1,sd)

C7.2=rowMeans(Table[,c("G3", "G4")]-Blank)
C7.2.sd=apply((Table[,c("G3", "G4")]-Blank), 1,sd)

C8.2=rowMeans(Table[,c("H3", "H4")]-Blank)
C8.2.sd=apply((Table[,c("H3", "H4")]-Blank), 1,sd)

C9.2=rowMeans(Table[,c("A5", "A6")]-Blank)
C9.2.sd=apply((Table[,c("A5", "A6")]-Blank), 1,sd)

C10.2=rowMeans(Table[,c("B5", "B6")]-Blank)
C10.2.sd=apply((Table[,c("B5", "B6")]-Blank), 1,sd)

D5.2=rowMeans(Table[,c("C5", "C6")]-Blank)
D5.2.sd=apply((Table[,c("C5", "C6")]-Blank), 1,sd)

D6.2=rowMeans(Table[,c("D6", "D5")]-Blank)
D6.2.sd=apply((Table[,c("D5", "D6")]-Blank), 1,sd)

D7.2=rowMeans(Table[,c("E5", "E6")]-Blank) 
D7.2.sd=apply((Table[,c("E5", "E6")]-Blank), 1,sd)

D8.2=rowMeans(Table[,c("F5", "F6")]-Blank)
D8.2.sd=apply((Table[,c("F5", "F6")]-Blank), 1,sd)

D10.2=rowMeans(Table[,c("G5", "G6")]-Blank)
D10.2.sd=apply((Table[,c("G5", "G6")]-Blank), 1,sd)

E5.2=rowMeans(Table[,c("H5", "H6")]-Blank)
E5.2.sd=apply((Table[,c("H5", "H6")]-Blank), 1,sd)

E6.2=rowMeans(Table[,c("A7", "A8")]-Blank)
E6.2.sd=apply((Table[,c("A7", "A8")]-Blank), 1,sd)

E7.2=rowMeans(Table[,c("B7", "B8")]-Blank)
E7.2.sd=apply((Table[,c("B7", "B8")]-Blank), 1,sd)

E8.2=rowMeans(Table[,c("C7", "C8")]-Blank)
E8.2.sd=apply((Table[,c("C7", "C8")]-Blank), 1,sd)

E9.2=rowMeans(Table[,c("D7", "D8")]-Blank) 
E9.2.sd=apply((Table[,c("D7", "D8")]-Blank), 1,sd)

E10.2=rowMeans(Table[,c("E7", "E8")]-Blank)
E10.2.sd=apply((Table[,c("E7", "E8")]-Blank), 1,sd)

F5.2=rowMeans(Table[,c("F7", "F8")]-Blank)
F5.2.sd=apply((Table[,c("F7", "F8")]-Blank), 1,sd)

F6.2=rowMeans(Table[,c("G7", "G8")]-Blank)
F6.2.sd=apply((Table[,c("G7", "G8")]-Blank), 1,sd)

F7.2=rowMeans(Table[,c("H7", "H8")]-Blank)
F7.2.sd=apply((Table[,c("H7", "H8")]-Blank), 1,sd)

F8.2=rowMeans(Table[,c("A9", "A10")]-Blank)
F8.2.sd=apply((Table[,c("A9", "A10")]-Blank), 1,sd)

F9.2=rowMeans(Table[,c("B9", "B10")]-Blank)
F9.2.sd=apply((Table[,c("B9", "B10")]-Blank), 1,sd)

F10.2=rowMeans(Table[,c("C9", "C10")]-Blank)
F10.2.sd=apply((Table[,c("C9", "C10")]-Blank), 1,sd)

G5.2=rowMeans(Table[,c("D9", "D10")]-Blank) 
G5.2.sd=apply((Table[,c("D9", "D10")]-Blank), 1,sd)

G6.2=rowMeans(Table[,c("E9", "E10")]-Blank)
G6.2.sd=apply((Table[,c("E9", "E10")]-Blank), 1,sd)

G7.2=rowMeans(Table[,c("F9", "F10")]-Blank)
G7.2.sd=apply((Table[,c("F9", "F10")]-Blank), 1,sd)

G8.2=rowMeans(Table[,c("G9", "G10")]-Blank)
G8.2.sd=apply((Table[,c("G9", "G10")]-Blank), 1,sd)

G9.2=rowMeans(Table[,c("H9", "H10")]-Blank)
G9.2.sd=apply((Table[,c("H9", "H10")]-Blank), 1,sd)

G10.2=rowMeans(Table[,c("A11", "A12")]-Blank)
G10.2.sd=apply((Table[,c("A11", "A12")]-Blank), 1,sd)

H5.2=rowMeans(Table[,c("B11", "B12")]-Blank)
H5.2.sd=apply((Table[,c("B11", "B12")]-Blank), 1,sd)

H6.2=rowMeans(Table[,c("C11", "C12")]-Blank)
H6.2.sd=apply((Table[,c("C11", "C12")]-Blank), 1,sd)

H7.2=rowMeans(Table[,c("D11", "D12")]-Blank) 
H7.2.sd=apply((Table[,c("D11", "D12")]-Blank), 1,sd)

H8.2=rowMeans(Table[,c("E11", "E12")]-Blank)
H8.2.sd=apply((Table[,c("E12", "E12")]-Blank), 1,sd)

H9.2=rowMeans(Table[,c("F11", "F12")]-Blank)
H9.2.sd=apply((Table[,c("F12", "F12")]-Blank), 1,sd)

H10.2=rowMeans(Table[,c("G11", "G12")]-Blank)
H10.2.sd=apply((Table[,c("G11", "G12")]-Blank), 1,sd)

F2.1=rowMeans(Table[,c("H11","H12")]-Blank)
F2.1.sd=apply((Table[,c("H11","H12")]-Blank),1,sd)


#graphs-chain length graph
library(reshape2)
library(ggplot2)
# chain_table=data.frame(Time,A6.2,A7.2,A8.2,A9.2,A10.2,B5.2,B6.2,B7.2,B8.2,B9.2,B10.2)
                       
chain_table = data.frame(Time, A6.2,A7.2,A8.2,A9.2,A10.2,B5.2,B6.2,B7.2,B8.2,B9.2,B10.2,C5.2,C6.2,C7.2,C8.2,C9.2,C10.2,D5.2,D6.2,D7.2,D8.2,D10.2,E5.2,E6.2,E7.2,E8.2,E9.2,E10.2,F5.2,F6.2,F7.2,F8.2,F9.2,F10.2,G5.2,G6.2,G7.2,G8.2,G9.2,G10.2,H5.2,H6.2,H7.2,H8.2,H9.2,H10.2)

#D2,G2,F2,B3,D3,E3,F3,G3,H3,A4,B4,C4,D4,E4,F4,G4,H4,A5,B5,C5,D5,E5,F5,G5,H5,A6,B6
chain.melt=melt(chain_table, id.vars="Time")
head(chain.melt)

write_csv(chain.melt, "data/C6_control/2019-09-09_whole_cell_JGI_2_melted_rep_2.csv")

# chain_SD=data.frame(Time,D2.sd)
chain_SD.melt=melt(chain_SD, id.vars="Time")

Chain.ymin=chain.melt$value - chain_SD.melt$value
Chain.ymax=chain.melt$value + chain_SD.melt$value

cbbPalette <- c("#DCDCDC", "#FF0000", "#00FF00", "#FFFF00", "#FF00FF", "#00FFFF", "#00FA9A", "#FF8080", "#0066CC","#008080", "#CC99FF", "#FFCC00", "#666699", "#FFB6C1", "#DB7093", "#333333", "#00FFFF", "#CCCCFF", "#C0C0C0", "#808000", "#FFE4E1", "#008000", "#800080", "#99CC00", "#FF6600","#7FFFD4", "#0000FF", "#A52A2A", "#7FFF00")

y.axis <- expression(paste("nmol ", italic("p"), "-NP"))


ggplot(chain.melt, aes(x=Time, y=value, color=variable)) +
  labs( y=y.axis, x="Time (minutes)")+
  geom_point(size=3)+
  theme(legend.title=element_blank(), axis.line=element_line(color="black"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        panel.background=element_blank(),
        text = element_text(size = 20),
        legend.key= element_rect(fill=NA, color=NA),
        legend.position=c(0.20,0.8))+
  scale_color_manual(labels=c("A6",
                              "A7",
                              "A8",
                              "A9",
                              "A10",
                              "B5",
                              "B6",
                              "B7",
                              "B8",
                              "B9",
                              "B10"),
                     values=cbbPalette)+
  coord_cartesian(xlim=c(0,60), ylim=(c(-2, 10)))





print(pls )

"D2",
"G2",
"F2",
"B3",
"D3",
"E3",
"F3",
"G3",
"H3",
"A4",
"B4",
"C4",
"D4",
"E4",
"F4",
"G4",
"H4",
"A5",
"B5",
"C5",
"D5",
"E5",
"F5",
"G5",
"H5",
"A6",
"B6",
"X.c OleA"
