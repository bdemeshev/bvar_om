# bvar. OM
if (Sys.info()[7]=="boris") {
  Sys.setenv(X13_PATH ="/usr/local/bin/x13tramo")
  setwd("~/Documents/bvar/mal_art")
} else {
  Sys.setenv(X13_PATH="D:/Oxana/Science/Projects/ProjectXIV/x13as/")
  setwd("D:/Oxana/Science/Projects/ProjectXIV/mal_art/mal_art/")
}

Sys.setenv(LANGUAGE="en")
library("knitr")
opts_chunk$set(echo=FALSE, message=FALSE, warning=FALSE)

library("pander")
library("stringr")
library("reshape2")

library("BMR")
library("bvarr")
library("MSBVAR")
library("vars")

library("seasonal")
checkX13()

library("ggplot2")
library("mvtsplot")

library("data.table")
library("dplyr")
library("lubridate")

library("ftsa")

