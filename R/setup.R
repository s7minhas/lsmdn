rm(list=ls())

if(Sys.info()['user']=='janus829' | Sys.info()['user']=='s7m'){
	inPath='~/Dropbox/Research/NetworkEvolution/Data/inData/';
	outPath='~/Dropbox/Research/NetworkEvolution/Data/outData/';
	graphicsPath='~/Dropbox/Research/NetworkEvolution/Graphics/';
	rFuncs='~/Research/NetworkEvolution/Code/helpers/';
}

if(Sys.info()['user']=='maxgallop' ){
	inPath='~/Dropbox/NetworkEvolution/Data/inData/';
	outPath='~/Dropbox/NetworkEvolution/Data/outData/';
	graphicsPath='~/Dropbox/NetworkEvolution/Graphics/';
	rFuncs='~/Documents/NetworkEvolution/Code/helpers/';
}

if(Sys.info()['user']=='cassydorff'){
	inPath='~/Dropbox/Research/NetworkEvolution/Data/inData/';
	outPath='~/Dropbox/Research/NetworkEvolution/Data/outData/';
	graphicsPath='~/Dropbox/Research/NetworkEvolution/Graphics/';
	rFuncs='~/ProjectsGit/NetworkEvolution/Code/helpers/';
}


# General functions/libraries
loadPkg=function(toLoad){
	for(lib in toLoad){
	  if(!(lib %in% installed.packages()[,1])){ 
	    install.packages(lib, repos='http://cran.rstudio.com/') }
	  suppressMessages( library(lib, character.only=TRUE) )
	}
}

toLoad=c(
	'RMySQL', 
	'ggplot2', 'network', 'grid',
	'reshape2', 'magrittr', 'doBy',
	'countrycode',
	'xtable', 'tikzDevice',
	'foreach', 'doParallel'
	)
loadPkg(toLoad)

# Set a theme for gg
theme_set(theme_bw())

# Other functions
mysqlSetup = function(user=NULL, pw=NULL, db=NULL, host=NULL) {
	tryCatch(conn <<- dbConnect(MySQL(), user=user, password=pw, 
		dbname=db, host=host), 
	error=function(e) warning("MySQL connection does not work") )
}

# Global params
seed=6886
set.seed(seed)

# Helpful functions
char = function(x){ as.character(x) }
num = function(x){ as.numeric(char(x)) }
pasteVec = function(x,y){ as.vector(outer(x,y,paste0)) }

# panel dataset for matching country observations
load(paste0(inPath, 'panel.rda'))