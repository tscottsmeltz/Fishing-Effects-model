### Edited on 10/11/2017...

library(foreign)
library(reshape2)
library(rgdal)
library(maptools)



# Fishing Effects Model function
FishingEffectsModel = function(I.prime_a, rho.prime_a, H_prop_0){
  nTimeSteps = dim(I.prime_a)[1]
  nGrid = dim(I.prime_a)[2]
  nHab = dim(I.prime_a)[3]
  
  #Make array to hold H
  H_prop = array(dim = c(nTimeSteps, nGrid, nHab))  
  
  for(t in 1:nTimeSteps){
      
      if(t == 1){      # First time step use H_prop_0 for t-1
        prior_state = H_prop_0
      } else{                 
        prior_state = H_prop[t-1,,]
      }
      
      H_from_H = (1-I.prime_a[t,,])*prior_state  # undisturbed remaining undisturbed
      H_from_h = (1-prior_state) * (rho.prime_a[t,,]) # disturbed recovered to undisturbed
      H_prop[t,,] = H_from_H + H_from_h  # Total proportion disturbed
      
    }
  
  return(H_prop)
}  # end function


# Suceptibility 

suscept.f = function(SuscTables, SuscCat, HabList, SuscCodes){
	gear.q = matrix(NA, nrow = length(SuscCat), ncol = length(HabList))
	i = 1
	for(gear in SuscCat){
	  gear.m = as.matrix(SuscTables[[2]][which(SuscTables[[1]] == gear)][[1]])
							  

	  gear.m = gear.m[,HabList]
	  
	  nSuscCodes = nrow(SuscCodes)
	  
	  for(column in 1:ncol(gear.m)){
	     for(j in 1:nrow(SuscCodes)){
		      gear.m[gear.m[,column] %in% SuscCodes[j,1], column] = 
		         runif(sum(gear.m[,column] %in% SuscCodes[j,1]), min = SuscCodes[j,2], max = SuscCodes[j,3])
		
	     }
	  }
	  
	  gear.q[i,] = colMeans(gear.m, na.rm=T)
	  i = i+1
	  
	}

	gear.q.df = data.frame(SuscCat = SuscCat, gear.q)
	names(gear.q.df)[-1] = HabList



	q_m = as.matrix(gear.q.df[,HabList])
	q_m[is.na(q_m)] = 0
	return(q_m)
}


# Recovery function



recovery.f = function(RecoveryTable, HabList, RecoveryCodes){
  RecoveryTable = RecoveryTable[,HabList]
  tau_m = RecoveryTable[,HabList] # Make sure sediments are in correct order

	for(column in 1:ncol(tau_m)){
	  for(j in 1:nrow(RecoveryCodes)){
	    tau_m[RecoveryTable[,column] %in% RecoveryCodes[j,1], column] = 
		  runif(sum(tau_m[,column] %in% RecoveryCodes[j,1]), min = RecoveryCodes[j,2], max = RecoveryCodes[j,3])
	  }
	}


	tau_v = colMeans(tau_m, na.rm=T) # Average recovery over all habitat features

	rho_v = 1 / tau_v
	
	return(rho_v)

}


#' Fishing Effects model runner
#'
#' This function runs the Fishing Effects model
#' @param fullgrid SpatialPolygonsDataFrame object that contains all grid cells (or other polygon shapes) with GridID and Area attributes
#' @param habitats Dataframe object of habitat proportions in each grid cell.  
#' @param f.effort Dataframe object with the following columns: "GridID", "Gear", "t", "Sum" 
#' @param GearTable Dataframe object with "Gear", "SuscTab", "ContactMin", "ContactMax" columns.
#' @param RecoveryTable Dataframe object with 
#' @param SuscTables List object with dataframe in each list element corresponding to susceptibility table used 
#' @param RecoveryCodes Dataframe object with "Code", "Min", and "Max" columns 
#' @param SuscCodes Dataframe object with "Code", "Min", and "Max" columns
#' @param InitCond Dataframe object with "GridID" and "H" corresponding to percent intact habitat in each grid cell at start of model run.   
#' @param outPath
#' @param outFile
#' @param logFile
#' @param logText

#' @keywords fishing effects
#' @export
#' @examples
#' cat_function()





	
FEM <- function(
		fullgrid,
		habitats,
		f.effort,
		GearTable,
		RecoveryTable,
		SuscTables,
		RecoveryCodes,
		SuscCodes,
		InitCond,
		max.time = NULL,
		outPath,
		outFile,
		logFile = TRUE,
		logText = "" ){
		


# Import data
if(is.null(max.time)) {
	t.tot = max(f.effort$t)
	}else {
	t.tot = max.time
	}
	
HabList = names(habitats)[-1]


# Merge fishing effort and gear table
fe = merge(f.effort, GearTable,
	   by.x = "GearID", by.y = "GearID", all.x = T)
fe$SuscCat = factor(fe$SuscCat)

fe$adjArea = fe$nomArea * runif(nrow(fe), min = fe$caMin, max = fe$caMax) ## Random uniform contact adj from min and max

		
## aggregate adj effort based on susceptibiltiy group
fe.agg = aggregate(adjArea ~ GridID + t + SuscCat, data = fe, sum)

grid_order = sort(unique(fe$GridID)) # establish consistent order of GridID


# Create habitat profile matrix for model.  Make sure grid order is same as I_a 
# and keep only grid cells with fishing effort
habProps = as.matrix(habitats[ match(grid_order, habitats$GridID),])



# Set parameters
nGrid = length(unique(fe$GridID))
SuscCat = levels(fe$SuscCat)
nSuscCat = length(SuscCat)
nHab = length(HabList)



eg = expand.grid(GridID=unique(fe$GridID), SuscCat = unique(fe$SuscCat)) # create all combos of SuscCat and grid cells

m = merge(x = fe.agg, y = fullgrid[,c("GridID", "Area")], by = "GridID", all.x = T)

m$prop = m$adjAre/m$Area


# Populate Fishing effort array
F_a = array(NA, dim = c(t.tot, nGrid, nSuscCat)) #create empty array


for(i in 1:t.tot){
  mym = subset(m, t == i)

    mym = merge(x = eg, y = mym, 
                by = c("GridID","SuscCat"), 
                all.x=T)
    
    mym[is.na(mym$prop),]$prop = 0
    
    mym.x = dcast(GridID ~ SuscCat, data=mym, value.var = "prop", fun.aggregate = function(x) sum(x))
    
    mym.x = mym.x[order(mym.x$GridID),]
    mym.x = mym.x[,SuscCat]
    
    F_a[i, ,] = as.matrix(mym.x)
    
    
  }








#Fishing impacts (I') 
I.prime_a = array(NA, dim = c(t.tot, nGrid, nHab))

for(t in 1:t.tot){
	q_m = suscept.f(SuscTables = SuscTables, SuscCat = SuscCat, HabList = HabList, SuscCodes = SuscCodes)  # Get new susceptibility table for each month
    I_m = F_a[t,,] %*% q_m
    I.prime_a[t,,] = 1-exp(-I_m)
  }




## Recovery
rho.prime_a = array(NA, dim = c(t.tot, nGrid, nHab))


for(t in 1:t.tot){
	rho_v = recovery.f(RecoveryTable = RecoveryTable, HabList = HabList, RecoveryCodes = RecoveryCodes)  # Get new recovery values for each month
    for(i in 1:nGrid){
		rho.prime_a[t,i,] = 1-exp(-rho_v)
    }
}

##### Run model
H_prop_0 = as.matrix(InitCond[match(grid_order, InitCond$GridID), -1])
H_tot = FishingEffectsModel(I.prime_a, rho.prime_a, H_prop_0 = H_prop_0)


undistProps = matrix(NA, ncol = t.tot, nrow = length(grid_order))

habitats_m = habitats[match(grid_order, habitats$GridID),-1]
for(t in 1:t.tot){
    undistProps[,t] = rowSums(H_tot[t,,]*habitats_m)
  }

undistProps = data.frame(GridID = grid_order, undistProps)


disturbProps = data.frame(GridID = undistProps[,1], 
                          apply(undistProps[,-1],2,
                                function(x) 1 - x))								

all.grid = data.frame(GridID = fullgrid$GridID)
disturbProps = merge(x = all.grid, y = disturbProps, by = "GridID", all.x = T)								
disturbProps[is.na(disturbProps)] = 0		
names(disturbProps)[-1] = paste("t", 1:(ncol(disturbProps)-1), sep = "")
return(disturbProps)							
								
}















