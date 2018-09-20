#see https://www.ini.uzh.ch/~acardona/fiji-tutorial/#s1 for help!!
#labels can be set using roiManager
#use with genFIJI_ROIs.m
from ij import IJ
from ij.gui import Roi, OvalRoi
from ij.plugin.frame import RoiManager
from ij.io import OpenDialog
import os

roiManager = RoiManager.getRoiManager()
fileSelection = OpenDialog("Select the positions file")
fileString = IJ.openAsString(fileSelection.getPath())
dataRows = fileString.split("\n")

#print fileString
#print fileSelection.getPath(), 
#print fileSelection.getFileName()
#Use names as labels in roiManager menu (more>>>)
#Save manually as well
imp = IJ.getImage()

rowIter = 0 

for row in dataRows:
	if len(row)>0:
		spotDataList = row.split("\t")
		spotData = [float(x) for x in spotDataList]
		topLeftX = spotData[1] - spotData[4]
		topLeftY = spotData[2] - spotData[4]
		diameter = 2*spotData[4]
		spotRoi = OvalRoi(topLeftX,topLeftY,diameter,diameter)
		roiManager.add(imp,spotRoi,rowIter)
		roiManager.rename(rowIter,"%4.0f"%spotData[0])
		rowIter = rowIter+1
	





#IJ.makeOval(900,900,50,80)
#a = OvalRoi(500,500,100,200)
#b = OvalRoi(200,200,300,100)

#roiManager.add(imp,a,0)
#roiManager.add(imp,b,0)
#roiManager.rename(0,"Spot 1")
#roiManager.rename(1,"Spot 2")
#roiManager.runCommand("Save","savedrois.zip")
#RoiManager.addRoi(a,0)

#print os.getcwd()

#imp.setRoi(test2)