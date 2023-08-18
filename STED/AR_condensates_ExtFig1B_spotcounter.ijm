
// Based on RB_220319_spotcounter.ijm by Rene Buschow (buschow@molgen.mpg.de), modified for STED FLIM by Shaon Basu
//// This is the Start////

// attention this part is meant to run directly after another one

run("Set Measurements...", "area mean standard perimeter shape redirect=None decimal=3");

dir=getDirectory("Choose a Directory");

cellresults =dir+"result_cells";
File.makeDirectory(cellresults);
if (!File.exists(cellresults))
      exit("Sorry this is already existent!");

cellrois =dir+"cell_rois";
File.makeDirectory(cellrois);
if (!File.exists(cellrois))
      exit("Sorry this also already exists!");

list=getFileList(dir);

for(k=0;k<list.length;k++){ 


open(""+dir+list[k]+"");

run("Set Scale...", "distance=1 known=0.045 pixel=1 unit=�m");
image=getTitle();
run("Duplicate...", "title=cell_mask");
run("Smooth");
run("Smooth");
run("Smooth");
run("Smooth");
run("Subtract Background...", "rolling=4 sliding disable");
setThreshold(40, 250);
run("Convert to Mask");
run("Watershed");
run("Watershed");

run("Analyze Particles...", "size=0.005-1 circularity=0.60-1.00 show=Nothing exclude clear add");

rois=roiManager("count");

if(rois>0)

{
roiManager("Save",dir+"//cell_rois//"+image+"_ROIs.zip");
selectWindow(image);
roiManager("Deselect");
roiManager("Show All");
roiManager("Measure");
roiManager("Delete");
selectWindow("Results");
saveAs("Text",dir+"//result_cells//"+image+".txt");
run("Clear Results");
selectWindow(image);
run("Close");
selectWindow("cell_mask");
run("Close");

			}

		else	{

run("Clear Results");
selectWindow(image);
run("Close");
selectWindow("cell_mask");
run("Close");


}



}







