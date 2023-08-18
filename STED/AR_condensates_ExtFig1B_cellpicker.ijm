

//// This is the START////

// Based on RB_211124_cell_picking by Rene Buschow (buschow@molgen.mpg.de), modified for STED FLIM by Shaon Basu

// regular bugs be aware of that the desired folder and image names are not allowed to contain free position e.g. "a b" should be "a_b"
// set intensity threshold for overall cell picking

lo_limit=1
up_limit=250



//setBatchMode(true)

//specify directory < where single channel tiff images are in

dir=getDirectory("Choose a Directory");
print(dir)

// get image or subset image
run("Image Sequence...", "open="+dir+" number=1 starting=1 increment=1 scale=100 file=[] or=[] sort");
run("8-bit");
run("Set Scale...", "distance=1 known=1 pixel=1 unit=�m");
//makeRectangle(815, 3114, 486, 501);
//run("Crop");
rename("raw_image");


// get nuclei segmented
run("Duplicate...", "title=prim_object duplicate range=1-1");
run("Smooth");
run("Smooth");
run("Smooth");
run("Smooth");
run("Subtract Background...", "rolling=10 sliding disable");
setThreshold(2, 150);
run("Convert to Mask");
run("Close-");
run("Close-");
run("Close-");
run("Close-");
run("Close-");
run("Close-");
run("Close-");
run("Close-");
run("Close-");
run("Close-");
run("Close-");
run("Close-");
run("Close-");
run("Close-");
run("Close-");




run("Dilate");
run("Dilate");
run("Dilate");
run("Dilate");
run("Dilate");



run("Fill Holes");

run("Watershed")

// get ROIs specify cell size and circularity
run("Analyze Particles...", "size=35000-1000000 circularity=0-1.00 show=Nothing display exclude clear add");
run("Set Measurements...", "  bounding redirect=None decimal=3");

// save ROIs
maskdir =dir+"ROIs";
File.makeDirectory(maskdir);
roiManager("Deselect");
roiManager("Save", ""+dir+"/ROIs/RoiSet.zip");
if (!File.exists(maskdir))
      exit("Fool this folder was already analyzed");


// get single cell images > bounding rectangle > specify band width
celldir =dir+"single_cells";
File.makeDirectory(celldir);
if (!File.exists(celldir))
      exit("Strange that you made it until here, check your folder.");

// Adjust number of desired cells 

selectWindow("prim_object");
run("Close");
selectWindow("raw_image");
run("Close");


// get image or subset image

run("Image Sequence...", "open="+dir+"Condensate"+" number=1 starting=1 increment=1 scale=100 file=[] or=[] sort");

run("8-bit");
run("Set Scale...", "distance=1 known=1 pixel=1 unit=�m");
//makeRectangle(815, 3114, 486, 501);
//run("Crop");
rename("condensate_image");
selectWindow("condensate_image");

rois=roiManager("count");
		for (j=0;j<rois;j++) 
{
			
		roiManager("Select", j);

		run("Set Measurements...", "mean bounding redirect=None decimal=3");
		run("Measure");
		par=getResult("Mean");
		run("Clear Results");

			if((lo_limit<=par)&&(par<=up_limit))

			{
			//run("Set Measurements...", "  bounding redirect=None decimal=3");
			//run("Measure");
			//boundx=getResult("BX");
			//boundy=getResult("BY");
			//boundwidth=getResult("Width");
			//boundheight=getResult("Height");
			selectWindow("condensate_image");
			//makeRectangle(boundx, boundy, boundwidth, boundheight);
			run("Duplicate...", "title=cell duplicate range=1-3");
			run("Clear Outside");
			saveAs("Tiff", ""+dir+"single_cells/cell_"+j+".tif");
			run("Close");
			run("Clear Results");

			}

		else	{}

}


selectWindow("condensate_image");

//oiManager("Delete");
//roiManager("Delete");

//dir=getDirectory("Choose a Directory");


run("Image Sequence...", "open="+dir+"Condensate"+" number=1 starting=1 increment=1 scale=100 file=[] or=[] sort");

//roiManager("Open", ""+dir+"/ROIs/RoiSet.zip");

run("Set Measurements...", "area mean min bounding redirect=None decimal=3");
roiManager("Measure");
saveAs("Results", ""+dir+"/ROIs/Results.csv");
roiManager("Delete");
roiManager("Delete");


run("Close");

//// This is the END////


