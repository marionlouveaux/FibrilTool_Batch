//////////////////////////////////////////////////////////////////////////////////
///////////////////////////FibrilTool Batch macro/////////////////////////////////
///////////////////////Automation of FibrilTool macro/////////////////////////////
/////////////////////////////Marion Louveaux//////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

// This macro is an automated version of FibrilTool macro written by Arezki Boudaoud (see Boudaoud et al., Nature Protocols 2014)
// Put .png (or .tif) and RoiSet.zip in a unique directory (as many as you want) = "dir"
// .png and RoiSet.zip must share the same name: ImgName.png and ImgName_RoiSet.zip
// One RoiSet.zip = one .png
// FibrilTool is executed on each .roi of the RoiSet
// Results, raw data (new RoiSet) and Flattened images are saved in dir




// the log output gives the average properties of the region
// 0) 	image title
//		cell number
// 1) 	x-coordinate of region centroid (scaled)
// 		y-coordinate of region centroid (scaled)
// 	 	area (scaled)
// 2)  nematic tensor
//		average orientation (angle in -90:90 in degrees)
// 		quality of the orientation (score between 0 and 1)
// The results are drawn on an overlay
// 3)  coordinates of polygon vertices for record



//threshold for flat regions
var thresh = 2;

//default for font size
var fsize = 15;

// number for cells
var num;

//default for width of lines
var lwidth = 2;

var pi = 3.14159265;


macro "Fibril Tool - C059L11eeL14beL41ebL178eL71e8L1a5eLa1e5" {

dir = getDirectory("Choose a folder.")

///Ask the user about fibril tool parameters
  Dialog.create("Fibril Tool");
  Dialog.addMessage("Choose the channel\n(Red/Green/Blue)");
  Dialog.addChoice("Channel for fibrils?\t", newArray("G", "R", "B"));
  Dialog.addChoice("Channel for drawing\t", newArray( "R", "G", "B", "No"));
  Dialog.addNumber("Multiply linelength by\t", 1);
  Dialog.addChoice("Input image format?\t", newArray(".png", ".tiff", ".tif"));
  Dialog.addChoice("Numbering of fibrils?\t", newArray("no", "yes"));
  Dialog.show();
  fib = Dialog.getChoice();
  drw = Dialog.getChoice();
  norm_constant = Dialog.getNumber(); // scaling factor for drawing of segments
  format =  Dialog.getChoice();
  numbering = Dialog.getChoice();
///compute MTs orientation for each projection of the directory
list = getFileList(dir);



for (j=0; j<list.length; j++){

	if(endsWith (list[j], format)){
		nom2=substring(list[j],0,indexOf(list[j], format));

		c=0;
		if( File.exists(dir+File.separator+nom2+".roi") ){
			c=1; //one ROI
		}else if (File.exists(dir+File.separator+nom2+"_RoiSet.zip") ){
			c=2; //two ROIs or more
		}

		if (c!=0){ //at least one ROI
	open( dir+File.separator+list[j] );
			if (c==1){
				open( dir+File.separator+nom2+".roi" );
				roiManager("add"); //a single ROI needs to be added to overlay
			}else if (c==2){
				open( dir+File.separator+nom2+"_RoiSet.zip" );
			}


//
id = getImageID(); 
title = getTitle();

getPixelSize(unit,pixelWidth,pixelHeight);
if (pixelWidth != pixelHeight) exit("Rectangular pixels!");
scale = pixelWidth;


//properties of selection
num++ ;
selectImage(id);

///ROI selection
  n = roiManager("count");
  for (cell=0; cell<n; cell++) {
      roiManager("select", cell);
      ROI_nm = Roi.getName(); 

getSelectionCoordinates(vertx, verty);
c = polygonCentre(vertx,verty);
c0s = c[0]*scale ;
c1s = c[1]*scale ;
getRawStatistics(area);
areas = area*scale*scale;
pr = 2;
sortie = title+"\t"+ROI_nm;
sortie = sortie+"\t"+d2s(c0s,pr)+"\t"+d2s(c1s,pr)+"\t"+d2s(areas,pr);

//extract fibril signal
selectImage(id);
run("Duplicate...", "title=Temp");
run("Crop"); 
getSelectionCoordinates(vertxloc, vertyloc);
if (fib == "R") setRGBWeights(1,0,0);
	else if (fib == "G") setRGBWeights(0,1,0); 
	else if (fib =="B") setRGBWeights(0,0,1);
	else exit("Fibril color undefined");
run("8-bit");


//compute x-gradient in "x"
selectWindow("Temp");
run("Duplicate...","title=x");
run("32-bit");
run("Translate...", "x=-0.5 y=0 interpolation=Bicubic");
run ("Duplicate...","title=x1");
run("Translate...", "x=1 y=0 interpolation=None");
imageCalculator("substract","x","x1");
selectWindow("x1");
close();

//compute y-gradient in "y"
selectWindow("Temp");
run ("Duplicate...","title=y");
run("32-bit");
run("Translate...", "x=0 y=-0.5 interpolation=Bicubic");
run ("Duplicate...","title=y1");
run("Translate...", "x=0 y=1 interpolation=None");
imageCalculator("substract","y","y1");
selectWindow("y1");
close();


//compute norm of gradient in "g"
selectWindow("x");
run("Duplicate...","title=g");
imageCalculator("multiply","g","x");
selectWindow("y");
run("Duplicate...","title=gp");
imageCalculator("multiply","gp","y");
imageCalculator("add","g","gp");
selectWindow("gp");
close();
selectWindow("g");
w = getWidth(); h = getHeight();
for (y=0; y<h; y++) {
	for (x=0; x<w; x++){
		setPixel(x, y, sqrt( getPixel(x, y)));
	}
}
//set the effect of the gradient to 1/255 when too low ; threshold = thresh
selectWindow("g");
for (y=0; y<h; y++) {
	for (x=0; x<w; x++){
		if (getPixel(x,y) < thresh) 
			setPixel(x, y, 255);
	}
}

//normalize "x" and "y" to components of normal
imageCalculator("divide","x","g");
imageCalculator("divide","y","g");


//compute nxx
selectWindow("x");
run("Duplicate...","title=nxx");
imageCalculator("multiply","nxx","x");
//compute nxy
selectWindow("x");
run("Duplicate...","title=nxy");
imageCalculator("multiply","nxy","y");
//compute nyy
selectWindow("y");
run("Duplicate...","title=nyy");
imageCalculator("multiply","nyy","y");

//closing
selectWindow("Temp");
close();
selectWindow("x");
close();
selectWindow("y");
close();
selectWindow("g");
close();

//averaging nematic tensor
selectWindow("nxx");
makeSelection("polygon",vertxloc,vertyloc);
getRawStatistics(area,xx);
selectWindow("nxx");
close();

selectWindow("nxy");
makeSelection("polygon",vertxloc,vertyloc);
getRawStatistics(area,xy);
selectWindow("nxy");
close();

selectWindow("nyy");
makeSelection("polygon",vertxloc,vertyloc);
getRawStatistics(area,yy);
selectWindow("nyy");
close();

//eigenvalues and eigenvector of texture tensor
m = (xx + yy) / 2;
d = (xx - yy) / 2;
v1 = m + sqrt(xy*xy + d*d);
v2 = m - sqrt(xy*xy + d*d);
//direction
tn = - atan((v2 - xx) / xy);
//score
scoren = abs((v1-v2) / 2 / m);

//output
tsn=tn*180/pi;
// nematic tensor tensor
sortie = sortie+"\t"+d2s(tsn,pr)+"\t"+d2s(scoren,2*pr);

//polygon coordinates
np = vertx.length;
for (i=0; i<np; i++){
xp = vertx[i]; yp = verty[i];
sortie = sortie+"\t"+d2s(xp,pr)+"\t"+d2s(yp,pr);
}



//
//print output
print(sortie);


//
//drawing of directions and cell contour
selectImage(id);
run("Add Selection...", "stroke=yellow width="+lwidth);


// drawing nematic tensor
if ( drw != "No" ) {
u1 = norm_constant*sqrt(area)*cos(tn)*scoren + c[0];
v1 = - norm_constant*sqrt(area)*sin(tn)*scoren + c[1];
u2 = - norm_constant*sqrt(area)*cos(tn)*scoren + c[0];
v2 =  norm_constant*sqrt(area)*sin(tn)*scoren + c[1];
if (drw == "R") stroke = "red";
	else if (drw == "G") stroke = "green"; 
	else if (drw =="B") stroke = "blue";
	else exit("Drawing color undefined");
makeLine(u1,v1,u2,v2);
run("Add Selection...", "stroke="+stroke+" width"+lwidth);
}


//print number at center
if (numbering == "yes") { makeText(ROI_nm,c[0],c[1]);
run("Add Selection...", "stroke="+stroke+" font="+fsize+" fill=none");
}


} //end of ROI selection


selectWindow("ROI Manager");
run("Close"); //Close the ROI manager
run("To ROI Manager");
roiManager("Save", dir+File.separator+nom2+"_RoiSet_MTs_contour.zip" );
roiManager("Show None"); 
roiManager("Show all without labels"); //OR roiManager("Show All");

 //Identification of ROI of cell contour
N = roiManager("count");

if (numbering == "yes") {
	step = 3;
}else{
	step = 2;
}

a1 = newArray(N/step);
for (i=0; i<a1.length; i++){
  a1[i] = i*step;
}


//Delete ROI(s)
roiManager("select", a1);
roiManager("Delete");

//Save microtubules orientation
roiManager("Save", dir+File.separator+nom2+"_RoiSet_MTs.zip" );
roiManager("Show all without labels");
run("Flatten"); //Flatten the overlay
saveAs("Tiff", dir+File.separator+nom2+"_MTs.tif" );
run("Close"); //Close the .tif image

newImage("Untitled", "RGB black", 1024, 1024, 1); //Project microtubules on a black background
roiManager("Show None"); 
roiManager("Show all without labels"); //OR roiManager("Show All");
run("Flatten");
saveAs("Tiff", dir+File.separator+nom2+"_MTs_blackBG.tif" );
close();
close();
close();
selectWindow("ROI Manager");
run("Close");

}//end of if( c!=0) i.e. at least one ROI

} //end of if endswith (list[j], format)

	} // end of for j=0; j<list.length; j++

//To save the final log and close it
selectWindow("Log");
saveAs("text", dir+"Full_Log.txt" );
run("Close");



// centroid of a polygon
function polygonCentre(x,y){
     n =x.length;
     area1 = 0;
     xc = 0; yc = 0;
     for (i=1; i<n; i++){
		  inc = x[i-1]*y[i] - x[i]*y[i-1];
         area1 += inc;
		  xc += (x[i-1]+x[i])*inc; 
		  yc += (y[i-1]+y[i])*inc;
     }
     inc = x[n-1]*y[0] - x[0]*y[n-1];
     area1 += inc;
     xc += (x[n-1]+x[0])*inc; 
     yc += (y[n-1]+y[0])*inc;    
     area1 *= 3;
     xc /= area1;
     yc /= area1;
     return newArray(xc,yc);
}



//distance between two points (x1,y1) et (x2,y2)
function distance(x1,y1,x2,y2) {
	return sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
   }


function printArray(a) {
  print("");
  for (i=0; i<a.length; i++)
      print(i+": "+a[i]);
}



