// @File(Label="Select data path", style="directory") dir
print("\\Clear");
print("Working on directory: "+dir);
dir_p=dir+File.separator+"processed";
File.makeDirectory(dir_p);

number=newArray("1631");
match=".*still.*";
max="750";

setBatchMode(true);
for (n = 0; n < number.length; n++) {
	print("Selected number: "+number[n]);
	print("Matching: "+match);
	
	files=getFileList(dir);
	for(i=0; i<files.length; i++){
		if (startsWith(files[i],number[n])){
		if (matches(files[i], match)){
				print("Working on file: "+files[i]);
				print("Opening...");
			open(dir+File.separator+files[i]);
				print("Doing my thing...");
			run("Split Channels");
			
			selectWindow("C1-"+files[i]);
			name=substring(files[i], 0, files[i].length-4);
			saveAs("Tiff", dir_p+File.separator+name+"_PC"+".tif");
			
			selectWindow("C2-"+files[i]);
			run("Subtract Background...", "rolling=50 stack");
			setMinAndMax(0, max);
			name=substring(files[i], 0, files[i].length-4);
			saveAs("Tiff", dir_p+File.separator+name+"_mNG"+".tif");

			selectWindow("C3-"+files[i]);
			run("Subtract Background...", "rolling=50 stack");
			setMinAndMax(15, max/5);
			name=substring(files[i], 0, files[i].length-4);
			saveAs("Tiff", dir_p+File.separator+name+"_mScI"+".tif");
	
			run("Close All");
				print("Done!");
				print(" ");
		}
		}
		}
}
setBatchMode(false);
print("Done with all files!");