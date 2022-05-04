cases = newArray('PSP1_2a1_1', 'PSP1_2a1_2', 'PSP1_2a1_3','PSP1_2a1_4', 'PSP1_2a1_5', 
				 'PSP2_2a1_1', 'PSP2_2a1_2', 'PSP2_2a1_3','PSP2_2a1_4', 'PSP2_2a1_5', 'PSP2_2a1_6', 
				 'PSP3_2a1_2', 'PSP3_2a1_3','PSP3_2a1_4', 'PSP3_2a1_5', 'PSP3_2a1_6', 'PSP3_2a1_7', 'PSP3_2a1_8'
				 'PSP4_2a1_1', 'PSP4_2a1_2', 'PSP4_2a1_3','PSP4_2a1_4', 'PSP4_2a1_5', 'PSP4_2a1_6',

				'CBD1_2a1_1', 'CBD1_2a1_2', 'CBD1_2a1_3', 'CBD1_2a1_4', 'CBD1_2a1_5', 'CBD1_2a1_6', 
				'CBD2_2a1_1', 'CBD2_2a1_2', 'CBD2_2a1_3', 'CBD2_2a1_3', 'CBD2_2a1_5', 'CBD2_2a1_6',
				'CBD3_2a1_1', 'CBD3_2a1_2', 'CBD3_2a1_3', 'CBD3_2a1_4', 'CBD3_2a1_5', 'CBD3_2a1_6',
				'CBD4_2a1_1', 'CBD4_2a1_2', 'CBD4_2a1_3', 'CBD4_2a1_4', 'CBD4_2a1_5', 'CBD4_2a1_6'
				);

for (i=0; i <lengthOf(cases); ++i) {
	case_no = cases[i];
	open("/home/nes/projs/scATACseq_MLL/Validation_IF/CTSD/2nd_batch/merged/Composite_" + case_no + ".tif");
	run("Z Project...", "projection=[Standard Deviation]");
	Stack.setDisplayMode("color");
	run("Window/Level...");
	Stack.setChannel(1);
	setMinAndMax(12.62, 121.72);
	run("Grays");
	Stack.setChannel(2);
	setMinAndMax("14.85", "128.08");
	run("Cyan");
	Stack.setChannel(3);
	setMinAndMax("8.84", "77.66");
	run("Magenta");
	Stack.setChannel(4);
	setMinAndMax(2.91, 74.20);
	run("Yellow");
	Stack.setDisplayMode("composite");
	saveAs("Tiff", "/home/nes/projs/scATACseq_MLL/Validation_IF/CTSD/2nd_batch/merged/Composite_mpi_" + case_no + ".tif");
	close();
	close();
	
}
