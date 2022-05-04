cases = newArray('PSP1_2a1', 'PSP1_2a2', 'PSP1_2a3','PSP1_2a4', 'PSP1_2a5', 
				 'PSP2_2a1', 'PSP2_2a2', 'PSP2_2a3','PSP2_2a4', 'PSP2_2a5', 
				 'PSP3_2a1', 'PSP3_2a2', 'PSP3_2a3','PSP3_2a4', 
				 'PSP4_2a1', 'PSP4_2a2', 'PSP4_2a3','PSP4_2a4', 'PSP4_2a5', 

				'CBD1_1a1', 'CBD1_1a2', 'CBD1_1a3', 'CBD1_1a4','CBD2_1a1', 'CBD2_1a2', 'CBD2_1a3',
				'CBD3_1a1', 'CBD3_1a2', 'CBD3_1a3', 'CBD3_1a4','CBD4_1a1', 'CBD4_1a2', 'CBD4_1a3', 'CBD4_1a4');

for (i=0; i <lengthOf(cases); ++i) {
	case_no = cases[i];
	open("/home/nes/projs/scATACseq_MLL/Validation_IF/MAP3K8/merged/" + case_no + ".tif");
	run("Z Project...", "projection=[Standard Deviation]");
	Stack.setDisplayMode("color");
	run("Window/Level...");
	Stack.setChannel(1);
	setMinAndMax("3.51", "36.24");
	run("Grays");
	Stack.setChannel(2);
	setMinAndMax("14.85", "128.08");
	run("Cyan");
	Stack.setChannel(3);
	setMinAndMax("8.84", "77.66");
	run("Magenta");
	Stack.setChannel(4);
	setMinAndMax("10.67", "74.69");
	run("Yellow");
	Stack.setDisplayMode("composite");
	saveAs("Tiff", "/home/nes/projs/scATACseq_MLL/Validation_IF/MAP3K8/merged/Composite_mpi_" + case_no + ".tif");
	close();
	close();
	
}
