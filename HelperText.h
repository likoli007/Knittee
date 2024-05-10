#pragma once

#include <qstring.h>


namespace HelperText
{

	const QString constraintText = QString("Constraints mode enabled.\nPress 'C' while hovering over the mesh to add a constraint.\n") 
		.append("Press 'Ctrl-Z' to undo the added constraint edges in the current constraint.\n")
		.append("Press 'Delete' when on a face with a constrained edge to delete the entire constraint. \n")
		.append("Press '+'/'-' when hovering over a face with an associated constraint to change its time value.\n")
		.append("Press 'ENTER' to finish the current constraint loop, and press the 'Cancel Constraint Mode' button to finish");
	const QString constraintEndText = QString("Constraints mode disabled.")
		.append("Press 'Constraints Mode' to add more constraints or Press 'Interpolate' to interpolate the mesh");
	const QString interpolateText = QString("Interpolating mesh, this may take a while...");
	const QString project3DText = QString("3D Project setup successful! Please proceed in KnitGraph creation.\n")
		.append("If this is your first time using this program, please press 'Help' to read the instructions.");

	const QString resetText = QString("Project reset to original mesh, feel free to experiment with new parameters and start again\n");

	const QString helpWindowText3D = QString("Welcome to Knittee 2D / 3D knitted objects modelling software!\n")
		.append("Instructions for 3D modelling:\n")
		.append("\nStep 1: Add at least 2 constraint chains to your mesh using the 'Constraints Mode' button\n")
		.append("Step 2: Tinker with the stitch parameters, usually it is best to leave the stitch width/height as is and change the unit length\n")
		.append("\t Caution: do not select too big of a unit length, as this causes the algorithm to run significantly slower\n")
		.append("Step 3: Press the 'Interpolate' button to embed your constraints into your mesh and linearly interpolate between them\n")
		.append("Step 4: Press the 'Step' button to step through the process of creating a KnitGraph from your mesh\n")
		.append("\t Note: you can also use the 'Auto Step' button to step through many iterations at once\n")
		.append("Step 5: Once the KnitGraph is created, press the 'Trace' button to trace the KnitGraph and create basic stitch information\n")
		.append("Step 6: If you are unsatisfied with your result, press the 'Reset' button to reset the algorithm and select different parameters\n")
		.append("Step 7: Once satisfied, press the 'Generate Knitout' button to generate machine independent Knitout instructions\n")
		.append("Step 8: Once Knitout instructions are generated, you can press 'Export Instructions' to convert Knitout into machine-specific language\n")
		.append("\nNote: The project saves itself automatically once certain steps are finished, when it is loaded again it will keep all built data structures");

	const QString helpWindowText2D = QString("Welcome to Knittee 2D / 3D knitted objects modelling software!\n")
		.append("Instructions for 2D modelling:\n");

	const QString aboutWindowText = QString("Knittee is a 2D sheet and 3D object to Knit instruction CAD software\n") +
		QString("Designed by Alojz Holubek as part of his bachelor's thesis at Zhejiang University\n");


}


