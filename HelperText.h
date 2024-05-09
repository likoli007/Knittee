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
	const QString project3DText = QString("3D Project setup sucessfull! Please proceed in KnitGraph creation\n");
/*
Press 'Delete' when on a face with a constrained edge to delete the entire constraint. \n \
Press 'ENTER' to finish the current constraint loop, and press the 'Cancel Constraint Mode' button to finish";
*/
}


