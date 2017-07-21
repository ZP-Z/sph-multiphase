#include "SolverGUI.h"

extern SolverGUI solverGUI;



int main ( int argc, char **argv )
{
	
	solverGUI.Initialize(argc, argv);

	solverGUI.SetupSolver();

	solverGUI.Run();

	solverGUI.Exit();
	
	return 0;
}

