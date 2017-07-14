#include "SolverGUI.h"

extern SolverGUI solverGUI;



int main ( int argc, char **argv )
{
	
	solverGUI.Initialize(argc, argv);

	solverGUI.SetupSolver();

	solverGUI.Run();

	solverGUI.Exit();
	
	/*cTime st;
	st.printTime();

	getchar();

	cTime st2;
	printf("%d\n", st2-st);
	st2.printTime();*/
	return 0;
}

