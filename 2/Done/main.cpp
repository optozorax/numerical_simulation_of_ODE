#include "pump.h"

int main() {

	pump pump1, pump2, pump3;

	pump1.solveDE_Euler_Explicit();
	pump1.outputResults("1.txt");

	pump2.solveDE_Euler_Nonexplicit();
	pump2.outputResults("2.txt");

	pump3.solveDE_Runge_Kutta4();
	pump3.outputResults("3.txt");
}