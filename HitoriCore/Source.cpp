#include <iostream>
#include "HitoriSolverCore.h"

double h(HitoriSolverCore::State) { return 0; }

int main() {
  HitoriSolverCore::HitoriSolver solver =
      HitoriSolverCore::HitoriSolver("18x18.txt");

  solver.PreProccess();

  HitoriSolverCore::State result = solver.GreedyBfs(solver.InitialState(), h);
  solver.PrintState(result);

  return 0;
}