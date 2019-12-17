
#include "HitoriSolverCore.h";
#include <iostream>;

double h(HitoriSolverCore::State) { return 0; }

int main() {
  HitoriSolverCore::HitoriSolver solver = HitoriSolverCore::HitoriSolver("input.txt");

  solver.PreProccess();

  HitoriSolverCore::State result = solver.GreedyBfs(HitoriSolverCore::State(15), h);

  solver.PrintState(result);

  return 0;
}