#include "main.hpp"

int32_t main(int argc, char** argv) {
    string fileName = argv[1];
    int b = strtol(argv[2], NULL, 10);

    Solver G1(fileName, b);
    vector<int> solution = G1.solve();

    for (auto node : solution)
        cout << node << " ";
    cout << endl;

    Checker checkG1(fileName);
    if (checkG1.check(solution, 1))
        cout << BOLD(FGRN("Verified!")) << endl;
    else
        cout << BOLD(FRED("NOT Verified!")) << endl;
}
