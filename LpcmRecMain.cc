#include "LpcmRec/LpcRun.hh"

#include <string>

int main(int argc, char **argv) {

    std::string parameterFileName = "lpcParameters.dat";
    if (argc > 1) {parameterFileName = argv[1];}

    LpcRun run(parameterFileName);
    run.run();

    return 0;

}
