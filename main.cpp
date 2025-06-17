#include <iostream>
#include <string>
#include <algorithm>
#include "src/Algorithm/Algorithm.h"
#include "src/CplexFlowDataTransferModel/CplexFlowDataTransferModel.h"
#include "src/DownlinkCount/DownlinkCount.h"
#include "src/Instance/Instance.h"
#include "src/Logger/Logger.h"
#include "src/Options/Options.h"
#include "src/Simulator/Simulator.h"

#include "tclap/CmdLine.h"

int main(int argc, char *argv[]) {
    try {
        TCLAP::CmdLine cmd("Command description message", ' ', "0.9");
        TCLAP::UnlabeledValueArg<std::string> instance_file("f", "Path to the instance file", true, "", "string");
        cmd.add(instance_file);
        TCLAP::ValueArg<int> verbosityCMD("v", "log",
                                       "log level (0:ERROR,1:WARNING,2:DEBUG)", false, 0, "int");
        TCLAP::ValueArg<double> timeLimitCMD("t", "time",
                                       "time limit (seconds)", false, 600.0, "double");

        cmd.add(verbosityCMD);
        cmd.add(timeLimitCMD);
        cmd.parse(argc, argv);
        logger::Logger logger((logger::Verbosity) verbosityCMD.getValue());
        double timeLimit = (double) timeLimitCMD.getValue();
        dataflow::Instance instance(instance_file.getValue(), logger);

        simulation::Algorithm algorithm(instance, logger);
        double rmax = algorithm.solve(timeLimit);

    } catch(TCLAP::ArgException &e) {
        std::cerr << "[ARG ERROR]: " << e.error() << " for arg " << e.argId() << std::endl;
    }
}
