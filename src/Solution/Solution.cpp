//
// Created by Julien Rouzot on 20/05/24.
//

#include "Solution.h"

namespace dataflow {

    FillRateSolution::FillRateSolution(const std::vector <std::vector<Event>> &fillRates) : mFillRates(fillRates) {}

    TransferSolution::TransferSolution() : mExists(false) {}
    TransferSolution::TransferSolution(double rmax,
                                       std::vector<std::vector<double>> &handover,
                                       double resolutionTime) :
            mExists(true), mRmax(rmax), mHandovers(handover), mResolutionTime(resolutionTime) {}

    void TransferSolution::Print() {

        if(mExists) {
            std::cout << "***\nSolution:\n";

            std::cout << "Objective: " << mRmax << std::endl;

            for(int w{0}; w < mHandovers.size(); ++w) {
                std::cout << "Handovers for window " << w << ": ";
                for(int b{0}; b < mHandovers[w].size(); ++b) {
                    std::cout << mHandovers[w][b] << " ";
                }
                std::cout << std::endl;
            }

            std::cout << "***\n";
        } else {
            std::cout << "***\nNo solution\n";
        }
    }
}
