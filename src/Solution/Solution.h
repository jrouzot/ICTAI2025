//
// Created by Julien Rouzot on 20/05/24.
//

#pragma once

#include <vector>
#include "../Instance/Instance.h"

namespace dataflow {

    class FillRateSolution {
    public:
        explicit FillRateSolution(const std::vector <std::vector<Event>> &fillRates);

    private:
        std::vector <std::vector<Event>> mFillRates;
    };

    class TransferSolution {
    public:
        TransferSolution();
        TransferSolution(
                double rmax,
                std::vector <std::vector<double>> &handover,
                double resolutionTime
                );

        void Print();
        double getRmax() { return mRmax; };
        double getResolutionTime() { return mResolutionTime; };
        bool isSolution() { return mExists; };
        const std::vector <std::vector <double>> &getHandovers() { return mHandovers; };


    private:
        bool mExists;
        double mRmax;
        double mResolutionTime;
        std::vector <std::vector <double>> mHandovers;
    };

}
