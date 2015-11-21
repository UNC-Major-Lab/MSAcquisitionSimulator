//
// Created by Dennis Goldfarb on 11/20/15.
//

#ifndef MSACQUISITIONSIMULATOR_IONTINYDAO_H
#define MSACQUISITIONSIMULATOR_IONTINYDAO_H


#include "Peptide.h"

class IonTinyDAO {
private:

public:
    IonTinyDAO(int neutrons, int charge, double abundance, double rt,
           double rt_start, double rt_end, double mz, const Peptide* peptide) :
            neutrons(neutrons), charge(charge), abundance(abundance),
            rt(rt), rt_start(rt_start), rt_end(rt_end), peptide(peptide),
            mz(mz) {};

    int neutrons;
    int charge;
    double abundance;
    double rt;
    double rt_start;
    double rt_end;
    double mz;
    const Peptide* peptide;

    static bool less_rt_start(const IonTinyDAO& a, const IonTinyDAO& b);
};


#endif //MSACQUISITIONSIMULATOR_IONTINYDAO_H
