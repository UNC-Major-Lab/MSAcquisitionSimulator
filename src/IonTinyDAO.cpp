//
// Created by Dennis Goldfarb on 11/20/15.
//

#include "IonTinyDAO.h"

bool IonTinyDAO::less_rt_start(const IonTinyDAO &a, const IonTinyDAO &b) {
    return a.rt_start < b.rt_start;
}