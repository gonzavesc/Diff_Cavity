#include <iostream>
#include <fstream>
#include <vector>
#include "import.hpp"

#ifndef INCLUDE_POS
#define INCLUDE_POS
class positions{
    private:
        double dx,dy, n, m;
        std::vector<double> Xp, Yp, Xu, Yv, Dxpr, Dxpl, Dypu, Dypd, Dxur, Dxul, Dyvu, Dyvd;

    public:
        positions(const double& a, const double& b);
        std::vector<double> get_Xp();
        std::vector<double> get_Yp();
        std::vector<double> get_Yv();
        std::vector<double> get_Xu();
        std::vector<double> get_Dxpr();
        std::vector<double> get_Dxpl();
        std::vector<double> get_Dypu();
        std::vector<double> get_Dypd();
        std::vector<double> get_Dxur();
        std::vector<double> get_Dxul();
        std::vector<double> get_Dyvu();
        std::vector<double> get_Dyvd();
        double get_Dxpr(const int& a);
        double get_Dxpl(const int& a);
        double get_Dypu(const int& a);
        double get_Dypd(const int& a);
        int get_m();
        int get_n();

};
#endif