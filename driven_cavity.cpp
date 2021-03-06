#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <omp.h>
#include "positions.hpp"
#include "import.hpp"
#include "velocity.hpp"
#include "boundary.hpp"
#include "export.hpp"
#include "util.hpp"
#include "poiss.hpp"
int main()
{
    std::vector<double> v;    
    std::vector <std::vector<double>> Rpu, Rnu, Rpv, Rnv, up, vp, Vm, C,Vp;
    std::vector<int> pos;
    double Pr, Ray, deltatc, deltatd, deltat, runtime, total_time(0),M(0), MM(0),diff;
    int n(0);
    double k(0.01);
    v = readfiledat();
    Pr = v[2];
    Ray = v[5];
    runtime = v[4];
    Poisson poiss(v[3]);
    Solver solv(v[3]);
    positions mesh(v[0], v[1]);
    Velocity Xvel(mesh.get_m() + 2, mesh.get_n() + 1);
    Velocity Yvel(mesh.get_m() + 1, mesh.get_n() + 2);
    
    Pressure P(mesh.get_m() + 2, mesh.get_n() + 2);
    Temperature T(mesh.get_m() + 2, mesh.get_n() + 2);
    
    std::vector<Velocity> V;
    
    
    V.push_back(Xvel); V.push_back(Yvel); 
    T.set_Tn(mesh, V, T);
    
    V[0].set_Vp(mesh, V[0]);  V[1].set_Vp(mesh, V[1]);
    V[0].set_Vpc(mesh, V[0],V);  V[1].set_Vpc(mesh, V[1], V);
    Rpu = get_Ru(V, mesh, Pr);
    
    Rpv = get_Rv(V, mesh, T, Pr, Ray);
    
    set_boundary(V, P, T);

    V[0].set_Vp(mesh, V[0]); V[1].set_Vp(mesh, V[1]);
    V[0].set_Vpc(mesh, V[0], V); V[1].set_Vpc(mesh, V[1], V);
    deltatc = get_deltatc(V, mesh);
    deltatd = get_deltatd(Pr, mesh);
    //std::cout << "deltatd= " << deltatd << " deltatc= " << deltatc << std::endl;
    deltat = std::min(deltatd,deltatc);
    Rnu = get_Ru(V, mesh, Pr);
    Rnv = get_Rv(V, mesh, T, Pr, Ray);
    C = substract(Rnu, Rpu);
    diff = get_max(C);
    pos = get_maxpos(C);
    std::cout << diff << ", " << pos[0] << ", " << pos[1] << ", " << total_time << ", ";
    up = get_up(V, Rnu, Rpu, deltat);
    vp = get_vp(V, Rnv, Rpv, deltat);
    set_boundary(up, vp);
    poiss.set_P(P, up, vp, mesh,  deltat);    
    poiss.set_V(P, mesh, V, up, vp, deltat);
    
    // The velocity at this point is known, with this information, the temperature field can be computed
    solv.get_T(T, mesh, V, deltat);
    T.set_Tn(mesh, V, T);
    total_time += deltat;
    std::string save = "";
    std::ostringstream strs;
    strs << total_time;
    save = "Results/Rnv" + strs.str() + ".out";
    exportarMatriu(Rnv, save.c_str());
    save = "Results/Rnu" + strs.str() + ".out";
    exportarMatriu(Rnu, save.c_str());

    diff = v[3] + 10;
    //while(total_time < runtime)
    while(diff > v[3])
    //while(n<0)
    {
        n+=1;
        V[0].set_Vp(mesh, V[0]); V[1].set_Vp(mesh, V[1]);
        V[0].set_Vpc(mesh, V[0], V); V[1].set_Vpc(mesh, V[1], V);
        copy_matrix(Rpu, Rnu); copy_matrix(Rpv, Rnv); copy_matrix(Vp,V[0].get_V());
        deltatc = get_deltatc(V, mesh);
        deltatd = get_deltatd(Pr, mesh);
        deltat = std::min(deltatd,deltatc);
        Rnu = get_Ru(V, mesh, Pr);
        Rnv = get_Rv(V, mesh, T, Pr, Ray);
        //C = substract(Rnu, Rpu);
        //diff = get_max(C);
        //pos = get_maxpos(C);
        std::cout << diff << ", " << pos[0] << ", " << pos[1] << ", " << total_time << ", ";
        up = get_up(V, Rnu, Rpu, deltat);
        vp = get_vp(V, Rnv, Rpv, deltat);
        set_boundary(up, vp);
        poiss.set_P(P, up, vp, mesh,  deltat);  
        //std::cout << ", " << total_time << std::endl;  
        poiss.set_V(P, mesh, V, up, vp, deltat);
        C = substract(Vp, V[0].get_V());
        diff = get_max(C);
        pos = get_maxpos(C);
        solv.get_T(T, mesh, V, deltat);
        T.set_Tn(mesh, V, T);
        total_time += deltat;
        if (total_time > k * runtime)
        {
            std::ostringstream strs;
            k+= 0.01;
            V[0].set_Vp(mesh, V[0]); V[1].set_Vp(mesh, V[1]);
            V[0].set_Vpc(mesh, V[0], V); V[1].set_Vpc(mesh, V[1], V);
            strs << total_time;
            save = "Results/vel_x" + strs.str() + ".out";
            exportarMatriu(V[0].get_V(),save.c_str());
            save = "Results/vel_y" + strs.str() + ".out";
            exportarMatriu(V[1].get_V(),save.c_str());
            save = "Results/Pressure" + strs.str() + ".out";
            exportarMatriu(P.get_P(), save.c_str());
            save = "Results/vel_Xp" + strs.str() + ".out";
            exportarMatriu(V[0].get_Vp(), save.c_str());
            save = "Results/vel_Yp" + strs.str() + ".out";
            exportarMatriu(V[1].get_Vp(), save.c_str());
            save = "Results/Temp" + strs.str() + ".out";
            exportarMatriu(T.get_T(), save.c_str());

        }
        M = 0;
        for (int i = 1; i <= mesh.get_m(); i++)
        {
            for (int j = 1; j <= mesh.get_n(); j++)
            {
                MM = (V[1].get_V(i, j) - V[1].get_V(i - 1, j)) * (mesh.get_Dxpr()[1] + mesh.get_Dxpl()[1]) + (V[0].get_V(i, j) - V[0].get_V(i, j - 1)) * (mesh.get_Dypu()[1] + mesh.get_Dypd()[1]) ;
                MM = std::abs(MM);
                M = std::max(M,MM);
            }
        }
        //std::cout << V[0].get_V(mesh.get_m() - 2, std::floor(mesh.get_n() / 2)) << ", " << total_time << std::endl;
        //std::cout << M << std::endl;
        //std::cin >> a;
    }





    char name3[] = "vel_xp.out";
    char name4[] = "vel_yp.out";

    char name6[] = "up.out";
    char name7[] = "vp.out";
    char name8[] = "Rnu.out";
    char name9[] = "Rnv.out";
    exportarMatriu(V[0].get_Vp(),name3);
    exportarMatriu(V[1].get_Vp(),name4);

    exportarMatriu(up, name6);
    exportarMatriu(vp, name7);
    exportarMatriu(Rnu, name8);
    exportarMatriu(Rnv, name9);
    std::cout << "End!" << std::endl;
    return 0;
}