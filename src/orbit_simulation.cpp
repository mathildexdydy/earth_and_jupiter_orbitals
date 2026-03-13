/*
PY2105 : Introduction to Computational Physics
Final Exam
Sun – Earth – Jupiter simulation
RK4 integration
*/

#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>

using namespace std;

/* ---------------------- CONSTANTS ---------------------- */

const double G = 6.67430e-11;

const double m_sun = 1.989e30;
const double m_earth = 5.972e24;
const double m_jupiter = 1.898e27; // 318 Earth masses

const double AU = 1.496e11;

const double r_earth = AU;
const double r_jupiter = 5.2 * AU;

const double T = 374016960.0; // Jupiter orbital period (s)

int time_steps[3] = {3600, 86400, 604800}; // 1h, 1d, 1week


/* ---------------------- ACCELERATIONS ---------------------- */

void compute_accelerations(
    double xe, double ye,
    double xj, double yj,
    double &axe, double &aye,
    double &axj, double &ayj)
{
    double re = sqrt(xe*xe + ye*ye);
    double rj = sqrt(xj*xj + yj*yj);

    double dx = xj - xe;
    double dy = yj - ye;
    double rej = sqrt(dx*dx + dy*dy);

    /* Sun forces */

    axe = -G*m_sun*xe/pow(re,3);
    aye = -G*m_sun*ye/pow(re,3);

    axj = -G*m_sun*xj/pow(rj,3);
    ayj = -G*m_sun*yj/pow(rj,3);

    /* Earth–Jupiter interaction */

    axe += G*m_jupiter*dx/pow(rej,3);
    aye += G*m_jupiter*dy/pow(rej,3);

    axj -= G*m_earth*dx/pow(rej,3);
    ayj -= G*m_earth*dy/pow(rej,3);
}


/* ---------------------- RK4 INTEGRATOR ---------------------- */

void rk4_step(
    double &xe, double &ye, double &vxe, double &vye,
    double &xj, double &yj, double &vxj, double &vyj,
    double dt)
{

    double axe1, aye1, axj1, ayj1;
    double axe2, aye2, axj2, ayj2;
    double axe3, aye3, axj3, ayj3;
    double axe4, aye4, axj4, ayj4;

    compute_accelerations(xe,ye,xj,yj,axe1,aye1,axj1,ayj1);

    double k1xe = vxe*dt;
    double k1ye = vye*dt;
    double k1vxe = axe1*dt;
    double k1vye = aye1*dt;

    double k1xj = vxj*dt;
    double k1yj = vyj*dt;
    double k1vxj = axj1*dt;
    double k1vyj = ayj1*dt;


    compute_accelerations(
        xe + 0.5*k1xe,
        ye + 0.5*k1ye,
        xj + 0.5*k1xj,
        yj + 0.5*k1yj,
        axe2,aye2,axj2,ayj2);

    double k2xe = (vxe + 0.5*k1vxe)*dt;
    double k2ye = (vye + 0.5*k1vye)*dt;
    double k2vxe = axe2*dt;
    double k2vye = aye2*dt;

    double k2xj = (vxj + 0.5*k1vxj)*dt;
    double k2yj = (vyj + 0.5*k1vyj)*dt;
    double k2vxj = axj2*dt;
    double k2vyj = ayj2*dt;


    compute_accelerations(
        xe + 0.5*k2xe,
        ye + 0.5*k2ye,
        xj + 0.5*k2xj,
        yj + 0.5*k2yj,
        axe3,aye3,axj3,ayj3);

    double k3xe = (vxe + 0.5*k2vxe)*dt;
    double k3ye = (vye + 0.5*k2vye)*dt;
    double k3vxe = axe3*dt;
    double k3vye = aye3*dt;

    double k3xj = (vxj + 0.5*k2vxj)*dt;
    double k3yj = (vyj + 0.5*k2vyj)*dt;
    double k3vxj = axj3*dt;
    double k3vyj = ayj3*dt;


    compute_accelerations(
        xe + k3xe,
        ye + k3ye,
        xj + k3xj,
        yj + k3yj,
        axe4,aye4,axj4,ayj4);

    double k4xe = (vxe + k3vxe)*dt;
    double k4ye = (vye + k3vye)*dt;
    double k4vxe = axe4*dt;
    double k4vye = aye4*dt;

    double k4xj = (vxj + k3vxj)*dt;
    double k4yj = (vyj + k3vyj)*dt;
    double k4vxj = axj4*dt;
    double k4vyj = ayj4*dt;


    xe += (k1xe + 2*k2xe + 2*k3xe + k4xe)/6.0;
    ye += (k1ye + 2*k2ye + 2*k3ye + k4ye)/6.0;

    vxe += (k1vxe + 2*k2vxe + 2*k3vxe + k4vxe)/6.0;
    vye += (k1vye + 2*k2vye + 2*k3vye + k4vye)/6.0;


    xj += (k1xj + 2*k2xj + 2*k3xj + k4xj)/6.0;
    yj += (k1yj + 2*k2yj + 2*k3yj + k4yj)/6.0;

    vxj += (k1vxj + 2*k2vxj + 2*k3vxj + k4vxj)/6.0;
    vyj += (k1vyj + 2*k2vyj + 2*k3vyj + k4vyj)/6.0;
}


/* ---------------------- ENERGY ---------------------- */

double total_energy(
    double xe,double ye,double vxe,double vye,
    double xj,double yj,double vxj,double vyj)
{

    double re = sqrt(xe*xe + ye*ye);
    double rj = sqrt(xj*xj + yj*yj);

    double rej = sqrt(pow(xe-xj,2) + pow(ye-yj,2));

    double KE =
        0.5*m_earth*(vxe*vxe + vye*vye) +
        0.5*m_jupiter*(vxj*vxj + vyj*vyj);

    double PE =
        -G*m_sun*m_earth/re
        -G*m_sun*m_jupiter/rj
        -G*m_earth*m_jupiter/rej;

    return KE + PE;
}


/* ---------------------- MAIN PROGRAM ---------------------- */

int main()
{

for(int k=0;k<3;k++)
{

    double dt = time_steps[k];
    int n_steps = T/dt;

    double xe = r_earth;
    double ye = 0;

    double xj = -r_jupiter;
    double yj = 0;


    double vxe = 0;
    double vye = sqrt(G*m_sun/r_earth);

    double vxj = 0;
    double vyj = sqrt(G*m_sun/r_jupiter);


    ofstream orbit("orbit_"+to_string((int)dt)+".csv");
    ofstream energy("energy_"+to_string((int)dt)+".csv");

    orbit<<"time_days,xe,ye,xj,yj\n";
    energy<<"time_days,E\n";


    for(int i=0;i<n_steps;i++)
    {

        double time_days = (i*dt)/86400.0;

        orbit<<time_days<<","
             <<xe<<","<<ye<<","
             <<xj<<","<<yj<<"\n";

        double E = total_energy(xe,ye,vxe,vye,xj,yj,vxj,vyj);

        energy<<time_days<<","<<E<<"\n";


        rk4_step(
            xe,ye,vxe,vye,
            xj,yj,vxj,vyj,
            dt);
    }

    orbit.close();
    energy.close();
}

return 0;
}
