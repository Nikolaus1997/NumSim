#pragma once

#include <iostream>
#include <array>

/** All settings that parametrize a simulation run.
 */
struct Settings
{
  std::array<int,3> nCells{10,10,10};          //< number of cells in x and y direction
  std::array<double,3> physicalSize{1.,1.,1.}; //< physical size of the domain
  double re =50.;                 //< reynolds number
  double endTime = 10.0;             //< end time of the simulation
  double tau = .6;                  //< safety factor for time step width
  double maximumDt = 10.0;            //< maximum time step width
  double L_lbm = 10.;
  std::array<double,2> g{0., 0.};    //< external forces

  bool useDonorCell = false;         //< if the donor cell scheme schould be used
  double alpha = 0.5;                //< factor for donor-cell scheme
  int deltawrite_ = 5;
  double rhoRight = 1. ;
  bool rightBcPressure = false;
  std::array<double,3> dirichletBcBottom{0.,0.,0.};  //< prescribed values of u,v at bottom of domain
  std::array<double,3> dirichletBcTop{0.1,0.,0.};     //< prescribed values of u,v at top of domain
  std::array<double,3> dirichletBcLeft{0.,0.,0.};    //< prescribed values of u,v at left of domain
  std::array<double,3> dirichletBcRight{0.,0.,0.};   //< prescribed values of u,v at right of domain
  std::array<double,3> dirichletBcFront{0.,0.,0.};   //< prescribed values of u,v at right of domain
  std::array<double,3> dirichletBcBack{0.,0.,0.};   //< prescribed values of u,v at right of domain


  std::string pressureSolver = "SOR";      //< which pressure solver to use, "GaussSeidel" or "SOR"
  double omega = 1.0;                //< overrelaxation factor
  double epsilon = 1e-5;             //< tolerance for the residual in the pressure solver
  int maximumNumberOfIterations = 1e5;    //< maximum number of iterations in the solver

  //! parse a text file with settings, each line contains "<parameterName> = <value>"
  void loadFromFile(std::string filename);

  //! output all settings to console
  void printSettings();
};