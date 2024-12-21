//                                Allen-Cahn equation
//
// Compile with: make allencahn
//
// Sample runs:  allencahn
//               allencahn -m ../data/square-disc.mesh
//               allencahn -m ../data/star.mesh
//               allencahn -m ../data/ref-square.mesh
//
//               allencahn -m ../data/square-disc-surf.mesh
//               allencahn -m ../data/star-surf.mesh
//               allencahn -m ../data/mobius-strip.mesh
//
//               allencahn -m ../data/fichera.mesh
//               allencahn -m ../data/ref-prism.mesh
//
// Note:         For visualization, a glvis server should be active when running the code.
//               As the initial condition is randomly generated, 
//               it is recommended to run the code with the same mesh at least twice.


#include "mfem.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;

// Computes J(u) = the value of the functional J at GridFunction u
double J(FiniteElementSpace *fes, GridFunction &u, double eps);

int main(int argc, char *argv[])
{
   // 1. Parse command-line options.
   const char *mesh_file = "../data/square-disc.mesh";
   bool visualization = true;

   double eps = 0.04;   // eps = epsilon from the model
   double step0 = 0.05; // Initial step guess. We take step0 of order eps
   double beta = 0.5;   // Armijo parameters
   double gamma = 0.5;  //
   double tol = 0.05;   // Exit tolerance in gradient descent
   int vis_step = 10;   // Visualization updates every vis_step number of iterations

   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh",
                  "Mesh file to use.");
   args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                  "--no-visualization",
                  "Enable or disable GLVis visualization.");
   args.Parse();
   if (!args.Good())
   {
      args.PrintUsage(cout);
      return 1;
   }
   args.PrintOptions(cout);

   // 3. Read the mesh from the given mesh file.
   Mesh mesh(mesh_file, 1, 1);
   int dim = mesh.Dimension();

   // 4. Refine the mesh to increase the resolution. We do 'ref_levels'
   //    of uniform refinement. We choose 'ref_levels' to be the largest number
   //    that gives a final mesh with no more than 15,000 elements.    
   int ref_levels = (int)floor(log(15000. / mesh.GetNE()) / log(2.) / dim);
   for (int l = 0; l < ref_levels; l++) {
      mesh.UniformRefinement();
   }

   // 5. Define a finite element space on the mesh (order = 1) and the
   //    solution GridFunction u_prev
   FiniteElementCollection *fec;
   fec = new H1_FECollection(1, dim);
   FiniteElementSpace fespace(&mesh, fec);
   cout << "Number of finite element unknowns: "
        << fespace.GetTrueVSize() << endl;
   Array<int> ess_tdof_list; // this list remains empty for Neumann b.c.
   GridFunction u_prev(&fespace);

   // 6. Start with random initial data in the range [-0.95, 0.95]
   u_prev.Randomize();
   u_prev -= 0.5;
   u_prev *= 1.9;

   // 7. Save the mesh and the initial data. This output can be viewed later
   //    using GLVis: "glvis -m refined.mesh -g sol_init.gf".
   ofstream mesh_ofs("refined.mesh");
   mesh_ofs.precision(8);
   mesh.Print(mesh_ofs);
   ofstream sol_ofs("sol_init.gf");
   sol_ofs.precision(8);
   u_prev.Save(sol_ofs);
   
   // 8. Set up Armijo-type gradient descent step sizes and pre-assemble
   // the corresponding biliear forms
   int max_num_steps = 6;
   double step[max_num_steps];
   step[0] = step0;
   for (int l = 1; l < max_num_steps; l++)
      step[l] = step[l-1] * beta;
   BilinearForm *a[max_num_steps];
   ConstantCoefficient one(1.0);
   for (int l = 0; l < max_num_steps; l++) {
      a[l] = new BilinearForm(&fespace);
      a[l]->AddDomainIntegrator(new MassIntegrator(one));
      ConstantCoefficient steptimeseps(step[l]*eps);
      a[l]->AddDomainIntegrator(new DiffusionIntegrator(steptimeseps));
      a[l]->Assemble();
   }

   // 9. Set up visualization
   socketstream sout;
   if (visualization) {
      char vishost[] = "localhost";
      int  visport   = 19916;
      sout.open(vishost, visport);
      if (!sout) {
         cout << "Unable to connect to GLVis server at "
              << vishost << ':' << visport << endl;
         visualization = false;
         cout << "GLVis visualization disabled.\n";
      }
      else {
         sout.precision(8);
         sout << "solution\n" << mesh << u_prev;
         sout << "pause\n";
         sout << flush;
         cout << "GLVis visualization paused."
              << " Press space (in the GLVis window) to resume it.\n";
      }
   }

   // 10. Preparation before starting the gradient descent
   int iter_counter = 0;
   GridFunction u_next(&fespace);
   double J_prev = J(&fespace, u_prev, eps);
   double J_next;
   double DJu_Norm;
   SparseMatrix A;
   Vector B, X;

   // 11. Running the gradient descent
   do {
      int l = 0; // Armijo step counter
      do {
         // Descending with step=step[l]. Initializing the right hand side
         LinearForm b(&fespace);
         GridFunctionCoefficient u_gfc(&u_prev);
         b.AddDomainIntegrator(new DomainLFIntegrator(u_gfc));
         PowerCoefficient u3(u_gfc, 3);
         ProductCoefficient RHS1(-step[l]/eps, u3);
         ProductCoefficient RHS2(step[l]/eps, u_gfc);
         b.AddDomainIntegrator(new DomainLFIntegrator(RHS1));
         b.AddDomainIntegrator(new DomainLFIntegrator(RHS2));
         b.Assemble();
         // Solving the linear system, the result is u_next
         a[l]->FormLinearSystem(ess_tdof_list, u_next, b, A, X, B);
         GSSmoother M(A);
         PCG(A, M, B, X, 0, 200, 1e-12, 0.0);
         a[l]->RecoverFEMSolution(X, b, u_next);
         J_next = J(&fespace, u_next, eps);
         DJu_Norm = u_next.ComputeL2Error(u_gfc) / step[l]; // Norm of DJ(u_prev)
         l++;
      } while (J_next > J_prev - gamma * step[l-1] * DJu_Norm*DJu_Norm &&
               l < max_num_steps); // Checking Armijo condition
      
      // Printing functional values at the intermediate solution.
      // Sending visualization of the intermediate solution
      if (iter_counter % vis_step == 0 || !(DJu_Norm > tol)) {
         cout << "This is iteration " << iter_counter << endl;
         cout << "J(u) = " << J_prev << endl;
         cout << "|| DJ(u) || = " << DJu_Norm << endl;
         if (visualization) {
            sout << "solution\n" << mesh << u_next;
         }
      }
      // Updating iteration
      u_prev = u_next;
      J_prev = J_next;
      iter_counter++;
   } while (DJu_Norm > tol && iter_counter <= 2000);

   // 12. Save the final solution. This output can be viewed later
   //     using GLVis: "glvis -m refined.mesh -g sol_final.gf".
   ofstream sol_ofs2("sol_final.gf");
   sol_ofs2.precision(8);
   u_prev.Save(sol_ofs2);
   
   // 13. Free the used memory.
   delete fec;

   return 0;
}

double J(FiniteElementSpace *fes, GridFunction &u, double eps)
{
   // Computing J(u)
   GridFunction one_gf(fes);
   one_gf = 1.0;
   BilinearForm GradNorm(fes);
   LinearForm NonlinNorm(fes);
   ConstantCoefficient one(1.0);
   GradNorm.AddDomainIntegrator(new DiffusionIntegrator(one));
   GradNorm.Assemble();

   GridFunctionCoefficient u_gfc(&u);
   PowerCoefficient u4(u_gfc, 4);
   PowerCoefficient u2(u_gfc, 2);
   ProductCoefficient u2m(-2, u2);
   NonlinNorm.AddDomainIntegrator(new DomainLFIntegrator(u4));
   NonlinNorm.AddDomainIntegrator(new DomainLFIntegrator(u2m));
   NonlinNorm.AddDomainIntegrator(new DomainLFIntegrator(one));
   NonlinNorm.Assemble();
   return eps/2 * GradNorm.InnerProduct(u, u) + 1.0/4.0 / eps * NonlinNorm(one_gf);
}
