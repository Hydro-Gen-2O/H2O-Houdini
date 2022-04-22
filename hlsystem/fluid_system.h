/*
  FLUIDS v.1 - SPH Fluid Simulator for CPU and GPU
  Copyright (C) 2008. Rama Hoetzlein, http://www.rchoetzlein.com

  ZLib license
  This software is provided 'as-is', without any express or implied
  warranty.  In no event will the authors be held liable for any damages
  arising from the use of this software.

  Permission is granted to anyone to use this software for any purpose,
  including commercial applications, and to alter it and redistribute it
  freely, subject to the following restrictions:

  1. The origin of this software must not be misrepresented; you must not
     claim that you wrote the original software. If you use this software
     in a product, an acknowledgment in the product documentation would be
     appreciated but is not required.
  2. Altered source versions must be plainly marked as such, and must not be
     misrepresented as being the original software.
  3. This notice may not be removed or altered from any source distribution.
*/

#ifndef DEF_FLUID_SYS
	#define DEF_FLUID_SYS

	#include <vector>
	#include "fluid.h"
	#include <iostream>
	
	// Physical constants
	#define GRAVITY glm::dvec3(0, 0, -9.8)
	#define GRAVITY_ON 1

	// Tunable(ish) parameters
	#define m_DT 0.0083
	#define SPH_RADIUS 0.1
	#define REST_DENSITY 6378.0
	#define MAX_NEIGHBOR 50
	#define RELAXATION 600.0


	// Vector params
	#define SPH_VOLMIN glm::ivec3(-4, -4, 0)
	#define SPH_VOLMAX glm::ivec3(4, 4, 10)
	#define SPH_INITMIN	glm::ivec3(-2, -2, 3)
	#define SPH_INITMAX	glm::ivec3(2, 2, 9)

	class FluidSystem {
	public:
		FluidSystem ();
		virtual void Run ();

		void SPH_CreateExample(int n, int nmax, glm::dvec3 start, glm::dvec3 end);
		void setParameters(int ite, double visc, double vor, double tensile);
		void cleanUp();
	private:
		glm::dvec3 scaledMin = glm::dvec3(SPH_VOLMIN) * SPH_RADIUS;
		glm::dvec3 scaledMax = glm::dvec3(SPH_VOLMAX) * SPH_RADIUS;

		// the axis between the bounds of the fluid https://en.wikipedia.org/wiki/Space_diagonal
		glm::ivec3 gridSpaceDiag = glm::ivec3((scaledMax - scaledMin) / SPH_RADIUS);

		int totalGridCells = gridSpaceDiag.x * gridSpaceDiag.y * gridSpaceDiag.z;

		void PredictPositions();
		void FindNeighbors();
		void ComputeDensity();
		void ComputeLambda();
		void ComputeCorrections();
		void ApplyCorrections();
		void Advance();

		double PolyKernel(double dist);
		void SpikyKernel(glm::dvec3 &r);

		glm::ivec3 GetGridPos(const glm::dvec3 &pos);
		// get index in grid space
		int GetGridIndex(const glm::ivec3 &gridPos);
	public:
		std::vector<std::unique_ptr<Fluid>> fluidPs;
	private:
		// grid maps indexSpace To vector of fluid there
		std::vector<std::vector<int>> grid;
		std::vector<std::vector<int>> neighbors;

		int myIteration;
		double viscConst;
		double vortConst;
		double kCorr;
	};
#endif
