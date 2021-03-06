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

#include <glm/gtx/norm.hpp>
#include <glm/gtx/string_cast.hpp>

#include "fluid_system.h"

FluidSystem::FluidSystem() :
	myIteration(2),
	viscConst(0.01),
	vortConst(0.0003),
	kCorr(0.0001)
{}

double FluidSystem::PolyKernel(double dist) {
	if (dist > SPH_RADIUS || dist == 0) {
		return 0.0;
	}
	double c = SPH_RADIUS * SPH_RADIUS - dist * dist;
	double wPoly = 315.0 / (64.0 * 3.141592 * pow(SPH_RADIUS, 9));
	return c * c * c * wPoly;
}

void FluidSystem::SpikyKernel(glm::dvec3& r) {
	double dist = glm::length(r);
	if (dist > SPH_RADIUS || dist == 0) {
		r = glm::dvec3(0.0);
		return;
	}
	r *= -45.0 / (3.141592 * pow(SPH_RADIUS, 6));
	r *= (SPH_RADIUS - dist) * (SPH_RADIUS - dist) / dist;
}

void FluidSystem::setParameters(int ite, double visc, double vor, double tensile)
{
	myIteration = ite;
	viscConst = visc;
	vortConst = vor;
	kCorr = tensile;
}

void FluidSystem::cleanUp()
{
	if (fluidPs.size() > 0)
	{
		fluidPs.clear();
	}
	if (grid.size() > 0)
	{
		grid.clear();
	}
	if (neighbors.size() > 0)
	{
		neighbors.clear();
	}
}

// SET SPH_RADIUS BEFORE THIS!
void FluidSystem::SPH_CreateExample(std::vector<glm::dvec3> p) {
	cleanUp();

	scaledMin = glm::dvec3(SPH_VOLMIN) * SPH_RADIUS;
	scaledMax = glm::dvec3(SPH_VOLMAX) * SPH_RADIUS;
	gridSpaceDiag = glm::ivec3((scaledMax - scaledMin) / SPH_RADIUS);
	totalGridCells = gridSpaceDiag.x * gridSpaceDiag.y * gridSpaceDiag.z;

	for (auto& pt : p) {
		fluidPs.push_back(std::make_unique<Fluid>(pt * SPH_RADIUS));
	}

	grid.reserve(totalGridCells);
	for (int _ = 0; _ < totalGridCells; ++_) {
		std::vector<int>& gridIndices = std::vector<int>();
		gridIndices.reserve(MAX_NEIGHBOR);
		grid.push_back(gridIndices);
	}

	for (int i = 0; i < fluidPs.size(); ++i) {
		neighbors.push_back(std::vector<int>());
	}
}

glm::ivec3 FluidSystem::GetGridPos(const glm::dvec3& pos) {
	return glm::ivec3(pos / SPH_RADIUS - SPH_VOLMIN);
}

int FluidSystem::GetGridIndex(const glm::ivec3& gridPos) {
	return gridPos.z * gridSpaceDiag.y * gridSpaceDiag.x + gridPos.y * gridSpaceDiag.x + gridPos.x;
}

void FluidSystem::Run() {
	PredictPositions();
	FindNeighbors();
	for (int _ = 0; _ < myIteration; ++_) {
		ComputeDensity();
		ComputeLambda();
		ComputeCorrections();
		ApplyCorrections();
	}
	Advance();
}

void FluidSystem::PredictPositions() {
	for (int i = 0; i < fluidPs.size(); ++i) {
		std::unique_ptr<Fluid>& p = fluidPs.at(i);
		glm::dvec3 deltaVel = FORCE * m_DT; // a * dt = change in v

		// apply force to velocity (gravity)
		p->vel += (double)GRAVITY_ON * deltaVel;

		p->predictPos = p->pos + (p->vel * m_DT);


		// Perform collision detection and response
		if (p->predictPos.y < scaledMin.y) { p->vel.y = 0.0; p->predictPos.y = scaledMin.y + 0.001; }
		if (p->predictPos.y > scaledMax.y) { p->vel.y = 0.0; p->predictPos.y = scaledMax.y - 0.001; }

		if (p->predictPos.x < scaledMin.x) { p->vel.x = 0.0; p->predictPos.x = scaledMin.x + 0.001; }
		if (p->predictPos.x > scaledMax.x) { p->vel.x = 0.0; p->predictPos.x = scaledMax.x - 0.001; }

		if (p->predictPos.z < scaledMin.z) { p->vel.z = 0.0; p->predictPos.z = scaledMin.z + 0.001; }
		if (p->predictPos.z > scaledMax.z) { p->vel.z = 0.0; p->predictPos.z = scaledMax.z - 0.001; }
	}
}

void FluidSystem::FindNeighbors() {
	for (int i = 0; i < totalGridCells; ++i) {
		grid.at(i).clear();
	}

	//equivalinet of insertgrid / update grid finding the postns within the grid
	for (int i = 0; i < fluidPs.size(); ++i) {
		glm::ivec3 gridPos = GetGridPos(fluidPs.at(i)->predictPos);
		int gIndex = GetGridIndex(gridPos);

		// this if shouldn't be necessary?
		if (0 <= gIndex && gIndex < grid.size()) {
			grid.at(gIndex).push_back(i); // maybe set a limit? (see MAX_NEIGHBOR)
		}
	}

	// equiv. Finding the Neighbors
	for (int i = 0; i < fluidPs.size(); ++i) {
	neighbor_loop:
		neighbors.at(i).clear(); // clear neighbors from prev.

		std::unique_ptr<Fluid>& p = fluidPs.at(i);
		glm::ivec3 gridPos = GetGridPos(p->predictPos);

		int SEARCH_SIZE = 1;
		// 2x2 neighborhood.
		for (int x = -SEARCH_SIZE; x <= SEARCH_SIZE; x++) {
			for (int y = -SEARCH_SIZE; y <= SEARCH_SIZE; y++) {
				for (int z = -SEARCH_SIZE; z <= SEARCH_SIZE; z++) {
					glm::ivec3 n = gridPos + glm::ivec3(x, y, z);
					if (0 <= n.x && n.x < gridSpaceDiag.x &&
						0 <= n.y && n.y < gridSpaceDiag.y &&
						0 <= n.z && n.z < gridSpaceDiag.z) {
						int gIndex = GetGridIndex(n);
						for (int pIndex : grid.at(gIndex)) { // each 
							std::unique_ptr<Fluid>& pcurr = fluidPs.at(pIndex);
							double lenR = glm::length(p->predictPos - pcurr->predictPos);
							if (lenR <= SPH_RADIUS) {
								neighbors.at(i).push_back(pIndex);
								//if (neighbors.at(i).size() >= MAX_NEIGHBOR) { goto neighbor_loop; }
							}
						}
					}
				}
			}
		}
	}
}

void FluidSystem::ComputeDensity() {
	for (int i = 0; i < fluidPs.size(); ++i) {
		std::unique_ptr<Fluid>& p = fluidPs.at(i);
		p->density = 0.0;
		for (int j : neighbors.at(i)) { // for each neighbor
			p->density += PolyKernel(glm::length(p->predictPos - fluidPs.at(j)->predictPos));
		}
	}
}

void FluidSystem::ComputeLambda() {
	for (int i = 0; i < fluidPs.size(); ++i) {
		std::unique_ptr<Fluid>& p = fluidPs.at(i);
		double sumGradients = 0.0;
		glm::dvec3 pGrad = glm::dvec3(0.0);
		for (int j : neighbors.at(i)) { // for each neighbor
			std::unique_ptr<Fluid>& pcurr = fluidPs.at(j);

			// Spiky Kernel - modifies r by ref
			glm::dvec3 r = (p->predictPos - pcurr->predictPos);
			SpikyKernel(r);
			r /= REST_DENSITY;
			// End Spiky Kernel
			sumGradients += glm::length2(r);
			pGrad += r; // -= r; ?? - i think += b/c -45
		}
		sumGradients += glm::length2(pGrad);
		double constraint = p->density / REST_DENSITY - 1.0; // real scale constraint
		p->lambda = -constraint / (sumGradients + RELAXATION); // maybe + 500 or so
	}
}

void FluidSystem::ComputeCorrections() {
	double polyDen = PolyKernel(0.2 * SPH_RADIUS);
	for (int i = 0; i < fluidPs.size(); ++i) {
		std::unique_ptr<Fluid>& p = fluidPs.at(i);
		p->deltaPos = glm::dvec3(0.0);
		for (int j : neighbors.at(i)) { // for each neighbor
			std::unique_ptr<Fluid>& pcurr = fluidPs.at(j);
			//---------Calculate SCORR-----
			double frac = PolyKernel(glm::length(p->predictPos - pcurr->predictPos)) / polyDen;
			double sCorr = -kCorr * frac * frac * frac * frac;
			//------------End SCORR calculation-------

			// Spiky Kernel - modifies r by ref
			glm::dvec3 grad = (p->predictPos - pcurr->predictPos);
			SpikyKernel(grad);
			grad /= REST_DENSITY;
			// End Spiky Kernel

			//p->deltaPos += r * (p->lambda + pcurr->lambda
			p->deltaPos += grad * (p->lambda + pcurr->lambda + sCorr);
		}
	}
}

void FluidSystem::ApplyCorrections() {
	for (std::unique_ptr<Fluid>& p : fluidPs) {
		p->predictPos += p->deltaPos;
	}
}

void FluidSystem::Advance() {
	//update all velocities
	for (int i = 0; i < fluidPs.size(); ++i) {
		std::unique_ptr<Fluid>& p = fluidPs.at(i);
		//std::cout << "for particle: " << i << std::endl;
		//std::cout << "p pos: " << glm::to_string(p->pos) << std::endl;
		//std::cout << "p pred: " << glm::to_string(p->predictPos) << std::endl;
		//std::cout << "p delP: " << glm::to_string(p->deltaPos) << std::endl;
		//std::cout << "p vel: " << glm::to_string(p->vel) << std::endl;
		//std::cout << "p dens: " << p->density << std::endl;
		//std::cout << "p lamb: " << p->lambda << std::endl << std::endl;
		p->vel = (p->predictPos - p->pos) / m_DT;

		// vorticity confinement here?
		p->pos = p->predictPos;
	}

	// VORTICITY CONFINEMENT
	for (int i = 0; i < fluidPs.size(); ++i) {
		std::unique_ptr<Fluid>& p = fluidPs.at(i);

		glm::dvec3 omega = glm::dvec3(0.0f);
		glm::dvec3 eta = glm::dvec3(0.0f);
		for (int j : neighbors.at(i)) {
			std::unique_ptr<Fluid>& pcurr = fluidPs.at(j);
			glm::dvec3 grad = (p->predictPos - pcurr->predictPos);
			SpikyKernel(grad);

			eta += grad;
			omega += glm::cross((pcurr->vel - p->vel), grad); // eqn 15 in pbf
		}
		eta *= glm::length(omega);

		glm::dvec3 vortForce = glm::cross(glm::normalize(eta), omega) * vortConst; // eqn 16
		p->vel += vortForce * m_DT;
	}
	// END VORTICITY CONFINEMENT

	// VISCOSITY
	for (int i = 0; i < fluidPs.size(); ++i) {
		std::unique_ptr<Fluid>& p = fluidPs.at(i);
		glm::dvec3 acc(0.0, 0.0, 0.0);
		for (int j : neighbors.at(i)) {
			std::unique_ptr<Fluid>& pcurr = fluidPs.at(j);
			acc += (pcurr->vel - p->vel) * PolyKernel(glm::length(p->predictPos - pcurr->predictPos));
		}
		p->tmp = acc;
	}

	for (std::unique_ptr<Fluid>& p : fluidPs) {
		p->vel += viscConst * p->tmp * m_DT;
	}
	// END VISCOSITY
}