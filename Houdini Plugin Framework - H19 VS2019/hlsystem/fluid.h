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

#ifndef DEF_FLUID
	#define DEF_FLUID

	#include <glm/glm.hpp>
	#include "common_defs.h"

	class Fluid {
	public:
		Fluid(const glm::dvec3 &pos, DWORD d) : 
			pos(pos), clr(d), 
			predictPos(glm::dvec3(0.0)), deltaPos(glm::dvec3(0.0)),
			vel(glm::dvec3(0.0)), 
			density(0.0), lambda(0.0), vorticity(glm::dvec3(0.0))
		{}
		glm::dvec3		predictPos;
		glm::dvec3		pos;			// Basic particle (must match Particle class)
		DWORD			clr;
		glm::dvec3		vel;

		double density;	
		double lambda;
		glm::dvec3 deltaPos;
		glm::dvec3 vorticity;
	};

#endif
